package pdb::pdbsws::querier;
use Moose::Role;

requires 'getACsFromPDBCodeAndChainID';
requires 'getIDsFromPDBCodeAndChainID';
requires 'mapResSeq2SwissProtNum';

package pdb::pdbsws::Remote;
use Moose;
use Carp;
use LWP::UserAgent;
use UNIPROT;

with 'pdb::pdbsws::querier';

has '_ua' => (
    is       => 'rw',
    isa      => 'LWP::UserAgent',
    required => 1,
    lazy     => 1,
    default  => sub {LWP::UserAgent->new()},
);

has '_URLBase' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
    lazy => 1,
    default => "http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=",
);

sub getACsFromPDBCodeAndChainID {
    my $self = shift;
    my ($pdbCode, $chainID) = @_;
    my $queryString = "pdb&id=$pdbCode&chain=$chainID";
    my $responseString = $self->_query($queryString);
    return map {$_->{AC}}
        $self->_getEntryHashesFromResponseString($responseString)
}

sub getIDsFromPDBCodeAndChainID {
    my $self = shift;
    my ($pdbCode, $chainID) = @_;
    my $queryString = "pdb&id=$pdbCode&chain=$chainID";
    my $responseString = $self->_query($queryString);
    return map {$_->{ID}}
        $self->_getEntryHashesFromResponseString($responseString)
}

sub mapResSeq2SwissProtNum {
    my $self     = shift;
    my $pdbCode  = shift;
    my $chainID  = shift;
    my $targetAC = shift;

    my $queryString    = "pdb&id=$pdbCode&chain=$chainID&all=yes";
    my $responseString = $self->_query($queryString);

    my @residueMappings = _parseMapFromResponseString($responseString);
    return map {$_->[0] => $_->[2]} grep {$_->[1] eq $targetAC} @residueMappings;
}

sub _parseMapFromResponseString {
    my $string = shift;
    _removeHeader(\$string);
    return map { [_parseMapLine($_)] } split(/\n/, $string);
}

sub _parseMapLine {
    my $line = shift;
    my ($pdbCode, $chainID, $chainSeq, $resName, $resSeq, $ac, $r1lc, $spNum)
        = split (/\s+/, $line);
    return ($resSeq, $ac, $spNum);
}

sub _removeHeader {
    my $stringRef = shift;
    ${$stringRef} =~ s/-.*-\n//xms;
}

sub _query {
    my $self         = shift;
    my $queryString  = shift;
    my $url          = $self->_URLBase. "$queryString";
    my $response     = $self->_ua->get($url);
    confess "user agent get request was not successful.\nURL: $url"
        if ! $response->is_success();

    return $response->decoded_content();
}

sub _getEntryHashesFromResponseString {
    my $self   = shift;
    my $string = shift;
    my @entryStrings = $self->_getEntryStringsFromResponseString($string);
    return map { {$self->_getHashFromEntryString($_)} } @entryStrings;
};

sub _getEntryStringsFromResponseString {
    my $self = shift;
    my $string = shift;
    # Grep to avoid trailing newline at end of string
    return grep {! /^\n$/} split(m{//}, $string);
}

sub _getHashFromEntryString {
    my $self = shift;
    my $string = shift;
    return $string =~ /(\S+?): \s+ (\S+?)\n+/gxms;
}

package pdb::pdbsws::Local;

use Moose;
use types;
use DBI;
use TCNPerlVars;
use Carp;

with 'pdb::pdbsws::querier';

has 'dbhost' => (
    isa => 'Str',
    is => 'rw',
    default => $TCNPerlVars::pghost,
);

has 'dbname' => (
    isa => 'Str',
    is => 'rw',
    default => 'pdbsws',
);

has 'dbh' => (
    isa => 'Ref',
    is => 'ro',
    lazy => 1,
    builder => '_get_dbh', 
);

sub _get_dbh {
    my $self = shift;

    my $dbname = $self->dbname();
    my $dbhost = $self->dbhost();
    
    my $dbh = DBI->connect("dbi:Pg:dbname=$dbname;host=$dbhost");
    croak "Could not connect to pdbsws" if ! $dbh;

    return $dbh;
}

sub getIDsFromPDBCodeAndChainID {
    my $self     = shift;
    my $pdbCode  = shift;
    my $chainID  = shift;
    my @acs      = $self->getACsFromPDBCodeAndChainID($pdbCode, $chainID);
    return map { $self->_getSwissProtIDFromAC($_) } @acs;
}

sub getACsFromPDBCodeAndChainID {
    my $self     = shift;
    my $pdbCode  = shift;
    my $chain    = shift;
    
    my $sql = "SELECT ac
               FROM pdbsws
               WHERE pdb = '$pdbCode'
               AND chain = '$chain'
               AND valid = 't'
               AND aligned = 't'
               AND ac != 'SHORT'
               AND ac != 'DNA'
               AND ac != 'ERROR';";

    my $sth = $self->dbh->prepare($sql);

    my @ac = ();
    
    if ($sth->execute) {
        while ( my ($pdb_ac) = $sth->fetchrow_array ){
            push(@ac, $pdb_ac);
        }
    }
    return @ac;
}

sub get_ac {
    my $self = shift;
    croak "get_ac is no longer supported! please replace with getSwissProtACFromPDBCodeAndChainID"
        . " and note that pdb code and chain id are now passed separately"
            . " \ne.g. getSwissProtACFromPDBCodeAndChainID('4hou', 'A')";
}

sub _getSwissProtIDFromAC {
    my $self     = shift;
    my $sprot_ac = shift;
    
    my $sql = "SELECT i.id FROM idac i, acac a WHERE a.altac = '$sprot_ac' AND i.ac = a.ac;";
    my $sprot_id = $self->dbh->selectrow_array($sql);
    return $sprot_id;
}

sub mapResSeq2SwissProtNum {
    my $self     = shift; 
    my $pdbCode  = shift;
    my $chainID  = shift;
    my $targetAC = shift;
    
    my $sql = "SELECT resid, pdbaa, ac, swsaa, swscount
               FROM alignment
               WHERE pdb = '$pdbCode'
               AND chain = '$chainID';";
    
    my $pdbswssth = $self->dbh->prepare($sql);

    my %resSeq2SprotResNum = ();
    if($pdbswssth->execute){
        while (my ($pdbResSeq, $pdbRes, $sprotAC, $sprotRes, $sprotResNum)
                   = $pdbswssth->fetchrow_array){
            if ($sprotAC eq $targetAC) {
                $resSeq2SprotResNum{$pdbResSeq} = $sprotResNum;
            }
        }
    }
    return %resSeq2SprotResNum;
}

package pdb::pdbsws::Factory;
use Moose;

has 'remote' => (
    is       => 'rw',
    isa      => 'Bool',
    required => 1,
    default  => 0,
);

sub getpdbsws {
    my $self = shift;
    my @args = @_;

    if ($self->remote) {
        return pdb::pdbsws::Remote->new(@args);
    }
    else {
        return pdb::pdbsws::Local->new(@args);
    }
}


1;
__END__

=head1 NAME

pdb::pdbsws - Perl extension for access to pdbsws, including methods for
              common searches

=head1 SYNOPSIS

   use pdb::pdbsws;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::pdbsws, 

Blah blah blah.

=head2 EXPORT

None by default.

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Tom Northey, E<lt>zcbtfo4@acrm18E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
