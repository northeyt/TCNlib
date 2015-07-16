package pdb::pdbsws;

use Moose;
use types;
use local::error;
use DBI;
use TCNPerlVars;

use Carp;

# Subtypes

# Attributes

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



# Methods

sub _get_dbh {
    my $self = shift;

    my $dbname = $self->dbname();
    my $dbhost = $self->dbhost();
    
    my $dbh = DBI->connect("dbi:Pg:dbname=$dbname;host=$dbhost");
    croak "Could not connect to pdbsws" if ! $dbh;

    return $dbh;
}

sub get_ac {
    my $self = shift;
    croak "get_ac: no pdbid supplied" if ! @_;

    my $pdbid = shift;
    
    croak "get_ac: pdbid '$pdbid' is not valid"
        unless $pdbid =~ /^(\d\w{3})([A-Z])$/i;

    my $pdb   = lc $1;
    my $chain = uc $2;

    my $sql
        = "SELECT DISTINCT ac FROM alignment WHERE pdb='$pdb' AND chain='$chain'";
    my $dbh = $self->dbh();
    
    my $sth = $dbh->prepare($sql);

    my @ac = ();
    
    if ($sth->execute) {
        while ( my ($pdb_ac) = $sth->fetchrow_array ){
            push(@ac, $pdb_ac);
        }
    }

    if (! @ac) {
        my $message
            = "get_ac: no accession code returned for pdbid '$pdbid'";

        my $error = local::error->new( message => $message,
                                       type => 'no_ac_for_pdbid',
                                       data => { pdbid => $pdbid,
                                                 sql => $sql,
                                                 pdbsws_object => $self, },
                                   );

        return $error;
    }
    return @ac;
}

sub seqFromAC {
    my $self = shift;
    my $ac   = shift;
    
    my $sql = "SELECT sequence FROM sprot WHERE ac = '$ac';";
    my $sequence = $self->dbh->selectrow_array($sql);
    if ($sequence){
        return $sequence;
    }
    else {
        croak "seqFromAC: No sequence returned for ac $ac";
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
