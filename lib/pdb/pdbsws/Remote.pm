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
    if (! $responseString) {
        croak "pdbsws failed to return anything for query pdb $pdbCode, chain $chainID\n"
    }
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

1;
