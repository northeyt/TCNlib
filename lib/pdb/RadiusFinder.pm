package pdb::RadiusFinder;
use Moose;
use Carp;
use TCNUtil::types;
use TCNPerlVars;

has 'radiiFile' => (
    is => 'rw',
    isa => 'FileReadable',
    default => $TCNPerlVars::radii_file,
    lazy => 1,
    required => 1,
);

has '_atomRadiusHref' => (
    is => 'rw',
    isa => 'HashRef',
    lazy => 1,
    required => 1,
    builder => '_buildAtomRadiusHref'
);

has '_defaultRadiusHref' => (
    is => 'rw',
    isa => 'HashRef',
    lazy => 1,
    required => 1,
    builder => '_buildDefaultRadiusHref'    
);

sub findRadiusOfAtomFromNames {
    my $self = shift;
    my $resName  = uc(shift);
    my $atomName = uc(shift);

    croak "No radius defined for atom $atomName, res $resName"
        if ! exists $self->_atomRadiusHref->{$resName}->{$atomName};
    
    return $self->_atomRadiusHref->{$resName}->{$atomName};
}

sub findRadiusOfAtomFromElement {
    my $self    = shift;
    my $element = shift;
    croak "Element is not defined!" if ! defined $element;
    croak "No radius defined for atom of element $element!"
        if ! exists $self->_defaultRadiusHref->{$element};
    return $self->_defaultRadiusHref->{$element};
}


sub _buildAtomRadiusHref {
    my $self = shift;

    my $inFile = $self->radiiFile();
    open(my $DATA, "<", $inFile) or die "Cannot open file $inFile, $!";

    my %atomRadius = ();

    my @data = <$DATA>;
    while (my $line = shift @data) {
        next if $line =~ /^#/;
        my ($resName, $numAtoms, $wholeSolvAcc, $sideChainSolvAcc)
            = _parseResidueLine($line);
        my @atomLines = _readNextXLines(\@data, $numAtoms);
        $atomRadius{$resName} = {map {_parseAtomLine($_)} @atomLines};
    }
    return \%atomRadius;
}

sub _parseAtomLine {
    my $line = shift;
    chomp $line;
    my ($atomName, $radius) = split(/\s+/, $line);
    $atomName =~ s/\.//g;
    return ($atomName, $radius);
}

sub _readNextXLines {
    my ($dataAref, $x) = @_;
    return map {shift @{$dataAref}} (1..$x);
}

sub _parseResidueLine {
    my $line = shift;
    chomp $line;
    return split(/\s+/, $line);
}

sub _buildDefaultRadiusHref {
    return {C  => 1.80, N  => 1.60, S  => 1.85, O  => 1.40, P  => 1.90,
            CA => 2.07, FE => 1.47, CU => 1.78, ZN => 1.39, MG => 1.73};
}

1;
