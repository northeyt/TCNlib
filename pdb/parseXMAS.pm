package pdb::parseXMAS;

use strict; 
use warnings;

use Moose;
use Carp;

### ATTRIBUTES #################################################################
################################################################################

has 'xmasData' => (
    is => 'rw',
    isa => 'ArrayRef',
    predicate => 'has_xmasData',
);

has 'ppHbDonor2AcceptorHref' => (
    is => 'rw',
    isa => 'HashRef',
);

has 'SSbondHref' => (
    is => 'rw',
    isa => 'HashRef',
);

has 'resID2secStructHref' => (
    is => 'rw',
    isa => 'HashRef'
);

### METHODS ####################################################################
################################################################################

sub parseXMAS {
    my $self = shift;

    my %sec_str = ();
    my %ppHbondDonor2Acceptor  = ();
    my %SSbonds = ();
    
    croak "xmasData must be assigned!" if ! $self->has_xmasData;
    
    my $inatoms_datatype = 0;
    my $inmolecules = 0;
    my $indata = 0;
    my $type = 0;
    my $chain = 'not defined';
    my $inSS_datatype = 0;     
    my $inHb_datatype = 0; 
        
    foreach my $line (@{$self->xmasData}) {
        chomp $line;	        
        
        # Clear flags on leaving atoms data section
        if ($inatoms_datatype && $line =~ /^<\/data/i){
            $inatoms_datatype = 0;
            $inmolecules = 0;
        }
        
        if ($line =~ /^<data type=molecules>$/i){
            $inmolecules = 1;
        }
        
        #Check for entering data atoms section
        if ($line =~ /<data type=atoms>/i){
            $inatoms_datatype = 1;
        }
        
        # Handle ATOMS data section
        if ($inatoms_datatype && $inmolecules){   
            if ($line =~ /<chain>/i){ 	
                $line =~ /<chain>(.*)<\/chain>/i;
                $chain = $1;
                $chain =~ s/\s//g;
                $chain = "" if $chain eq "|SP|";
            }
            elsif ($line =~ /<type>(.*)<\/type>/i){				
                $type = 0;
                $type = 1 if $line =~ /atom/i;
            }
        }
        
        if ($type){
            my @fields = parseResidueLine($line);
            if (@fields) {
                my $resnum = $fields[1];
                my $second_str = $fields[7];
                my $resid = "$chain.$resnum";
                $sec_str{$resid} = checkSSLabel($second_str);
            }
        }
        parseSSbondData($line, \$inSS_datatype, \%SSbonds);
        parseppHbondData($line, \$inHb_datatype, \%ppHbondDonor2Acceptor);
    }
    $self->ppHbDonor2AcceptorHref(\%ppHbondDonor2Acceptor);
    $self->SSbondHref(\%SSbonds);
    $self->resID2secStructHref(\%sec_str);
}

sub parseppHbondData {
    my ($line, $inHb_datatype, $ppHbondDonor2AcceptorHref) = @_;

    # Clear flags on leaving ppHbonds data section
    if ($$inHb_datatype && $line =~ /^<\/data/i){	
        $$inHb_datatype = 0;
    }
    
    # Handle PPHBONDS data section
    if($$inHb_datatype){
        my ($donorAtomNum, $acceptorAtomNum,
            $donorResnam, $donorChain, $donorResnum,
            $donorAtomnam, $acceptorResnam, $acceptorChain,
            $acceptorResnum, $acceptorAtomnam) = split(' ', $line);
        
        $ppHbondDonor2AcceptorHref->{$donorAtomNum} = $acceptorAtomNum;
    }
    
    #Check for entering data pphbonds section
    if ($line =~ /<data type=pphbonds>/i){
        $$inHb_datatype = 1;
    }
}

sub parseSSbondData {
    my ($line, $inSS_datatype, $SSbondsHref) = @_;

    # Clear flags on leaving ssbond data section
    if ($$inSS_datatype && $line =~ /^<\/data/i){	
        $$inSS_datatype = 0;
    }
    
    # Handle SSBOND data section
    if($$inSS_datatype){
        my ($atomnum1, $ch1, $resnum1,
            $atomnum2, $ch2, $resnum2) = split(' ', $line);
        
        $resnum1 =~ s/\.//g;
        $resnum2 =~ s/\.//g;
        
        $SSbondsHref->{$atomnum1} = $atomnum2;
        $SSbondsHref->{$atomnum2} = $atomnum1;
    }

    #Check for entering data ssbond section
    if ($line =~ /<data type=ssbond>/i){
        $$inSS_datatype = 1;
    }
}
    
sub parseResidueLine {
    my $line = shift;
    
    #parsing accessibility data from residue entry
    if ($line =~ /^<residue>\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+<\/residue>/i){
        my $resnam = $1;
        my $resnum = $2;
        my $molnum = $3;
        my $a1 = $4;
        my $a2 = $5;
        my $a3 = $6;
        my $a4 = $7;
        my $second_str = $8;
        
        $resnum =~ s/\.//g;     #removing dots
        return ($resnam, $resnum, $molnum, $a1, $a2, $a3, $a4, $second_str);
    }
    else {
        return ();
    }
}

# This function currently just translates the |SP| label to C. |SP| = space
# which indicates a lack of any assigned secondary structure, i.e. Coil 'C'.
sub checkSSLabel {
    my $label = shift;

    $label = 'C' if $label eq '|SP|';

    return $label;
}

1;
__END__

=head1 NAME

pdb::parseXMAS - provides functions to parse data from XMAS files

=head1 SYNOPSIS

   use pdb::parseXMAS;
   blah blah blah

=cut
