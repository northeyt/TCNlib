package pdb::pdbFunctions;

use strict;
use warnings;
use Carp;
use TryCatch;
use write2tmp;
use TCNPerlVars;
use Math::VectorReal;
use Scalar::Util qw(looks_like_number);

# This subroutine returns a filename of a pdb file from either a pdb or chain
# object. If the input is already a valid filename then the filename is returned
sub getPDBFile {
    my $input = shift;
    
    my $inputFile = "";
    
    # Is input a pdb or chain object?
    try {
        my $inputAtomArray = eval {$input->atom_array()};
        $inputFile = _atoms2tmp($inputAtomArray);
        
    }
    catch ($err) {
        # Is input a file path?
        try{
            if (-e $input) {
                $inputFile = $input;
            }
            else {
                croak "Input file $input does not exist!";
            }
        }
        catch ($err) {
            # Is input an array of atoms?
            if (_isAtomAref($input)) {
                $inputFile = _atoms2tmp($input);
            }
            else {
                croak "Could not get pdb file from input $input!";
            }
        }
    };
    return $inputFile;
}

sub _atoms2tmp {
    my $atomAref = shift;
    
    my @atomStrings = map {"$_"} @{$atomAref};
    my $w2t = write2tmp->new(data => [@atomStrings]);

    return $w2t->file_name();
}


sub _isAtomAref {
    my $Aref = shift;

    return 0 if ref $Aref ne 'ARRAY';
    
    my $check = 1;

    foreach my $ele (@{$Aref}) {
        $check = 0 if ref $ele ne 'atom';
    }

    return $check;
}


# This function takes a list of chain/pdb objects and returns a ref to an array
# of all atoms found within all the passed objects
sub generateAtomAref {
    my @atomContaining = @_;

    my @atomList = ();

    foreach my $obj (@atomContaining){
        try {
            push(@atomList, @{$obj->atom_array()});
        }
        catch {
            croak "Could not get an atom_array from object $obj";
        }
    }
    return \@atomList;
}

# This function returns a list of PDB codes for PDBs that containing antibody
# chains. The list is obtained by parsing a SACS xml file
sub parseSACSxml {
    my ($SACSxml) = @_;

    # Use latest SACS list if no file is supplied
    $SACSxml = $TCNPerlVars::SACSxml if ! $SACSxml;

    open(my $FH, "<", $SACSxml)
        or croak "Cannot open SACS xml file '$SACSxml'";

    my @pdbCodes = ();
    
    while (my $line = <$FH>) {
        my ($pdbCode) = $line =~ /antibody pdb="(\S+)"/g;
        if ($pdbCode) {
            push(@pdbCodes, $pdbCode);
        }
    }
    return @pdbCodes;
}

# This function will find the mean 3D point for the atoms passed to it
# Returns a 3-element array (x, y, z)
sub findAtomMean {
    my @atoms = @_;

    my ($xsum, $ysum, $zsum) = (0, 0, 0);
    
    map { $xsum += $_->x() } @atoms;
    map { $ysum += $_->y() } @atoms;
    map { $zsum += $_->z() } @atoms;

    my $n = scalar @atoms;
    
    return($xsum / $n,
           $ysum / $n,
           $zsum / $n);
}

# Given a rotation matrix and a ref to an array of atoms,
# this function performs the rotation upon the atoms
sub rotateAtoms {
    my $RM = shift;
    my $atomAref = shift;

    foreach my $atom (@{$atomAref}) {
        my $rVect = vector($atom->x, $atom->y, $atom->z) * $RM;
        
        foreach my $coord ('x', 'y', 'z') {
            $atom->$coord($rVect->$coord);
        }
    }
}
            
# This functions compares two resSeqs to order them, like inbuilt cmp
# -1 is returned if $rsA should before $rsB
#  1 is returned is $rsA should come after $rsB
sub compare_resSeqs {
    my $rsA = shift;
    my $rsB = shift;

    my $rS2OrdHref = shift;
    $rS2OrdHref = {} if ! ref $rS2OrdHref eq 'HASH';
    
    croak "Two resSeqs must be passed to compare_resSeqs!"
        if (! (defined $rsA && defined $rsB));

    # If a hash has been supplied that defines the order of resSeqs
    # (e.g. according to some ATOM fields occurance), and both resSeqs
    # appear in this hash, then use this ordering
    if (exists $rS2OrdHref->{$rsA} && $rS2OrdHref->{$rsB}) {
        return $rS2OrdHref->{$rsA} <=> $rS2OrdHref->{$rsB};
    }
    # If one does not exist in this hash, guess ordering
    else {
        if (looks_like_number($rsA) && looks_like_number($rsB)) {
            # Both are resSeqs are numbers and can be sorted with a simple <=>
            return $rsA <=> $rsB;
        }   
        else {
            my ($rsA_num) = $rsA =~ /(\d+)+/g;
            my ($rsB_num) = $rsB =~ /(\d+)+/g;
        
            my ($rsA_suffix) = $rsA =~ /([A-Z]+)$/g;
            my ($rsB_suffix) = $rsB =~ /([A-Z]+)$/g;

            # Set suffixes to empty strings if not initialized
            foreach my $suffref (\$rsA_suffix, \$rsB_suffix) {
                ${$suffref} = "" if ! ${$suffref};
            }
            
            # Sort first by resSeq number then suffix
            # If resSeqs are the same, sort by suffix. Else, sort by resSeq
            if ($rsA_num eq $rsB_num) {
                return $rsA_suffix cmp $rsB_suffix;
            }
            else {
                return $rsA_num <=> $rsB_num;
            }
        }
    }
}

# This functions compares two resids to order them, like inbuilt cmp
# -1 is returned if $riA should before $riB
#  1 is returned is $riA should come after $riB
sub compare_resids {
    my ($riA, $riB) = @_;

    # Attempt to split on "." to get chainID and resSeq
    my ($riAChainID, $rsA) = split(/\./, $riA);
    if ($riAChainID eq $riA) {
        # No "." to split on, so assume that first character is a chainID
        $riAChainID = substr($riA, 0, 1);
        $rsA = substr($riA, 1);
    }
    my ($riBChainID, $rsB) = split(/\./, $riB);
    if ($riBChainID eq $riB) {
        # No "." to split on, so assume that first character is a chainID
        $riBChainID = substr($riB, 0, 1);
        $rsB = substr($riB, 1);
    }

    # Sort first by chainID, then resSeq
    return($riAChainID cmp $riBChainID || compare_resSeqs($rsA, $rsB));
}

1;
