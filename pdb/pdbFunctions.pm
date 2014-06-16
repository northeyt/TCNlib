package pdb::pdbFunctions;

use strict;
use warnings;
use Carp;
use TryCatch;
use write2tmp;
use TCNPerlVars;
use Math::VectorReal;

# This subroutine returns a filename of a pdb file from either a pdb or chain
# object. If the input is already a valid filename then the filename is returned
sub getPDBFile {
    my $input = shift;
    
    my $inputFile = "";
    
    # Is input a pdb or chain object?
    try {
        $inputFile = $input->pdb_file();
        # If pdb file has not been assigned, create tmp file
        if (! $inputFile) {
            $inputFile = _atoms2tmp($input->atom_array());
        }
        
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

1;
