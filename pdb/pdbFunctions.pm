package pdb::pdbFunctions;

use strict;
use warnings;
use Carp;
use TryCatch;
use write2tmp;

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
            my @atomStrings = map {"$_"} @{$input->atom_array()};
            my $w2t = write2tmp->new(data => [@atomStrings]);
            $inputFile = $w2t->file_name();
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
            croak $err;
        }
    };
    return $inputFile;
}

1;
