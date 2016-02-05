#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl idabchain.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl idabchain.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/03 17:19:46

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use lib ('../..');
use Test::More qw( no_plan );
use Test::Exception;
use Test::Deep;

use strict;
use warnings;
use pdb;

BEGIN { use_ok( 'pdb::idabchain' ); }

my $class = 'pdb::idabchain';

my $testObj = new_ok($class);

### Test _runExec
$testObj->input("notFilePath");
throws_ok{ $testObj->_runExec() } qr/Input file .* not exist/,
    'File does not exist thrown with bad input file';

my $expString = <<EOF;
Chain 1, A: Antigen
Chain 2, B: Antigen
EOF

$testObj->input("1djs.pdb");
is($testObj->_runExec(), $expString, "_runExec works ok");

$testObj->input(createDummyFile());
# Commented for now because it idabchain sends nasty output
# directly to terminal

#throws_ok{ $testObj->_runExec() } qr/idabchain run failed/,
#    'Exception thrown if idabchain fails';

$testObj->input(getAntibodyPDB());

my $expHash = { A => 'Antigen',
                B => 'Antigen',
                L => 'Light',
                H => 'Heavy',
                M => 'Light',
                K => 'Heavy' };

cmp_deeply($testObj->_parseOutput($testObj->_runExec()), $expHash,
           "_parseOutput works oK");

throws_ok{$testObj->chainIs('Antigen')} qr/no chain_id!/,
    "isChain throws error if input is not a chain";

$testObj->input(getChain());
is($testObj->chainIs(), 'Antigen', "chainIs works okay");

$testObj->input(getscFvChain());
is($testObj->chainIs(), 'scFv', "chainIs works okay with scFv");
cleanUp();

### Subroutines ################################################################

sub getChain {
    return chain->new(pdb_file => getAntibodyPDB(), chain_id => 'A');
}

sub getscFvChain{
    return chain->new(pdb_file => getscFvPDB(), chain_id => 'A');
}

sub getAntibodyPDB{
    my $fName = "1afv.pdb";
    return $fName;
}

sub getscFvPDB{
    my $fName = "1qok.pdb";
    return $fName;
}

sub createDummyFile {
    my $fName = 'dummyFile.tmp';
    open(my $fh, ">", $fName)
        or die "Could not create dummy file $fName!";

    print {$fh} "Some data\n";
    close $fh;
    
    return $fName;
}

sub cleanUp{
    my $dummyFile = createDummyFile();
    `rm $dummyFile`;
}


#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


