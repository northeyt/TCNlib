#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl uniprot.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl uniprot.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/11/12 11:22:06

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'uniprot::uniprot' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $ac = 'P35762';
my $file_name  = $ac . '.txt';
my $class = 'uniprot::uniprot';

my @arg = ( file_name => $file_name );

my $test_obj = new_ok( $class => \@arg );

$test_obj->data_array();

print "Creating object without file name specified ...\n";

@arg = ( accession_code => $ac );

my $test_obj_2 = new_ok( $class => \@arg );

$test_obj_2->data_array();

cmp_deeply( [ @{ $test_obj->data_array() }[ 0 .. 3 ] ],
            [ @{ $test_obj_2->data_array()}[ 0 .. 3 ] ],
            "Data retrieval using LWP works okay" );

ok($test_obj->data_hash, "data_hash okay");

is( $test_obj->organism_NCBI_TaxID->[0], '10090',
    "organism_NCBI_TaxID okay" );

my $bad_code = 'B5RX19';

my $bad_obj = uniprot::uniprot->new( accession_code => $bad_code ); 

is(ref $bad_obj->organism_NCBI_TaxID, 'local::error',
    "error object set to 'organism_NCBI_TaxID' attribute with bad object");

