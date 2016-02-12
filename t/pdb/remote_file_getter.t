#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl remote_file_getter.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl remote_file_getter.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/14 14:17:41

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::remote_file_getter' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


my $test = new_ok("pdb::remote_file_getter");

$test->pdb_code("1qok");

subtest 'Test get_pdb_file_data' => sub {
    my $test_pdb_code = "1qok";
    my $test_obj      = pdb::remote_file_getter->new(pdb_code => $test_pdb_code);
    my $exp_data      = _get_exp_data();

    my $got_data = $test_obj->get_pdb_file_data();

    is($got_data, $exp_data, "get_pdb_file_data returns expected data");
};

sub _get_exp_data {
    my $file = "1qok.pdb";
    open(my $IN, "<", $file) or die "Cannot open file $file, $!";
    my $data;
    {
        local $/;
        $data = <$IN>;
    }
    return $data;
}
