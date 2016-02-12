#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl get_files.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl get_files.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/14 15:09:06

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::get_files' ); }

#########################

new_ok('pdb::get_files');

my $test_depo_dir  = "test-depo-dir";
my $test_cache_dir = "test-cache-dir";
subtest 'get file from local depo' => sub {
    my $test_obj      = pdb::get_files->new(pdb_code => '1nox',
                                            pdbdir => $test_depo_dir);

    my $expected_file = "test-depo-dir/pdb1nox.ent";
    my $got_file      = $test_obj->pdb_file();

    is($got_file, $expected_file, "file found in local depo as expected");
};

subtest 'get pqs file' => sub {
    my $test_obj =  pdb::get_files->new(pdb_code => "1ndm",
                                        pqsdir => $test_depo_dir);
    ok($test_obj->pqs_file());
};

subtest 'get file from local cache' => sub {
    my $test_obj      = pdb::get_files->new(pdb_code => '1afv',
                                            pdbdir => $test_depo_dir);
    $test_obj->local_cache->cache_dir($test_cache_dir);

    my $expected_file = "test-cache-dir/pdb1afv.ent";
    my $got_file      = $test_obj->pdb_file();

    is($got_file, $expected_file, "file found in local cache as expected");
};
