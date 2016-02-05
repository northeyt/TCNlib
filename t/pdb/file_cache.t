#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl file_cache.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl file_cache.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/14 14:35:51

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::file_cache' ); }

#########################

my $test_cache_dir = "test-cache-dir";
new_ok('pdb::file_cache', [cache_dir => $test_cache_dir]);

subtest 'Test pdb file found in cache' => sub {
    my $test_obj = pdb::file_cache->new(cache_dir => $test_cache_dir,
                                        pdb_code  => '1afv');
    my $expected_file = "test-cache-dir/pdb1afv.ent";

    my $got_file = $test_obj->get_file();
    is($got_file, $expected_file, "returned file patch matches expected");
    ok(-e $got_file, "and path leads to an existing file");
};

subtest 'Test pdb file added to cache from remote' => sub {
    my $test_obj = pdb::file_cache->new(cache_dir => $test_cache_dir,
                                        pdb_code  => '1nox');

    # Ensure that file to be added to cache is not present already
    my $expected_file = "test-cache-dir/pdb1nox.ent";
    unlink($expected_file);

    my $got_file = $test_obj->get_file();
    is($got_file, $expected_file, "returned file patch matches expected");
    ok(-e $got_file, "and path leads to an existing file");

    # Remove created file
    unlink($expected_file);
}
