#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl types.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl types.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/25 13:55:11

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use lib ( '..' );
use t::test_types;
use IO::All;

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $test = t::test_types->new();

# check subtype of subtypes
ok($test->evennum(2), 'Even num okay');
ok($test->evenby4(20), 'Div by 4 oay' );

my $io = io('test_file.txt');

print ref $io;
my $ref = ref $io;
print "\n";
if ($ref =~ /^IO::All/) {
    print "Regex is working\n"; 
}

ok( $test->test_IOAllFile( io('test_file.txt') ), 'IO::All::File type ok');
ok( $test->test_IOAllTemp( io( '?' ) ), 'IO::All::Temp type ok'  );

my $io = io('test_file.txt');
ok( $test->test_IOAllFile($io), 'Pass io as a named var');

ok( $test->test_IOAllFile( 'test_file.txt' ),
    'Coerce IOAllObject from str');

ok( $test->test_FileReadable( $io ),
    'FileReadable from io ok' );

ok( $test->test_FileReadable( 'test_file.txt' ),
    'FileReadable coerce from string ok');

my @data = qw( A B C D F );

ok($test->test_FileReadable( [@data] ),
   'FileReadable coerce from ArrayRef ok' );
