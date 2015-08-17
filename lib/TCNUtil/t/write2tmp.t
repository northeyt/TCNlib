#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl write2tmp.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl write2tmp.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/25 11:36:05

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use Test::More qw( no_plan );
use Test::Deep;

use lib ( '..' );
BEGIN { use_ok( 'TCNUtil::write2tmp' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my @data = qw ( Q W E R T Y );

my $tmp = write2tmp->new( suffix => '.dat',
                       data   => [@data], );

my $fname = $tmp->file_name;

is(ref [ ( values %{ write2tmp->Cache } ) ]->[0], 'File::Temp',
   'Cache catching File::Temp objects');


is( $tmp->retain_file, 1, "retain_file ok with no arg" );
is($tmp->retain_file( file_name => $tmp->file_name ), 1,
    "retain_file ok with file_name" );
is( $tmp->retain_file( all => 1 ), 1, "retain_file ok with 'all' arg" );

print "Temp file name: " . $tmp->file_name . "\n";

undef $tmp;

is(ref [ ( values %{ write2tmp->Cache } ) ]->[0], 'File::Temp',
   'Cache retains objects when write2tmp reference is undef');

write2tmp->Cache_Limit(5);

my $to_be_removed;

my %arg = ( data => \@data,
            suffix => '.test' );

for my $i ( 0 .. 5 ) {
    if ( ! $i ) {
        $to_be_removed = write2tmp->new(%arg);
        $to_be_removed->file_name();
    }
    else {
        my $tmp = write2tmp->new( data => [ qw ( some data ) ],
                                  suffix => '.test' );
        $tmp->file_name();   
    }
}

is( exists write2tmp->Cache->{$to_be_removed->file_name}, '',
    "Cache limit works okay" );
