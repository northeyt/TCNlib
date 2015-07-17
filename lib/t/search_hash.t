#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl search_hash.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl search_hash.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/23 09:44:51

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use Test::More qw( no_plan );
use Test::Deep;
use Test::Exception;

BEGIN { use_ok( 'search_hash' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


my %test_hash
    = ( A => { 1 => { Qwe => 0,
                      Wer => 1,
                      Ert => 2, },
               2 => { Rty => 3,
                      Tyu => 4,
                      Yui => 5, } },
        B => { 1 => { Uio => 6,
                      Opa => 7, },
               2 => { Asd => 8, } },
        C => { 2 => { Rty => 10,
                      Wer => 11,
                      Opa => 12 } }, );

print Dumper \%test_hash;

my %term_hash
    = ( A => { 2 => '' } );

print Dumper \%term_hash;

my $search = search_hash->new( hash        => \%test_hash,
                               depth       => 3,
                               term_hash => \%term_hash,
                           );

cmp_deeply( [ $search->search ],
            [ 3, 4, 5 ],
            "_search ok" );

%term_hash = ( '' => { 1 => '' } );

$search->term_hash(\%term_hash);

cmp_deeply( [$search->search],
            bag( 0, 1, 2, 6, 7 ),
            "branch from wildcard ok");

$test_hash{''} = 'some_val';

$search->hash(\%test_hash);

dies_ok( sub{ $search->search },
         'dies if wildcard val found in subj hash' );

$search->wildcard_value('%');

%term_hash = ( '%' => { 1 => '%' } );
$search->term_hash(\%term_hash);

cmp_deeply( [$search->search],
            bag( 0, 1, 2, 6, 7),
            'wildcard modified to non-subj hash key' );

%term_hash = ( '%' => { '%' => 'Rty' } );

$search->term_hash(\%term_hash);

cmp_deeply( [$search->search],
            bag(10, 3),
            'returns search term results from different branches' );

            
