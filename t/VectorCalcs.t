#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl rotate2pc.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl rotate2pc.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/10/17 17:06:41

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;
use Math::VectorReal qw(:all);
use Math::MatrixReal;

use Math::Trig;
use lib ( '..' );

use Test::More qw( no_plan );
use Test::Deep;

BEGIN { use_ok( 'TCNUtil::VectorCalcs' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


my @vector = ( vector(2, 2, 2), vector(-1,-1,-1) );

my @c_vector = TCNUtil::VectorCalcs::meancenter(@vector);

cmp_deeply( [ $c_vector[0]->x, $c_vector[0]->y, $c_vector[0]->z,
              $c_vector[1]->x, $c_vector[1]->y, $c_vector[1]->z, ],
            [ 1.5, 1.5, 1.5, -1.5, -1.5, -1.5 ],
            "meancenter works ok" );

is(TCNUtil::VectorCalcs::innerproduct(@vector), -1, "innerproduct works okay");

my $R = TCNUtil::VectorCalcs::RM_about_vector( Z, pi );

my $r_vector = X * $R;


# Treatment of element 2 of return array is due to very small fp number
# that differs very slightly on different bit-architectures
cmp_deeply( [ $r_vector->x,
              int (1*(10**19) * $r_vector->y ),
              $r_vector->z ],
            [ -1,        , -1224, 0 ],
            "RM_about_vector works okay" );

my @point = ( vector(1.5,0,0), vector(-1.5,0,0),
               vector(0,2,0), vector(0,-2,0),
               vector(0,0,1), vector(0,0,-1) );

my @pc = TCNUtil::VectorCalcs::get_eigenvectors(@point);

cmp_deeply( [ [ $pc[0]->x, $pc[0]->y, $pc[0]->z, ],
              [ $pc[1]->x, $pc[1]->y, $pc[1]->z, ] ],
            [ [ 0, 1, 0 ], [ 1, 0, 0 ] ],
            "PC1 and 2 correct" );

$R = TCNUtil::VectorCalcs::rotate2pc(@point);

@point = map { $_ * $R } @point;

@pc = TCNUtil::VectorCalcs::get_eigenvectors(@point);

cmp_deeply( [ [ $pc[0]->x, $pc[0]->y, $pc[0]->z ],
              [ $pc[1]->x, $pc[1]->y, $pc[1]->z ] ],
            [ [ 1, 0, 0],
              [ 0, 1, 0], ],
            "\nrotate2pc aligns pc1 to X and pc2 to Y" );

                  
