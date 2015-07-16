#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl getresol.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl getresol.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/11/25 16:43:40

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::getresol' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $class = 'pdb::getresol';
my $pdb_file = '/acrm/data/pdb/pdb4hou.ent';

my $getresol = new_ok( $class, [ 'pdb_file',  $pdb_file ] );

my $output = "crystal, 2.40A/22.70%\n";

cmp_deeply( [ pdb::getresol::_parse_output($output) ],
            [ 'crystal', '2.40', '22.70' ],
            "_parse_output works okay"  );

can_ok( $getresol, 'run' );

   

                     
