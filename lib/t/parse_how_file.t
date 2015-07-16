#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl parse_how_file.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl parse_how_file.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/12/20 15:05:54

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use lib ('..');

use Test::More qw( no_plan );
use Test::Deep;
use Data::Dumper;
use parse_how_file qw ( parse_how_file ) ;


#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $test_entry_file = 'test_entry.how';

open( my $fh, '<', $test_entry_file )
    or die "Cannot open $test_entry_file, $!";

my $test_entry;

{
    local $/;
    $test_entry = <$fh>;
}


my($antigen_seqlength, $pdb_code, $antigen_chid, $antigen_seq,
   $contacts, @antibody_chainids )
    = parse_how_file::_parse_entry($test_entry);

my $exp_seq = 'NTVAAYNLTWKSTNFKTILEWEPKPVNQVYTVQISTKSGDWKSKCFYTTDTECDLTDEIVKDVKQTYLARVFSYPAGNVEPLYENSPEFTPYLETNLGQPTIQSFEQVGTKVNVTVEDERTLVRRNNTFLSLRDVFGKDLIYTLYYWKSSSSGKKTAKTNTNEFLIDVDKGENYCFSVQAVIPSRTVNRKSTDSPVECMG';

my $exp_contacts = '............................................................................................................................................H.H.H.......HHLLLLLL..................H.H.H..HHHHH.HH.......';


is( $antigen_seqlength, 200, "Antigen seq length parsed ok" );
is( $pdb_code, '1JPS', "PDB code parsed ok" );
is( $antigen_chid, 'T', "chain id parsed ok" );
is( $antigen_seq, $exp_seq, "antigen seq parsed okay" );
is( $contacts, $exp_contacts, "contacts parsed okay" );
cmp_deeply( \@antibody_chainids, [ qw( H L )  ],
            "antibody chain ids parsed okay" );

my $test_record_file = 'test.how';

my @return_array = parse_how_file($test_record_file);

is( scalar @return_array, 12, "parse_how_file parses all records" );
