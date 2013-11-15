#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl automatc_patches.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl automatc_patches.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/16 17:35:35

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use lib ( '..' );

use pdb::pdb;

use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'automatic_patches' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


my %arg
    = ( radius => 8,
        patch_type => 'normal',
        pdb_code => '1djs',
    );

my $auto = automatic_patches->new(%arg);

my @summary = ();

my %summ_hash = ();

foreach my $patch ($auto->get_patches) {
    if ( ref $patch ne 'patch' ) {
        next;
    }
    else {
        $summ_hash{$patch->summary. "\n"} = 1;
    }    
}

my $exp_patch_file = 'automatic_patches.out';

open(my $fh, '<', $exp_patch_file)
    or die "Canot open file $exp_patch_file, $!\n"; 

my @exp_summary = <$fh>;

my %exp_hash = ();

foreach my $summ (@exp_summary) {
    $exp_hash{$summ} = 1;
}

cmp_deeply(\%summ_hash, \%exp_hash,
           "get_patches produces correct patch summaries");

# Test to see if tmp xmas file write works
$arg{pdb_code} = '1nox';
my $noxmas = automatic_patches->new(%arg);

# Directly access pdb_file and pdb code  attr for testing purposes
$noxmas->{pdb_code} = '1noxm'; # Will not be found in xmas dir
$noxmas->{pdb_file} = '1nox.pdb';

ok($noxmas->xmas_file, "xmas file created when not found in xmas dir");


print "Testing  BUILDARGS for when given a pdb object ...\n";

my $chain = chain->new( pdb_code => '1djs', chain_id => 'A',
                        pdb_file => '1djs.pdb', xmas_file => '1djs.xmas',
                    );

my $ap_from_pdb = new_ok('automatic_patches', [ pdb_object => $chain,
                                                radius => 8,
                                                patch_type => 'contact' ]
                                            );
