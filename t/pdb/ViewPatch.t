#!/usr/bin/env perl
#   Tom Northey <zcbtfo4@acrm18>     2016/03/31 15:37:59

#########################

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::ViewPatch' ); }

#########################

new_ok("pdb::ViewPatch");

subtest "patchDescriptions" => sub {
    my $tViewer = pdb::ViewPatch->new(patchDir => "t/pdb/test-patchdir/");
    my $tID = "1a2y:C:47";
    ok(exists $tViewer->patchDescriptionLookup->{$tID},
       "test patch ID is present in patchDescriptionLookup");
    my @expected = qw(C.47 C.44 C.45 C.46 C.48 C.49 C.50 C.52 C.59 C.61 C.70);
    my $got = $tViewer->patchDescriptionLookup->{$tID};
    cmp_deeply($got, \@expected);
};

subtest "pymolPipe" => sub {
    my $tViewer = pdb::ViewPatch->new();
    ok($tViewer->pymolPipe);
};

subtest "prepareObject" => sub {
    my $tViewer = pdb::ViewPatch->new();
    ok($tViewer->prepareChainObject("1a2y:C"));
}
