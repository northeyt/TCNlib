#!/usr/bin/env perl
use strict;
use warnings;
use pdb::ViewPatch;
use TCNUtil::ARFF::FileParser;
use Getopt::Long;

my $arffFile;
GetOptions("a=s", \$arffFile);

@ARGV or Usage();
my $patchDir = shift @ARGV;

my $viewer = pdb::ViewPatch->new(patchDir => $patchDir);
if ($arffFile) {
    my $arff = ARFF::FileParser->new(file => $arffFile, idAttributeIndex => 0)->parse();
    $viewer->arff($arff);
}

print "Please enter a patchID: ";
while (my $patchID = <>) {
    chomp $patchID;
    if (! $patchID) {
        print "Trying patch 1a2y:C:100\n\n";
        $patchID = "1a2y:C:100";
    }
    eval {$viewer->view($patchID)};
    print $@->{message} . "\n" if $@;
    print "\nPlease enter another patchID: ";
}

sub Usage {
    print <<EOF;
$0 [-a arff-file] patches-dir

viewPatch.pl allows the visualisation of patches in pymol. Given a directory of
patch files, the user can enter a patch id in order to visualise that patches.
That patch will appear in pymol imposed on top of its parent chain.

Args:
    patches-dir : directory containing patch files with names following the
                  convention pdbCodeChainID.patches. E.g, 1a2yC.patches

Opts:
    -a : optional arff file containing information about your patches.
         Feature values of a patch will be displayed in the terminal
         when that patch is selected.
EOF
    exit(1);
}
