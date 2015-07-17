#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl WEKA.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl WEKA.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/02/13 17:02:59

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
use Test::Deep;
use ARFF;

BEGIN { use_ok( 'WEKA' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest 'testStandardizeArffFiles' => sub {
    my $testArff = "testStdize.arff";
    my $tWEKA = WEKA->new();
    
    ok(! -z $tWEKA->standardizeArffFiles($testArff),
       "standardizeArffFiles (passed one arff) returns non-empty file");
    
    my $scndArff = "testStdizeSecondary.arff";

    my ($primaryStdized, $secondaryStdized)
        = $tWEKA->standardizeArffFiles($testArff, $scndArff);
    
    ok(! -z $primaryStdized && ! -z $secondaryStdized,
       "standardizeArffFiles (passed two arffs) returns two non-empty files");

    # If we can standardise the output of standardizeArffFiles then they must
    # be arff files!
    ok($tWEKA->standardizeArffFiles($primaryStdized)
           && $tWEKA->standardizeArffFiles($secondaryStdized),
       "and both files are arff files");
};

subtest 'testNominalToBinary' => sub {
    my $testArff = "testNom2Bin.arff";
    my $weka     = WEKA->new();

    my $nomAttrColNum = 5;
    my $outArffStr = $weka->nominalToBinary($testArff, $nomAttrColNum);

    my $aN = "secondary_str";
    ok($outArffStr =~ /$aN=E.*?
                       $aN=H.*?
                       $aN=C.*?
                       $aN=EH/xms,
       "nominal attribute replaced by binary attributes");
};

subtest 'testTrain' => sub {
    my $trainArff = 'WEKA.train.arff';
    my $rf = new_ok("WEKA::randomForest", [trainArff => $trainArff,
                                           model => 'test.model',
                                           removeAttribute => 1]);
    
    my $gotCmdStr = $rf->buildTrainCmd();
    like($gotCmdStr,
         qr/.* -Xmx.* -cp .* -t .*\.arff -F ".* -R 1" -c .* -W .* -- -I .* -K .*/,
     "buildTrainCmd ok");
    
    can_ok($rf, 'train');
};

subtest 'testTest' => sub {
    # Using training arff for test arff
    my $testArff = "WEKA.train.arff";
    my $model = "test.model";
    my $rf = new_ok("WEKA::randomForest", [testArff => $testArff,
                                           model => $model]);

    my $gotCmdStr = $rf->buildTestCmd();

    like($gotCmdStr,
         qr/.* -Xmx.* -cp .* -T .*\.arff -l $model .*/, "buildTestCmd ok");
};
