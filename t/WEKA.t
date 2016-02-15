#!/usr/bin/env perl
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
use TCNUtil::ARFF;

use File::Basename;
my ($testFile, $testDir, $suffix) = fileparse($0);
chdir($testDir);

BEGIN { use_ok( 'TCNUtil::WEKA' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

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

subtest 'parseTableFromOutput, CSV output' => sub {
    my $testObj   = WEKA->new(posLabel => 'I', negLabel => 'S', undefLabel => '?');
    no warnings 'qw';
    my $testLines = [qw(inst#,actual,predicted,error,prediction
                        1,1:I,1:I,,0.879,
                        2,1:I,2:S,+,0.638,
                        3,2:S,1:I,+,0.507,
                        4,2:S,2:S,,0.507,
                        5,1:?,1:I,,0.520)];
    my $table = $testObj->parseTableFromOutput($testLines, 'CSV');
    is($table->total,     4, "expected number of instances in parsed table");
    is($table->true_pos,  1, "expected number of true positives");
    is($table->true_neg,  1, "expected number of true negatvies");
    is($table->false_pos, 1, "expected number of false positives");
    is($table->false_neg, 1, "expected number of false negatives");
};

subtest 'parseTableFromOutput, default output' => sub {
    my $testObj   = WEKA->new(posLabel => 'I', negLabel => 'S', undefLabel => '?');
    no warnings 'qw';
    my $testLines = [('     1        1:I        1:I       0.896 (2dd8:S:437)',
                      '     2        2:S        2:S       0.758 (1rvf:1:201)',
                      '     3        1:I        2:S   +   0.514 (3grw:A:322)',
                      '     4        2:S        1:I   +   0.782 (2xtj:A:292)')];
    
    my $table = $testObj->parseTableFromOutput($testLines, 'DEF');
    is($table->total,     4, "expected number of instances in parsed table");
    is($table->true_pos,  1, "expected number of true positives");
    is($table->true_neg,  1, "expected number of true negatvies");
    is($table->false_pos, 1, "expected number of false positives");
    is($table->false_neg, 1, "expected number of false negatives");
};


subtest 'translating unlabelled class labels' => sub {
    my $testObj   = WEKA->new(posLabel   => 'I', negLabel => 'S',
                              undefLabel => '?', translateUndefLabelTo => 'posLabel');
    no warnings 'qw';
    my $testLines = [qw(inst#,actual,predicted,error,prediction
                        1,1:I,1:I,,0.879,
                        2,1:I,2:S,+,0.638,
                        3,2:S,1:I,+,0.507,
                        4,2:S,2:S,,0.507,
                        5,1:?,1:I,,0.520)];
    my $table = $testObj->parseTableFromOutput($testLines);
    is($table->total,     5, "expected number of instances in parsed table");
    is($table->true_pos,  2, "expected number of true positives");
    is($table->true_neg,  1, "expected number of true negatvies");
    is($table->false_pos, 1, "expected number of false positives");
    is($table->false_neg, 1, "expected number of false negatives");
};
