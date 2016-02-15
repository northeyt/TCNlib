#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl ARFF.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl ARFF.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/03/20 11:20:30

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
use Test::Deep;
use Test::Exception;

BEGIN { use_ok(qw(TCNUtil::ARFF)); };
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

use File::Basename;
my ($testFile, $testDir, $suffix) = fileparse($0);
chdir($testDir);

use TCNUtil::ARFF::FileParser;
use TCNUtil::ARFF::AttributeDescription;
use TCNUtil::ARFF::ListOfAttributeDescriptions;
use TCNUtil::ARFF::Attribute;
use TCNUtil::ARFF::Instance;

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest 'testParse' => sub {
    my $tArff = getTestArff();

    my @expInstanceLines
        = ("1a2y:C:100,-0.197333,-0.814082,-0.590363,-0.347633,1.359971,-0.21583,I",
           "1a2y:C:101,-0.197333,-0.814082,1.69354,-0.347633,1.359971,-0.722922,I",
           "1adq:A:266,-0.197333,-0.814082,1.69354,-0.347633,-1.030077,?,S",
           "1adq:A:267,-0.197333,-0.814082,1.69354,-0.347633,-1.091358,?,S");
    
    my @expAttributeLines
        = map {'@attribute ' . $_} ("patchID string", "propensity numeric",
                                    "hydrophobicity numeric", "planarity numeric",
                                    "SSbonds numeric", "Hbonds numeric",
                                    "fosta_scorecons numeric", "intf_class {I,S}"
                                );
    
    cmp_deeply([$tArff->getInstanceStrings], \@expInstanceLines,
               "instances parsed from arff ok");
    cmp_deeply([$tArff->getAttributeDescriptionStrings], \@expAttributeLines,
               "attribute headers parsed from arff ok");
        
    my ($tInstance) = $tArff->allInstances;
    my %expAttributeHash = ("patchID"         => '1a2y:C:100',
                            "propensity"      => '-0.197333',
                            "hydrophobicity"  => '-0.814082',
                            "planarity"       => '-0.590363',
                            "SSbonds"         => '-0.347633',
                            "Hbonds"          => '1.359971',
                            "fosta_scorecons" => '-0.21583',
                            "intf_class"      => 'I');
    
    cmp_deeply({_hashAttributeName2Value($tInstance->getAttributes)},
               \%expAttributeHash,
               "instance->getAttributes returns attributes ok");
    
    is($tInstance->id,    "1a2y:C:100",
       "instance->id returns instance id value");
    is($tInstance->class, "I",
       "instance->class returns instance class value");
};

sub _hashAttributeName2Value {
    my @attributes = @_;
    return map {$_->name => $_->value} @attributes;
}

subtest 'testAddAttributeDescription' => sub {
    my $tArff = getTestArff();
    my $tAttrDesc = _getTestAttrDesc();
    
    $tArff->addAttributeDescription($tAttrDesc);

    my @expAttributeLines
        = map {'@attribute ' . $_}
            ("patchID string", "propensity numeric",
             "hydrophobicity numeric", "planarity numeric",
             "SSbonds numeric", "Hbonds numeric",
             "fosta_scorecons numeric", "intf_class {I,S}",
             "testAdd numeric");

    cmp_deeply([$tArff->getAttributeDescriptionStrings], \@expAttributeLines,
               "addAttributeFromDescription works ok");

};

sub test_addAttributesWithDescriptionToAllInstances {
    my $tArff = getTestArff();
    my $tAttrDesc = _getTestAttrDesc();
    lives_ok {tArff->addAttributesWithDescriptionToAllInstances($tAttrDesc)}
        "addAttributesWithDescriptionToAllInstances runs";
    
    ok($tArff->doAllInstancesHaveAttributeWithName($tAttrDesc->name()),
       "and all instances now have an attribute with the test description");
}


sub _getTestAttrDesc {
    my $tAttrName = "testAdd";
    return  ARFF::AttributeDescription->new(name => $tAttrName,
                                            type => 'numeric');
}

subtest 'testRemoveAttribute' => sub {
    my $tArff       = getTestArff();
    my $nameOfAttr2Remove = "propensity";

    my @expAttributeLines
        = map {'@attribute ' . $_}
            ("patchID string", "hydrophobicity numeric",
             "planarity numeric", "SSbonds numeric",
             "Hbonds numeric", "fosta_scorecons numeric",
             "intf_class {I,S}");
    
    $tArff->removeAttributeWithName($nameOfAttr2Remove);
    
    cmp_deeply([$tArff->getAttributeDescriptionStrings], \@expAttributeLines,
               "removeAttribute removes attribute from attribute lines ...");
    
    ok($tArff->doAllInstancesLackAttributeWithName($nameOfAttr2Remove),
       "... and test instances no longer have removed attribute");
};

subtest 'testAddAttributeAndValuesToInstances' => sub {
    my $tArff = getTestArff();

    my $newAttrName = "testAttr";
    my $newAttrDesc = ARFF::AttributeDescription->new(name => "testAttr");

    my @orderedValues = qw(1 2 3 4);

    $tArff->addAttributeValuesToInstances($newAttrDesc, \@orderedValues);
    cmp_deeply([$tArff->getValuesForAttributeWithName($newAttrName)],
               \@orderedValues,
               "values added correctly from ordered array");
    
    my %id2value = ("1a2y:C:100" => 1, "1a2y:C:101" => 2, "1adq:A:266" => 3,
                    "1adq:A:267" => 4);

    $tArff = getTestArff();

    $tArff->addAttributeValuesToInstances($newAttrDesc, \%id2value);
    
    cmp_deeply([$tArff->getValuesForAttributeWithName($newAttrName)],
               \@orderedValues,
               "values added correctly from id => value hash");
};

subtest 'testTransAttributeWithNameToBinaryFromNominal' => sub  {
    my $tArff = nomArff();
    my $nomAttrName  = "secondary_str";
    my $nomAttrType =  "{E,H,C,EH}";

    my @binAttrNames = map {"secondary_str=" . $_} qw(E H C EH);
    my @binAttrNameAndTypes = map {"$_ numeric"} @binAttrNames; 
    $tArff->transAttributeWithNameToBinaryFromNominal($nomAttrName);

    my $tFlag = 1;
    foreach my $binAttrName (@binAttrNames) {
        $tFlag = 0
            if ! $tArff->doAllInstancesHaveAttributeWithName($binAttrName);
    }
    ok($tFlag, "All instances have new binary attributes ...");

    my @expAttributeLines
        = map {'@attribute ' . $_}
            ("patchID string", "propensity numeric", "hydrophobicity numeric",
             "planarity numeric", @binAttrNameAndTypes, "SSbonds numeric",
             "Hbonds numeric", "fosta_scorecons numeric",
             "blast_scorecons numeric", "intf_class {S,I}");

    cmp_deeply([$tArff->getAttributeDescriptionStrings], \@expAttributeLines,
               "... and original nominal attribute has been removed");

    my $tInst = $tArff->getInstance(1);
    my @expValues = qw(1 0 0 0);
    my @gotValues = map {$tInst->getValueForAttributeWithName($_)} @binAttrNames;
    cmp_deeply(\@gotValues, \@expValues,
               "binary attribute values match expected for test instance");
};

subtest 'testStandardizeAgainstOtherArff' => sub {
    my $tArff
        = ARFF::FileParser->new(file => "testStdize.arff",
                                idAttributeIndex => 'first')->parse();
    my $arff2StdAgainst
        = ARFF::FileParser->new(file => "testStdizeSecondary.arff",
                                idAttributeIndex => 'first')->parse();
    
    my $expInstanceStr
        = "1a14:N:100,-0.836009,0.25036,-0.609557,-1.287254,"
        . "-3280241.930748,159168821.453041,-0.758384,0.496494,-3.352247,?,"
        . "0.661901,S";

    testStandardization($tArff, $expInstanceStr, $arff2StdAgainst);
};

### SUBROUTINES ################################################################
################################################################################

sub testStandardization {
    my $tArff = shift;
    my $firstInstanceExpStr = shift;
    my $arff2StdAgainst = shift;

    my @expInstanceIDs = map {$_->id} $tArff->allInstances();
    $tArff->standardize(defined $arff2StdAgainst ? $arff2StdAgainst : ());

    ok($tArff->allInstances(), "arff has instances after standardizaion");
    
    my @gotInstanceIDs = map {$_->id} $tArff->allInstances();

    cmp_deeply(\@gotInstanceIDs, \@expInstanceIDs,
               "instance IDs match pre and post standardization");
    
    is([$tArff->getInstanceStrings()]->[0],
       $firstInstanceExpStr,
       "attribute values of first instance standardized as expected");
}

sub getTestArff {
    return ARFF::FileParser->new(file => "test.arff", idAttributeIndex => 'first',
                                 classAttributeIndex => 'last')->parse();
}

sub nomArff {
    return ARFF::FileParser->new(file => "testNom2Bin.arff", idAttributeIndex => 'first',
                                 classAttributeIndex => 'last')->parse();
}

sub _allInstancesHaveAttr {
    my $tArff = shift;
    my $attr  = shift;
    my @instanceHasAttrChecks
        = map {exists $_->attributeHash->{$attr} ? 1 : 0}
            @{$tArff->instances};
    
    return grep {$_ eq 0} @instanceHasAttrChecks ? 0 : 1;
}
