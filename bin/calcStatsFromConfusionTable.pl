#!/usr/bin/env perl
use strict;
use warnings;
use TCNUtil::confusion_table;
use Carp;
use Getopt::Long;

@ARGV or die "Supply four numbers: TP TN FP FN\n";

my $table = confusion_table->new(item_class => 'instance');
my ($TP, $TN, $FP, $FN) = @ARGV;

map {my $obj = bless {}, 'instance';
     my $datum = datum->new(object => $obj, prediction => 1, value => 1);
     $table->add_datum($datum)} (0..$TP - 1);

map {my $obj = bless {}, 'instance';
     my $datum = datum->new(object => $obj, prediction => 0, value => 0);
     $table->add_datum($datum)} (0..$TN - 1);

map {my $obj = bless {}, 'instance';
     my $datum = datum->new(object => $obj, prediction => 1, value => 0);
     $table->add_datum($datum)} (0..$FP - 1);

map {my $obj = bless {}, 'instance';
     my $datum = datum->new(object => $obj, prediction => 0, value => 1);
     $table->add_datum($datum)} (0..$FN - 1);

$table->print_all();
