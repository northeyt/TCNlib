#!/usr/bin/env perl
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found
use Test::Class::Moose::Load 'TestsFor';
use Test::Class::Moose::Runner;

Test::Class::Moose::Runner->new->runtests();
