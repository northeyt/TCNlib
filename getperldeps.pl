#!/usr/bin/env perl
use Module::Build;
my @deps = qw(
Moose
MooseX::ClassAttribute
IO::CaptureOutput
Parallel::ForkManager
Math::VectorReal
Math::MatrixReal
Statistics::PCA
TryCatch
LWP
XML::DOM
XML::Simple
Test::MockObject::Extends
LWP::Protocol::https
);
# Currently installed with versions specified - this can be
# changed in the future if needs be
my $build = Module::Build->new(
    module_name => 'TCNlib',
    requires => {map {$_ => 0} @deps},
    );
$build->dispatch('installdeps');
