use 5.006;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'tcnutils',
    license             => 'artistic_2',
    dist_author         => q{Tom Northey <northeyt@gmail.com>},
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => 0,
    },
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'Moose'             => 2.1604
        'Foo::Bar::Module'  => 5.0401,
    },
    add_to_cleanup     => [ 'TCNUtil-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
