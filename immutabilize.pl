#!/usr/bin/perl -w
# immutabilize.pl --- make all Moose packages immutable
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 21 Oct 2013
# Version: 0.01

use warnings;
use strict;

# Get all .pm files
my @module = `find -name '*.pm'`;

foreach my $file (@module) {

    chomp $file;
    
    my @modified = ();
    
    open(my $fh, '<', $file) or die "Cannot open file $file, $!";

    my $data;
    
    {
        local $/;
        $data = <$fh>;
    }

    my @package
        = $data
            =~ m{ ( (?: (?: \n|\A) package .*? (?= \n package|\n 1; )
                      | \n 1; .* \z ) ) }gxms;
    
    my @new_module = ();
    
    foreach (@package) {
        my $new_package = $_;
        if ( $_ =~ m{ \n use \s Moose; }xms ) {
            $new_package .= "\n__PACKAGE__->meta->make_immutable;\n\n";
        }
        push(@new_module, $new_package);
    }

    open(my $OUT, '>', $file) or die "Cannot open to $file to write";
    
    print {$OUT} @new_module;
}


__END__

=head1 NAME

immutabilize.pl - Describe the usage of script briefly

=head1 SYNOPSIS

immutabilize.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for immutabilize.pl, 

=head1 AUTHOR

Tom Northey, E<lt>zcbtfo4@acrm18E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
