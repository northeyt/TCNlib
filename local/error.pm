package local::error;

use Moose;
use Carp;

has 'message' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
);

has 'caller' => (
    isa => 'Str',
    is => 'ro',
    default => sub { caller(1) },
);

has 'type' => (
    isa => 'Str',
    is => 'rw',
    required => 1,
);

has 'data' => (
    is => 'rw',
);

has 'parent' => (
    is => 'rw',
    isa => 'local::error',
    predicate => 'has_parent',
);

__PACKAGE__->meta->make_immutable;

1;
__END__

=head1 NAME

local::error - Perl extension for blah blah blah

=head1 SYNOPSIS

   use local::error;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for local::error, 

Blah blah blah.

=head2 EXPORT

None by default.

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

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
