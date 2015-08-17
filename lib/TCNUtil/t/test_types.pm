package t::test_types;

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;

use Carp;

has 'evennum' => (
    is => 'rw',
    isa => 'enum'
);

has 'evenby4' => (
    is => 'rw',
    isa => 'divby4'
);

has 'test_IOAllFile' => (
    is => 'rw',
    isa => 'IOAllObject',
    coerce => 1
);

has 'test_IOAllTemp' => (
    is => 'rw',
    isa => 'IOAllObject',
    coerce => 1,
);

has 'test_FileReadable' => (
    is => 'rw',
    isa => 'FileReadable',
    coerce => 1,
);



__PACKAGE__->meta->make_immutable;


__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

t::test_types - Perl extension for blah blah blah

=head1 SYNOPSIS

   use t::test_types;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for t::test_types, 

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
package t::test_types;

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;

use Carp;

has 'evennum' => (
    is => 'rw',
    isa => 'enum'
);

has 'evenby4' => (
    is => 'rw',
    isa => 'divby4'
);

has 'test_IOAllFile' => (
    is => 'rw',
    isa => 'IOAllObject',
    coerce => 1
);

has 'test_IOAllTemp' => (
    is => 'rw',
    isa => 'IOAllObject',
    coerce => 1,
);

has 'test_FileReadable' => (
    is => 'rw',
    isa => 'FileReadable',
    coerce => 1,
);



__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

t::test_types - Perl extension for blah blah blah

=head1 SYNOPSIS

   use t::test_types;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for t::test_types, 

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

__PACKAGE__->meta->make_immutable;


__PACKAGE__->meta->make_immutable;

