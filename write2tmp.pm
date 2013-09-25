package write2tmp;

use Moose;
use Moose::Util::TypeConstraints;

use Carp;

use File::Temp;

# Subtypes

# Must be subtyped here to avoid problems with types module using write2tmp
#subtype 'TempDirectory',
#        as 'Str',
#    where { -d $_ },
#    message { "$_ is not a directory" };

# Attributes

# NOTE modified isa to avoid circular reference with types package
# need to replace these in future

has 'dir' => (
    is => 'rw',
    isa => 'Str', # Dir
    default => '/tmp',
);

has 'suffix' => (
    is => 'rw',
    isa => 'Str',
    default => '.tmp'
);

has 'data' => (
    is => 'rw',
    isa => 'ArrayRef[Str]',
    required => 1,
);

has 'file_name' => (
    is => 'ro',
    isa => 'Str', #FileReadable
    builder => '_write_file',
    lazy => 1,
);

sub _write_file {
    my $self = shift;

    my %arg = ( #DIR => $self->dir,
                SUFFIX => $self->suffix,
                UNLINK => 0 );

    my $tmp = File::Temp->new(%arg);

    croak "data array contains no data" if ! @{ $self->data };
    
    print $tmp @{ $self->data };

    my $fname = $tmp->filename;

    return $fname;
}


1;
__END__

=head1 NAME

write2tmp - Perl extension for blah blah blah

=head1 SYNOPSIS

   use write2tmp;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for write2tmp, 

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
