package write2tmp;

use Moose;
use Moose::Util::TypeConstraints;
use Carp;

use File::Temp;

use MooseX::ClassAttribute;

class_has 'Cache' =>
    ( is => 'rw',
      isa => 'ArrayRef',
      default => sub { [] },
  );

no MooseX::ClassAttribute;

# Subtypes

# Attributes

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
    lazy => 1,
    builder => '_write_file',
);

# Methods

sub _write_file {
    my $self = shift;
    
    my %arg = ( DIR => $self->dir,
                SUFFIX => $self->suffix,
                UNLINK => 1);

    my $tmp = File::Temp->new(%arg);

    croak "data array contains no data" if ! @{ $self->data };
    
    print $tmp @{ $self->data };

    my $fname = $tmp->filename;

    push( @{ write2tmp->Cache }, $tmp );
    return $fname;
}



__PACKAGE__->meta->make_immutable;


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
