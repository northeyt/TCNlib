package write2tmp;

use Moose;
use Moose::Util::TypeConstraints;
use Carp;

use File::Temp;

use MooseX::ClassAttribute;

class_has 'Cache' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
);

class_has '_Cache_Array' => (
    is => 'rw',
    isa => 'ArrayRef[write2tmp]',
    default => sub { [] },
);

class_has 'Cache_Limit' => (
    is => 'rw',
    isa => 'Int',
    default => 0,
);

class_has 'retainAll' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
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

has 'retain' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

# Methods

sub _write_file {
    my $self = shift;
    
    my %arg = ( DIR => $self->dir,
                SUFFIX => $self->suffix,
                UNLINK => ( write2tmp->retainAll || $self->retain ? 0 : 1 )
            );

    my $tmp = File::Temp->new(%arg);
                
    print $tmp @{ $self->data };
    $tmp->flush();
    
    my $fname = $tmp->filename;

    if ( exists write2tmp->Cache->{$fname} ) {
        croak "write2tmp has created a file with a file name that already "
            . "exists!";
    }
    else {
        write2tmp->Cache->{$fname} = $tmp;
        push( @{ write2tmp->_Cache_Array }, $tmp );
        if ( write2tmp->Cache_Limit
                 && @{ write2tmp->_Cache_Array } > write2tmp->Cache_Limit ) {
            until ( @{ write2tmp->_Cache_Array } == write2tmp->Cache_Limit) {
                my $fh = shift @{ write2tmp->_Cache_Array };
                delete write2tmp->Cache->{ $fh->filename };
            }
        }
    }
    return $fname;
}

sub retain_file {
    my $self = shift;

    my %arg = @_;
    
    if ( ! %arg ) {
        write2tmp->Cache->{ $self->file_name }->unlink_on_destroy( 0 );
        return 1;
    }
    elsif ( exists $arg{all} && $arg{all} ) {
        foreach my $fh ( values %{ write2tmp->Cache } ) {
            $fh->unlink_on_destroy( 0 );
        }
        return 1;
    }
    elsif ( exists $arg{file_name} ) {
        if ( exists write2tmp->Cache->{ $arg{file_name} } ){
            my $fh = write2tmp->Cache->{ $arg{file_name} };
            $fh->unlink_on_destroy( 0 );
            return 1;
        }
        else {
            croak "cannot retain file $arg{file_name}: file is not in Cache";
        }
    }
    else {
        croak "write2tmp: no file was specified to retain";
    }
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
