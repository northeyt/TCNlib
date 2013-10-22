package makepatch;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use pdb::pdb;

use Carp;

# subtypes

subtype 'ValidPatchType',
    as 'Str',
    where { $_ =~ m{ \A (?: contact|normal ) \s* \z }xms },
    message { "$_ is not a valid type - must be 'contact' or 'normal' " };

# attributes

has 'makepatch_file' => (
    is  => 'rw',
    isa => 'FileExecutable',
    required => 1,
);

has 'patch_type' => (
    is  => 'rw',
    isa => 'ValidPatchType',
    required => 1,
);

has 'radius' => (
    is  => 'rw',
    isa => 'Int',
    required => 1,
);

has 'central_atom' => (
    is => 'rw',
    isa => 'atom',
    required => 1,
);

has 'pdb_file' => (
    is => 'rw',
    isa => 'FileReadable',
    coerce => 1,
    required => 1,
);

has 'pdb_code' => (
    is  => 'rw',
    isa => 'Str',
    default => '????'
);

has 'output' => (
    is => 'ro',
    isa => 'ArrayRef[Str]',
    lazy => 1,
    builder => '_run_makepatch',
);


# Minimum absolute accessibility an atom must have if it is to be considered
# surface
has 'surf_min' => (
    is => 'rw',
    isa => 'Num',
    lazy => 1,
    default => 0,
);

# methods

sub _run_makepatch {
    my $self = shift;

    my $makepatch = $self->makepatch_file();

    if ( -l $makepatch ) {
        $makepatch = readlink $makepatch;
    }
    
    my $pdb = $self->pdb_file();
    my $radius = $self->radius();
    my $patch_type = $self->patch_type eq 'contact' ? '-c' : '';
    my $atomname = $self->central_atom->name();
    my $resspec
        = $self->central_atom->chainID() . $self->central_atom->resSeq();
    my $surf_min = $self->surf_min;
    
    my $cmd = "$makepatch -s -r $radius -m $surf_min"
              . " $patch_type $resspec $atomname $pdb";

    my @output = `$cmd`;

    croak "makepatch produced no output\ncommand run: $cmd"
        if ! @output;

    return [ @output ];
}

# methods

__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

run_make_patch - Perl extension for blah blah blah

=head1 SYNOPSIS

   use run_make_patch;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for run_make_patch, 

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
