package pdb::antigen;

use Moose::Role;
use Moose::Util::TypeConstraints;
use TCNUtil::types;

use Carp;

requires 'atom_index';

# Subtypes

subtype 'ArrayRefOfValidChars',
    as 'ArrayRef[ValidChar]';


coerce 'ArrayRefOfValidChars',
    from 'Str',
    via { [$_] };

subtype 'ValidResidueID',
    as 'Str',
    where { $_ =~ /^[A-Z]\d+$/i };

subtype 'ArrayRefOfResidues',
    as 'ArrayRef[ValidResidueID]';


# Attributes

has 'antigen_chain_ids' => (
    is => 'rw',
    isa => 'ArrayRefOfValidChars',
    coerce => 1,
    builder => '_attempt_copy',
    lazy => 1,
);

has 'epitope_residue_array' => (
    is => 'rw',
    isa => 'ArrayRefOfResidues',
    predicate => 'has_epitope_residues',
);


has 'epitope_atom_index' => (
    is => 'ro',
    isa => 'HashRef',
    lazy => 1,
    builder => '_build_epitope_atom_index',
);

# Methods

sub _attempt_copy {
    my $self = shift;

    if (ref $self eq 'chain' ) {
        $self->antigen_chain_ids($self->chain_id);
    }
    else {
        croak "You must set antigen_chain_ids";
    }
}


sub _build_epitope_atom_index {
    my $self = shift;
    
    croak "Attribute epitope_residue_array must be set before "
        . "epitope_atom_index can be built"
            if ! $self->has_epitope_residues;

    my @residues = @{ $self->epitope_residue_array };
    
    my %hash = ();
    
    foreach my $residue (@residues) {
        my($chain, $num) = $residue =~ /([A-Z])(\d+)/;

        croak "Cannot define an epitope residue on a non-antigen chain"
            if ! grep {/^$chain$/} @{ $self->antigen_chain_ids };
        
        # Grab hash from main atom index

        croak "resSeq $chain$num does not exist. pdb: " . $self->pdb_code
            if ! exists $self->atom_index->{$chain}->{$num};
 
        my $resSeq_h = $self->atom_index->{$chain}->{$num};
        
        if ( ! exists $hash{$chain} ) {
            $hash{$chain} = {};
        }
        $hash{$chain}->{$num} = $resSeq_h;
    }

    return \%hash;
}


1;
__END__

=head1 NAME

pdb::antigen - Perl extension for blah blah blah

=head1 SYNOPSIS

   use pdb::antigen;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::antigen, 

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
