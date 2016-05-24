package pdb::patch_desc;

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;

use Math::Trig;
use Math::VectorReal qw(:all);
use Math::MatrixReal;
use Statistics::PCA;

use TCNUtil::VectorCalcs qw(rotate2pc);
use pdb::pdbFunctions;

use TCNPerlVars;
use TCNUtil::GLOBAL qw( three2one_lc );
use Data::Dumper;
use TCNUtil::local::error;

use Carp;

# SubTypes

subtype 'ValidParent',
    as 'Ref',
    where { ref $_ eq 'pdb' || ref $_ eq 'chain' },
    message { "parent is not pdb or chain object" };

# Attributes

has 'patch' => (
    is => 'rw',
    isa => 'patch',
    required => 1,
);

has 'parent' => (
    is => 'rw',
    isa => 'ValidParent',
    predicate => 'has_parent',
);

# Default contact threshold - if one side of patch has < threshold
# contacts then it is considered a surface side and a patch order for
# that side will be returned
has 'threshold' => (
    is => 'rw',
    isa => 'Int',
    default => 4,
);

# Holds hash of orginal atom co-ordinates, so that co-ordinates can be reset
# if error occurs, etc
has 'orig_coords' => (
    is => 'rw',
    isa => 'HashRef',
    builder => '_build_orig_coords',
    lazy => 1,
);

# Methods

sub patch_order {
    my $self = shift;

    # If patch only contains one, two, or three residues, calculate simple
    # patch order
    if (scalar keys %{ $self->patch->resid_index } < 4){
        return $self->_simplePatchOrder($self);
    }

    croak "This method requires a parent pdb or chain object"
        if ! $self->has_parent;

    # Ensure parent and patch atom radii are (required for determining contacts)
    map {$_->readAtomRadii() if ! $_->atomRadiiAreAssigned}
        ($self->patch, $self->parent);
    
    # Get surface, non-solvent atoms of chain   
    my ($surface_errors, $surface_atoms) = $self->_surface_atoms;
    
    # Return error if any patch atoms do not have ASA defined
    my %no_ASA_atom = _getNoASAAtomsHash(@{$surface_errors});
   
    # Hash all atoms of parent
    my %all_atom
        = map { $_->serial => $_ }
            grep { ! $_->is_solvent() } @{ $self->parent->atom_array };
    
    # Hash all atoms of patch
    my %patch_atom
        = map { $_->serial => $all_atom{$_->serial} }
            @{ $self->patch->atom_array };

    # Return error if a patch atom has no ASA value
    my $ASACheck = $self->_checkPatchASAsError(\%patch_atom, \%no_ASA_atom);
    return $ASACheck if $ASACheck;
    
    # Get patch surface atoms
    my %patch_surf_atom
        = map { $_->serial => $_ }
            grep ( defined $patch_atom{$_->serial}, @{ $surface_atoms } );

    # Get all non patch atoms
    my %nonpatch_atom
        = map { $all_atom{$_}->serial => $all_atom{$_} }
            grep( ! defined $patch_atom{$_}, keys %all_atom );

    # Transform all points so that patch centre lies at origin and PC1 and PC2
    # of patch lay on x and y axis
    $self->centrePatch([values %all_atom], [values %patch_surf_atom]);
       
    # Get number of positive and negative-z contacts
    # (i.e. number of contacts on either face of patch)
    my($posz_contacts, $negz_contacts)
        = $self->_surface_sides(\%patch_surf_atom, \%all_atom);

    # Return error if _surface_sides returned error
    my $surfaceSidesError = $self->_checkSurfaceSidesError($posz_contacts);
    return $surfaceSidesError if $surfaceSidesError;

    # Get order of resids
    my @residue_order = $self->_residue_order;

    my @aacid_order = $self->residTo1lc(@residue_order);

    # Stringify patch order
    my $string = join('', @aacid_order);
    my $rev_string = reversePatchOrder($string);

    # Determine patch orders strings that represent the surface-facing
    # sides of patch
    my @surfaceFacingPOrders
        = $self->_surfaceFacingPOrders($string, $rev_string,
                                       $posz_contacts, $negz_contacts);

    # Reset parent co-ordinates
    $self->reset_parent_coords();
    
    return @surfaceFacingPOrders;
}

# This subroutine determines small patch orders (< 4 residues)
sub _simplePatchOrder {
    my $self = shift;
     
    my @simple_residue_order = ();
    
    my $central_atom = $self->patch->central_atom;
    
    push( @simple_residue_order, $central_atom->resName() );
        
    foreach my $resid ( keys %{ $self->patch->resid_index() } ) {
        
        my %atom_names = %{ $self->patch->resid_index->{$resid} };
        
        my $p_atom = [values %atom_names]->[0];
        
        next if $p_atom->resid eq $central_atom->resid;
        
        push(@simple_residue_order, $p_atom->resName());
    }
    
    my @aacid_order = ();
    
    foreach my $resName (@simple_residue_order){
        my $onelc = eval { three2one_lc($resName) };
            $onelc = 'X' if ! $onelc;
        push(@aacid_order, $onelc);
    }
    my $string = join( '', @aacid_order );
    return $string;
}

sub centrePatch {
    my $self = shift;
    my ($allAtomAref, $pSurfAtomAref) = @_;
    
    $self->centralAtom2Origin($allAtomAref);

    # Get rot matrix for transforming x and y unit vectors to patch PC1 and PC2
    my $RM = rotate2pc( map { vector($_->x, $_->y, $_->z) }
                            @{$pSurfAtomAref} );

    # Transform all atoms using rotation matrix
    pdb::pdbFunctions::rotateAtoms($RM, $allAtomAref);

}

# Gets non-solvent surface atoms of parent pdb/chain
# Returns refs to error array and surface atom array
# e.g. ($errorAref, $surfaceAtomAref) =  $self->_surface_atoms() 
sub _surface_atoms {
    my $self = shift;

    my @surface_atoms = ();
    my @errors = ();

    # Run read_ASA if it has not been run and check for any errors
    if (! $self->parent->has_read_ASA) {
        my $ret = $self->parent->read_ASA();

        if (ref $ret eq 'local::error') {
            push(@errors, $ret);
        }
        elsif (ref $ret eq 'ARRAY'){
            foreach (@{$ret}) {
                if (ref $_ eq 'local::error'){
                    push(@errors, $_);
                }
            }
        }
    }
    
    foreach my $atom (@{$self->parent->atom_array}) {
        if (! $atom->has_ASAb()){
            my $message
                = "atom " . $atom->serial() . " does not have ASAb assigned";
            
            my $error = local::error->new(message => $message,
                                          type => 'atom_no_ASAb',
                                          data => {atom => $atom});
            push(@errors, $error);
            next;
        }
        
        if ($atom->ASAb > 0 && ! $atom->is_solvent()) {
            push(@surface_atoms, $atom);
        }
    }
    return (\@errors, \@surface_atoms);
}

# Returns number of atoms contacting either side of patch
sub _surface_sides {
    my ($self, $ps_atom_h, $all_atom_h) = @_;
    my @p_atom_serial = map { $_->serial  } values %{ $ps_atom_h };
    my %limit = (xmin => 0, xmax => 0, ymin => 0, ymax => 0);
    
    # Get patch x and y limits
    foreach my $ps_atom ( values %{ $ps_atom_h } ) {
        my $radius = $ps_atom->radius();
        foreach my $axis ('x', 'y') {
            if ( ($ps_atom->$axis - $radius) < $limit{$axis.'min'} ) {
                $limit{$axis.'min'} = $ps_atom->$axis;
            }
            elsif ( ($ps_atom->$axis + $radius) > $limit{$axis.'max'} ) {
                $limit{$axis.'max'} = $ps_atom->$axis;
            }
        }
    }

    my @zpos = ();
    my @zneg = ();
    
    # Only consider non-patch atoms within x and y ranges to form min set
    foreach my $atom ( values %{ $all_atom_h } ) {
        my $serial = $atom->serial;
        next if grep ( /^$serial$/, @p_atom_serial );
        
        if (! $atom->has_radius  || ! defined $atom->radius) {
            my $message
                =  "non-patch atom " . $atom->serial . " has no radius set";

            my $error = local::error->new( message => $message,
                                           type => 'no_radius_attribute',
                                           data => { atom => $atom }, );
            return $error;
        }
        
        next unless   $atom->x < $limit{xmax}
                   && $atom->x > $limit{xmin}
                   && $atom->y < $limit{ymax}
                   && $atom->y > $limit{ymin};

        foreach my $ps_atom ( values %{ $ps_atom_h } ) {
            if (! $ps_atom->has_radius || ! defined $atom->radius) {
                my $message
                    =  "patch atom " . $atom->serial . " has no radius set";
                
                my $error
                    = local::error->new( message => $message,
                                         type => 'no_radius_attribute',
                                         data => { atom => $atom }, );
                return $error;
            }
            
            my $dist_thresh =  $atom->radius + $ps_atom->radius;
            my $distance = sqrt (   ( $atom->x - $ps_atom->x )**2
                                  + ( $atom->y - $ps_atom->y )**2
                                  + ( $atom->z - $ps_atom->z )**2 );
            
            if ($distance < $dist_thresh) {
                $atom->z > $ps_atom->z ? push(@zpos, $atom) 
                    : push(@zneg, $atom);
            }
        }
    }
    return (\@zpos, \@zneg);
}

# This method returns resids of patch, ordered by angle around the xy plane
# (the central residue is always the first in the array)
# e.g. ($centralResidue, @orderedPeripheralResidues) = $self->_residue_order()
sub _residue_order {
    my $self = shift;
    
    my @calpha
        = map { $self->parent->resid_index->{$_}->{CA} }
            keys %{ $self->patch->resid_index };
    
    my %angle = ();

    my $c_atom_resid = $self->patch->central_atom->resid();

    foreach my $p_atom (@calpha) {

        next if $p_atom->resid eq $c_atom_resid;
        
        my $vector = vector($p_atom->x, $p_atom->y, $p_atom->z);
        
        my $angle = $self->_true_angle_from_x($vector);
       
        $angle{$p_atom->resid} = $angle;      
    }

    my @sorted = sort { $angle{$a} <=> $angle{$b} } keys %angle;

    return ( $c_atom_resid, @sorted );
}

# Sub to deal with calulating obtuse angles by calculating angle on x y plane
sub _true_angle_from_x {
    my $self = shift;
    my($v) = @_;
    
    croak "$v is not a vector object" if ref $v ne 'Math::VectorReal';

    # Set z component to zero and normalise
    my $v2d = vector($v->x, $v->y, 0 );
    
    # Avoid attempts to divide by zero
    return 0 if $v2d->length == 0;
 
    my $inner_prod = ( $v2d . X ) / $v2d->length * 1; # Length X = 1

    my $angle = acos ($inner_prod);
    
    if ($v2d->y < 0) {
            $angle = (2 * pi) - $angle;
        }
    return $angle;    
}

sub reset_parent_coords {
    my $self = shift;
    my $orig_coords = $self->orig_coords();

    foreach my $atom ( @{ $self->parent->atom_array() } ) {
        foreach my $coord ( qw ( x y z ) ) {
            $atom->$coord( $orig_coords->{$atom->serial}->{$coord} );
        }
    }
}

sub residTo1lc {
    my $self = shift;
    
    my @resids = @_;
    my @oneLetterCodes = ();
    
    foreach my $res (@resids) {
        
        my @atoms = values %{$self->patch->resid_index->{$res}};

        my $resName = $atoms[0]->resName();
        
        my $onelc = eval { three2one_lc($resName) };
        $onelc = 'X' if ! $onelc;
        push(@oneLetterCodes, $onelc);
    }
    return @oneLetterCodes;
}

# This subroutine "flips" a patch order string, so that the order of the
# peripheral residues is reversed
sub reversePatchOrder{
    my $string = shift;

    my @residues = split("", $string);
    
    my $centralResidue = shift @residues;

    my $revString = $centralResidue . reverse (join('', @residues));

    return $revString;
}

sub _surfaceFacingPOrders {
    my $self = shift;
    my($pOrder, $revpOrder, $posZConts, $negZConts) = @_;

    my $threshold = $self->threshold();
    my @return = ();
    
    if (abs (scalar @{$posZConts} - scalar @{$negZConts}) < $threshold) {
        # Difference between contacts on either side of patch is not enough
        # to choose one side over the other: return both
        push(@return, ($pOrder, $revpOrder));
        
    }elsif (scalar @{$posZConts} < scalar @{$negZConts}) {
        push(@return, $pOrder);
    }
    else {
        push(@return, $revpOrder);
    }
    
    return @return;
}


# Transforms all atoms passed in array so that patch central atom lies at
# origin (0,0,0)
sub centralAtom2Origin {
    my $self = shift;
    my $allAtomAref = shift;

    # Transform all atom co-ordinates so patch central atom is at origin
    my %cent_atom_coord = (x => $self->patch->central_atom->x(),
                           y => $self->patch->central_atom->y(),
                           z => $self->patch->central_atom->z(),);

    foreach my $atom (@{$allAtomAref}) {
        foreach my $coord ( 'x', 'y', 'z') {
            $atom->$coord( $atom->$coord - $cent_atom_coord{$coord} ); 
        }
    }
}

# Check return value of _surface_sides for errors. If error has been returned,
# returns error. Returns 0 if not.
sub _checkSurfaceSidesError {
    my $self = shift;
    my $surfaceSidesReturn = shift;
    
    if ( ref $surfaceSidesReturn eq 'local::error' ) {
        my $message
            = "Could not determine surface side: "
                . $surfaceSidesReturn->message;
        
        my $error = local::error->new(message => $message,
                                      type => 'surface_side',
                                      parent => $surfaceSidesReturn);
        
        $self->reset_parent_coords();
        return $error;
    }
    else {
        return 0;
    }
}

sub _checkPatchASAsError {
    my $self = shift;
    my($pAtomHref, $noASAatomHref) = @_;

    foreach my $serial (keys %{$pAtomHref}) {
        if (exists $noASAatomHref->{$serial}) {
            my $message
                = "patch desc: patch atom " . $serial . " has no ASA value";
            
            my $atom = $pAtomHref->{$serial};
            
            my $error = local::error->new( message => $message,
                                           type => 'patch_atom_no_ASA',
                                           data => { atom => $atom }, );

            $self->reset_parent_coords();
            return $error;
        }
    }
    return 0;
}

sub _build_orig_coords {
    my $self = shift;

    my %orig_coords = ();
    foreach my $atom ( @{ $self->parent->atom_array() } ) {
        $orig_coords{ $atom->serial() }
            = { x => $atom->x(), y => $atom->y(), z => $atom->z() }; 
    }
    return \%orig_coords;
}

# Given an array of errors returned by _surface_atoms, this sub returns a hash
# of form serial -> atom where atoms have no ASAb value
sub _getNoASAAtomsHash {
    my @surface_atoms_errors = @_;

    my %no_ASA_atom = ();
    
    foreach my $error (@surface_atoms_errors) {
        if ( $error->type() eq 'atom_no_ASAb' ) {
            my $atom = $error->data->{atom};
            $no_ASA_atom{ $atom->serial() } = $atom;
        }
    }
    return %no_ASA_atom;
}

__PACKAGE__->meta->make_immutable;

1;
__END__

=head1 NAME

pdb::patch_desc - object to generate descriptions for patches

=head1 SYNOPSIS

   use pdb::patch_desc;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::patch_desc, 

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
