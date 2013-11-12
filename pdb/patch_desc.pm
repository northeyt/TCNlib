package patch_desc;

use Moose;
use Moose::Util::TypeConstraints;
use types;

use Math::Trig;
use Math::VectorReal qw(:all);
use Math::MatrixReal;
use Statistics::PCA;

use pdb::xmas2pdb;
use pdb::rotate2pc qw(:all);

use TCNPerlVars;
use GLOBAL qw( three2one_lc );
use Data::Dumper;
use local::error;

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


# Methods

sub patch_order {
    my $self = shift;

    # Default contact threshold - if one side of patch has < threshold
    # contacts then it is considered a surface side and a patch order for
    # that side will be returned
    my $threshold = 3;
    
    if (@_){
        croak "Args must be passed as a hash ref"
            if ref $_[0] ne 'HASH';

        my %arg = %{ $_[0] };

        if ( exists $arg{threshold} ){
            croak "Threshold must be a positive integer"
                unless abs( int $arg{threshold} ) eq $arg{threshold};

            $threshold = $arg{threshold};
        }
    };
    
    croak "This method requires a parent pdb or chain object"
        if ! $self->has_parent;

    croak "Parent object must have xmas_file attribute set"
        if ! $self->parent->has_xmas_file;

    my @surface_atoms = $self->_surface_atoms;

    # Hash surf atoms by serial

    my %all_atom = map { $_->serial => $_ } @{ $self->parent->atom_array };
    
    my %surf_atom = map { $_->serial => $_ } @surface_atoms;

    my %patch_atom = map { $_->serial => $all_atom{$_->serial} }
                               @{ $self->patch->atom_array };
    
    my %patch_surf_atom
        = map { $_->serial => $_ }
            grep ( defined $patch_atom{$_->serial}, @surface_atoms );
    
    my %nonpatch_atom
        = map { $all_atom{$_}->serial => $all_atom{$_} }
            grep( ! defined $patch_atom{$_}, keys %all_atom );

    # Set central atom to origin
    my $cent_atom = $self->patch->central_atom;
    
    foreach my $atom (values %all_atom) {
        foreach my $coord ( 'x', 'y', 'z') {
            $atom->$coord( $atom->$coord - $cent_atom->$coord ); 
        }
    }

    # Get rot matrix for transform  x and y ->  patch PC1 and PC2
    my $RM = rotate2pc::rotate2pc( map { vector($_->x, $_->y, $_->z) }
                            values %patch_surf_atom );

    # Transform all atoms
    foreach my $atom (values %all_atom) {
        my $rVect = vector($atom->x, $atom->y, $atom->z) * $RM;
        
        foreach my $coord ('x', 'y', 'z') {
            $atom->$coord($rVect->$coord);
        }
    }
    
    # Determine number of atomic contacts on either side of patch
    my($posz_contacts, $negz_contacts)
        = $self->_surface_sides( \%patch_atom, \%nonpatch_atom );

    if ( ref $posz_contacts eq 'local::error' ) {
        my $message
            = "Could not determine surface side: "
             . $posz_contacts->message;

        my $error = local::error->new( message => $message,
                                       type => 'surface_side',
                                       parent => $posz_contacts );

        return $error;
    }

    if (   scalar @{$posz_contacts} > $threshold
        && scalar @{$negz_contacts} > $threshold ) {
        # return an error
        my $message
            = "patch order: neither side of patch has contact number"
              . " less than threshold";

        my $error
            = local::error->new( message => $message,
                                 type => 'no_surface_sides',
                                 data => { positive_side =>
                                               scalar @{$posz_contacts},
                                           negative_side =>
                                               scalar @{ $negz_contacts },
                                           patch => $self->patch },
                             );

        return $error;
    }

    my @return = ();
    
    my @residue_order = $self->_residue_order;

    my @aacid_order
        = map { three2one_lc($_) }
            map { $self->patch->atom_array->[$_]->resName }
                map { $self->patch->resid_index->{$_}->{CA} } @residue_order;

    if ( @{$posz_contacts} < $threshold ) {
        push( @return, join( '', @aacid_order ) );
    }

    if ( @{$negz_contacts} < $threshold ) {
        my $rev_string
            =  shift (@aacid_order)
             . reverse ( join( '', @aacid_order ) );
                    
        push( @return, $rev_string );
    }
    return @return;
}

sub _surface_atoms {
    my $self = shift;

    my @surface_atoms = ();

    my $xmas2pdb
        = xmas2pdb->new(xmas_file => $self->parent->xmas_file,
                        radii_file => $TCNPerlVars::radii_file,
                        xmas2pdb_file => $TCNPerlVars::xmas2pdb,
                        form => 'monomer',);
    
    my $ret = $self->parent->read_ASA($xmas2pdb);

    my @no_ASA = ();

    # Capture atoms with no ASA defined from xmas file
    # (normally altLoc atoms)
    foreach my $err ( @{$ret} ) {
        push( @no_ASA, $err->data->{atom}->serial() );
    }
    
    foreach my $atom ( @{ $self->parent->atom_array } ) {

        my $serial = $atom->serial();
        
        next if grep { /^$serial$/} @no_ASA;
        
        croak "Monomer ASA is not defined for atom " . $atom->serial
            if ! $atom->has_ASAm;
        
        if ( $atom->ASAm > 0 ) {
            push(@surface_atoms, $atom);
        }
    }
    return ( @surface_atoms );
}

# Returns number of atoms contacting either side of patch
# 
sub _surface_sides {
    my ( $self, $p_atom_h, $np_atom_h ) = @_;

    my %limit = ( xmin => 0, xmax => 0, ymin => 0, ymax => 0 );

    # Get patch x and y limits
    foreach my $p_atom ( values %{ $p_atom_h } ) {
        my $radius = $p_atom->radius();

        foreach my $axis ('x', 'y') {
            if ( ($p_atom->$axis - $radius) < $limit{$axis.'min'} ) {
                $limit{$axis.'min'} = $p_atom->$axis;
            }
            elsif ( ($p_atom->$axis + $radius) > $limit{$axis.'max'} ) {
                $limit{$axis.'max'} = $p_atom->$axis;
            }
        }
    }

    my @zpos = ();
    my @zneg = ();
    
    # Only consider non-patch atoms within x and y ranges to form min set
    foreach my $np_atom ( values %{ $np_atom_h } ) {

        if ( ! $np_atom->has_radius  || ! defined $np_atom->radius ) {
            my $message
                =  "non-patch atom " . $np_atom->serial
                 . " has no radius set";

            my $error = local::error->new( message => $message,
                                           type => 'no_radius_attribute',
                                           data => { atom => $np_atom }, );
            return $error;
        }
        
        next unless   $np_atom->x < $limit{xmax}
                   && $np_atom->x > $limit{xmin}
                   && $np_atom->y < $limit{ymax}
                   && $np_atom->y > $limit{ymin};

        foreach my $p_atom ( values %{ $p_atom_h } ) {

            if ( ! $p_atom->has_radius || ! defined $np_atom->radius) {
                my $message
                    =  "patch atom " . $np_atom->serial
                        . " has no radius set";
                
                my $error
                    = local::error->new( message => $message,
                                         type => 'no_radius_attribute',
                                         data => { atom => $np_atom }, );
                return $error;
            }
            
            my $dist_thresh =  $np_atom->radius + $p_atom->radius;

            my $distance = sqrt (   ( $np_atom->x - $p_atom->x )**2
                                  + ( $np_atom->y - $p_atom->y )**2
                                  + ( $np_atom->z - $p_atom->z )**2 );
            
            if ($distance < $dist_thresh) {
                
                abs $np_atom->z == $np_atom->z ? push(@zpos, $np_atom) 
                                               : push(@zneg, $np_atom);
                last;
            }
        }
    }
    return ( \@zpos, \@zneg );
}

sub _residue_order {
    my $self = shift;
    
    my @calpha
        = map { $self->parent->atom_array->[$self->parent->resid_index->{$_}->{CA}] }
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
