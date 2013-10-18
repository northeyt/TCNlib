package patch_desc;

use Moose;
use Moose::Util::TypeConstraints;
use types;

with 'MooseX::Clone';

use Math::Trig;
use Math::VectorReal qw(:all);
use Math::MatrixReal;

use Statistics::PCA;
use pdb::xmas2pdb;
use TCNPerlVars;
use Data::Dumper;

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

    croak "This method requires a parent pdb or chain object"
        if ! $self->has_parent;

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

    print Dumper $all_atom{ $self->patch->central_atom->serial };

    # Move all atoms so that patch central atom is at origin
    $self->_centralatom2origin(\%all_atom);
    
    # Rotate all atoms so that x and y axis correspond to patch PC1 and PC2
    $self->_orientate2patch(values %patch_surf_atom);

    # Determine number of atomic contacts on either side of patch
    my($posz_contacts, $negz_contacts)
        = $self->_surface_sides( \%patch_atom, \%nonpatch_atom );
    
    my $contact_threshold = 35;
    if ( $posz_contacts < $contact_threshold )  {
        print "Positive z face:\n";
        $self->_residue_order(1);
    }

    if ( $negz_contacts < $contact_threshold ) {
        print "Negative z face:\n";
        foreach my $atom (values %patch_atom) {

            # Rotate pi rad around y then determine res order
            my $vector = [ $atom->x, $atom->y, $atom->z ];

            #if ($atom->name eq 'CA') {
            #    print $atom->resSeq . Dumper $vector;
            #}
            
            #my $newvector
            #    = _rotate( $vector, _rotation_matrix( pi , 'y' ) );
            
            #$atom->x( $newvector->[0] );
            #$atom->y( $newvector->[1] );
            #$atom->z( $newvector->[2] );
        }
        $self->_residue_order(-1);    
    }   
}

sub _surface_atoms {
    my $self = shift;

    my @surface_atoms = ();
    
    my $xmas2pdb
        = xmas2pdb->new(xmas_file => $self->parent->xmas_file,
                        radii_file => $TCNPerlVars::radii_file,
                        xmas2pdb_file => $TCNPerlVars::xmas2pdb,
                        form => 'monomer',);
    
    $self->parent->read_ASA($xmas2pdb);
    
    foreach my $atom ( @{ $self->parent->atom_array } ) {
        croak "Monomer ASA is not defined for atom " . $atom->serial
            if ! $atom->has_ASAm;
        
        if ( $atom->ASAm > 0 ) {
            push(@surface_atoms, $atom);
        }
    }
    return @surface_atoms;
}

# Transforms all atoms so that central atom of patch is at origin (0, 0, 0)
sub _centralatom2origin {
    my $self = shift;
    my $atom_h = shift;

    my $central_atom = $self->patch->central_atom;
    
    my($xi, $yi, $zi)
        = ( $central_atom->x(), $central_atom->y(), $central_atom->z() );
    
    # Subtract xi, yi and zi from all atoms
    foreach my $atom ( values %{ $atom_h } ) {
        $atom->x( $atom->x() - $xi );
        $atom->y( $atom->y() - $yi );
        $atom->z( $atom->z() - $zi );
    }
}


sub _orientate2patch_v3 {
    my $self = shift;
    my @patch_surf_atom = @_;

    my @pc = @{ $self->_get_eigenvectors(@patch_surf_atom) }[ 0 .. 2 ];
    
    # From vector objects from arrays
    @pc = map { vector( $_->[0], $_->[1], $_->[2] ) } @pc;
    
    my $rot_matrix = $self->_RMforvector2x($pc[0]);

    print "\nRotation Matrix:\n" , $rot_matrix ,"\n";
    
    # Rotate all points by matrices
    foreach my $point ( @{ $self->parent->atom_array } ) {
        my $vector = vector( $point->x(), $point->y(), $point->z() ) ;
        my $newvector = $vector * $rot_matrix;
            $point->x( $newvector->x );
            $point->y( $newvector->y );
            $point->z( $newvector->z );
    }
    
    # Also rotate pcs
    for my $j( 0 .. 2 ) {
        $pc[$j] = $pc[$j] * $rot_matrix; 
    }
    print "PCs after rotation: " . Dumper \@pc;

    #my $angle = atan( $pc[1]->y, $pc[1]->z );

    #$rot_matrix = $self->_rot_aboutx($angle);
    
    # pc2 will now be in z y plane, so rotate to reduce pc2->z to 0
    # Rotate all points by matrices
    #foreach my $point ( @{ $self->parent->atom_array } ) {
    #    my $vector = vector( $point->x(), $point->y(), $point->z() ) ;
    #    my $newvector = $vector * $rot_matrix;
    #    $point->x( $newvector->x );
    #    $point->y( $newvector->y );
    #    $point->z( $newvector->z );
    #}
    
}

sub _rot_aboutx {
    my $self = shift;
    my $angle = shift;

    my $R
        = Math::MatrixReal->new_from_rows(
            [ [ 1, 0          , 0,            ],
              [ 0, cos($angle), -sin($angle), ],
              [ 0, sin($angle),  cos($angle)  ],
          ],
        );

    return $R;
            
}


sub _new_orientate2patch {
    my $self = shift;
    my @patch_surf_atom = @_;

    croak "No patch surface atoms passed to orientate2patch"
        if ! @patch_surf_atom;

    my @pc = @{ $self->_get_eigenvectors(@patch_surf_atom) }[ 0 .. 2 ];

    # From vector objects from arrays
    @pc = map { vector( $_->[0], $_->[1], $_->[2] ) } @pc;

    my $rot_matrix = vector_matrix( @pc  );

    print "\nRotation Matrix:\n" , $rot_matrix ,"\n";
    
    # Rotate all points by matrices
    foreach my $point ( @{ $self->parent->atom_array } ) {
        my $vector = vector( $point->x(), $point->y(), $point->z() ) ;
        my $newvector = $vector * $rot_matrix;
            $point->x( $newvector->x );
            $point->y( $newvector->y );
            $point->z( $newvector->z );
    }
    
    # Also rotate pcs
    for my $j( 0 .. 2 ) {
        $pc[$j] = $pc[$j] * $rot_matrix; 
    }

    print "PCs after rotation: " . Dumper \@pc;
}

# Rotates all atoms so that x unit vector = Patch PCA vector 1
# and y unit vector = Patch PCA vector 2
sub _orientate2patch {
    my $self = shift;
    my @patch_surf_atom = @_;

    croak "No patch surface atoms passed to orientate2patch"
        if ! @patch_surf_atom;
    
    # Get eigenvectors for which axes are to be aligned to
    my @pc =  @{ $self->_get_eigenvectors(@patch_surf_atom) }[ 0 .. 1 ];
    
    # Rotations
    # Array form: pcomponent index, axis to rotate about, vector indices
    my @rot = ( [ 0, 'x', [ 2, 1 ] ],
                [ 0, 'z', [ 1, 0 ] ],
                [ 1, 'x', [ 2, 1 ] ],
            );
    
    for (my $i = 0 ; $i < @rot ; ++$i ) {

        my $axis = $rot[$i]->[1];
        
        print "\nVectors pre rotation about $axis\n"
            . Dumper \@pc;
    
        # Get rotation angle
        my $vector = $pc[ $rot[$i]->[0] ] ;
        
        my $angle = _true_rot_angle( $vector->[ $rot[$i]->[2]->[0] ],
                                     $vector->[ $rot[$i]->[2]->[1] ], );
    
        # Get rotation matrix
        my $matrix = _rotation_matrix( $angle, $rot[$i]->[1] );
    
        # Rotate all points by matrices
        foreach my $point ( @{ $self->parent->atom_array } ) {
            my $vector = [ $point->x(), $point->y(), $point->z() ];
            my $newvector = _rotate( $vector, $matrix );
            $point->x( $newvector->[0] );
            $point->y( $newvector->[1] );
            $point->z( $newvector->[2] );
        }
    
        # Also rotate pcs
        for my $j( 0 .. 1 ) {
            $pc[$j] = _rotate( $pc[$j], $matrix ); 
        }
    
        print "\nVectors after rot about $axis \n"
             . Dumper \@pc;
    }
}

sub _get_eigenvectors{
    my $self = shift;
    my @patch_surf_atom = @_;
    
    my $patch = $self->patch;
    my $parent = $self->parent;

    my $data = [ map { [ $_->x(), $_->y(), $_->z() ] } @patch_surf_atom  ];

    my $pca = Statistics::PCA->new();
    
    $pca->load_data( { format => 'table', data => $data } );
    $pca->pca();
    
    return $pca->results('eigenvector');
}

sub _true_rot_angle {
    my($opp, $adj) = @_;
    
    my $angle = 0;

    if ($adj != 0.0) {
        
	$angle = atan( - $opp / $adj);
        
	if($opp < 0.0 && $adj > 0.0){
	    $angle += 2*pi;
	}
        
	if($adj < 0.0){
	    $angle += pi;
	}
    }
    else{
	if($opp > 0.0) { # 1st->2nd quadrant boundary
	    $angle = pi / 2;
	}
	else { # 3rd -> 4th quadrant boundary
	    $angle = 3 * (pi / 2);
	}
    }

    if($opp == 0){
	if($adj > 0.0) { # 4th -> 1st quadrant boundary
	    $angle = 0;
	}
	else { # 2nd -> 3rd quadrant boundary
	    $angle = pi;
	}
    }
    return($angle);
}

# Give this subroutine the angle and the axis that you want to rotate about
sub _rotation_matrix {

    my($angle, $axis) = @_;
    $axis = lc($axis);

    if($axis eq 'x'){
	my @rotationXMatrix = ( 1, 0,           0,
                                0, cos($angle),  -sin($angle),
                                0, sin($angle),  cos($angle), );
        
	return \@rotationXMatrix;
    }
    elsif($axis eq 'y'){
	my @rotationYMatrix = ( cos($angle), 0, -sin($angle),
                                0,           1, 0,
                                sin($angle), 0, cos($angle), );

	return \@rotationYMatrix;

    }
    elsif($axis eq 'z') {
	my @rotationZMatrix = ( cos($angle), -sin($angle), 0,
                                sin($angle),  cos($angle), 0,
                                0,            0,           1, );

	return \@rotationZMatrix;

    }
    else {
	croak "ERROR: Axis '$axis' for rotation is not valid!\n";
    }
}

sub _rotate {

    my( $vector, $matrix_ref ) = @_;
    
    my ($x, $y, $z) = @{ $vector };
    
    my @matrix = @{ $matrix_ref };

    my $xd = ( $x*$matrix[0] + $y*$matrix[1] + $z*$matrix[2] );
    my $yd = ( $x*$matrix[3] + $y*$matrix[4] + $z*$matrix[5] );
    my $zd = ( $x*$matrix[6] + $y*$matrix[7] + $z*$matrix[8] );

    return [$xd, $yd, $zd];
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

    my $zpos_count = 0;
    my $zneg_count = 0;
    
    # Only consider non-patch atoms within x and y ranges to form min set
    foreach my $np_atom ( values %{ $np_atom_h } ) {
        next unless   $np_atom->x < $limit{xmax}
                   && $np_atom->x > $limit{xmin}
                   && $np_atom->y < $limit{ymax}
                   && $np_atom->y > $limit{ymin};

        foreach my $p_atom ( values %{ $p_atom_h } ) {
            my $dist_thresh =  $np_atom->radius + $p_atom->radius;

            my $distance = sqrt (   ( $np_atom->x - $p_atom->x )**2
                                  + ( $np_atom->y - $p_atom->y )**2
                                  + ( $np_atom->z - $p_atom->z )**2 );
            
            if ($distance < $dist_thresh) {
                #print  "Atom " . $np_atom->resSeq . $np_atom->name
                #     . " in contact with atom "   . $p_atom->resSeq
                #     . $p_atom->name . "\n";
                
                abs $np_atom->z == $np_atom->z ?
                    ++$zpos_count : ++$zneg_count;
            }
        }
    }
    print "zpos: $zpos_count zneg: $zneg_count\n";
    return ( $zpos_count, $zneg_count );
}

sub _residue_order {
    my $self = shift;
    
    my @calpha
        = map { $self->parent->atom_array->[$self->parent->resid_index->{$_}->{CA}] }
            keys %{ $self->patch->resid_index };
    
    my %angle = ();
    
    foreach my $p_atom (@calpha) {

        my $vector = vector($p_atom->x, $p_atom->y, $p_atom->z);
        
        my $angle = $self->_true_angle_from_x($vector);

        print $p_atom . "\n";
        print $p_atom->resSeq . ': ' .  $angle * (180 / pi) . "\n";
       
        $angle{$p_atom->resSeq} = $angle;      
    }

    my @sorted = sort { $angle{$a} <=> $angle{$b} } keys %angle;

    print "@sorted\n";
}

# Sub to deal with calulating obtuse angles by calculating angle on x y plane
sub _true_angle_from_x {
    my $self = shift;
    my($v) = @_;
    
    croak "$v is not a vector object" if ref $v ne 'Math::VectorReal';

    # Set z component to zero and normalise
    my $v2d = vector($v->x, $v->y, 0 );
    
    # Avoid attempts to divide by zero
    return 0 if length($v2d) == 0;
 
    my $inner_prod = ( $v2d . X ) / $v2d->length * 1; # Length X = 1

    my $angle = acos ($inner_prod);
    
    if ($v2d->y < 0) {
            $angle = (2 * pi) - $angle;
        }
    return $angle;    
}

# Given a vector, returns a rotation matrix to rotate vector to x 
sub _RMforvector2x {
    my $self = shift;
    my $vector = shift;

    croak "$vector is not a vector object"
        if ref $vector ne 'Math::VectorReal';

    my($x, $y, $z) = ( $vector->x, $vector->y, $vector->z );

    my $ang1;
    
    if ( 0 == $y ) {
        $ang1 = pi / 2;
    }
    else {
        $ang1 = atan ( - $z / $y );
    }

    my $ang2;
    
    if ( 0 == $x ) {
        $ang2 = pi / 2;
        
    }
    else {
        $ang2 = atan ( - ( cos($ang1)*$y + sin($ang1)*$z ) / $x );
    }
    
    my $RM
        = Math::MatrixReal->new_from_rows(
            [ [ cos($ang2), -sin($ang2)*cos($ang1),  sin($ang2)*sin($ang1) ],
              [ sin($ang2),  cos($ang2)*cos($ang1), -cos($ang2)*sin($ang1),],
              [ 0         ,  sin($ang1)           ,  cos($ang1),           ],
          ],
        );
    return $RM;
}



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
