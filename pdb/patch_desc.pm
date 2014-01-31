package patch_desc;

use Moose;
use Moose::Util::TypeConstraints;
use types;

use Math::Trig;
use Math::VectorReal qw(:all);
use Math::MatrixReal;
use Statistics::PCA;

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

    # If patch only contains one, two, or three residues
    if ( scalar keys %{ $self->patch->resid_index } < 4 ){

        my @simple_residue_order = ();
        
        my $central_atom = $self->patch->central_atom;

        push( @simple_residue_order, $central_atom->resName() );
        
        foreach my $resid ( keys %{ $self->patch->resid_index() } ) {
            
            my %atom_names = %{ $self->patch->resid_index->{$resid} };
            
            my $index = [ values %atom_names ]->[0];
            my $p_atom = $self->patch->atom_array->[$index];
            
            next if $p_atom->resid eq $central_atom->resid;
            
            push( @simple_residue_order, $p_atom->resName() );
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
    
    # Default contact threshold - if one side of patch has < threshold
    # contacts then it is considered a surface side and a patch order for
    # that side will be returned
    my $threshold = 4;
    
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


    # Keep record of all original parent atom co-ords to reset to later
    my %orig_coords = ();
    foreach my $atom ( @{ $self->parent->atom_array() } ) {
        $orig_coords{ $atom->serial() }
            = { x => $atom->x(), y => $atom->y(), z => $atom->z() }; 
    }
    
    # Avoid inclusion of solvent atoms in all hashes
    
    my ($surface_errors, $surface_atoms) = $self->_surface_atoms;

    my %no_ASA_atom = ();
    
    # Return error if any patch atoms do not have ASA defined
    foreach my $error ( @{$surface_errors} ) {
        if ( $error->type() eq 'atom_no_ASAm' ) {
            my $atom = $error->data->{atom};
            $no_ASA_atom{ $atom->serial() } = $atom;
        }
    }
    
    $surface_atoms = [ grep { ! $_->is_solvent() } @{ $surface_atoms } ];
   
    my %all_atom
        = map { $_->serial => $_ }
            grep { ! $_->is_solvent() } @{ $self->parent->atom_array };
        
    my %patch_atom
        = map { $_->serial => $all_atom{$_->serial} }
            @{ $self->patch->atom_array };

    foreach my $serial ( keys %patch_atom ) {
        if ( grep /^$serial$/, keys %no_ASA_atom ) {
            my $message
                = "patch desc: patch atom " . $serial . " has no ASA value";

            my $atom = $patch_atom{$serial};
            
            my $error = local::error->new( message => $message,
                                           type => 'patch_atom_no_ASA',
                                           data => { atom => $atom }, );

            $self->reset_parent_coords(\%orig_coords);
            return $error;
        }
    }
    
    my %patch_surf_atom
        = map { $_->serial => $_ }
            grep ( defined $patch_atom{$_->serial}, @{ $surface_atoms } );
    
    my %nonpatch_atom
        = map { $all_atom{$_}->serial => $all_atom{$_} }
            grep( ! defined $patch_atom{$_}, keys %all_atom );

    # Set central atom to origin
    my %cent_atom_coord = ( x => $self->patch->central_atom->x(),
                            y => $self->patch->central_atom->y(),
                            z => $self->patch->central_atom->z(), );
    
    foreach my $atom (values %all_atom) {
        foreach my $coord ( 'x', 'y', 'z') {
            $atom->$coord( $atom->$coord - $cent_atom_coord{$coord} ); 
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
    
    my($posz_contacts, $negz_contacts)
        = $self->_surface_sides( \%patch_surf_atom, \%all_atom );

    if ( ref $posz_contacts eq 'local::error' ) {
        my $message
            = "Could not determine surface side: "
             . $posz_contacts->message;

        my $error = local::error->new( message => $message,
                                       type => 'surface_side',
                                       parent => $posz_contacts );

        $self->reset_parent_coords(\%orig_coords);
        return $error;
    }

    my @return = ();
    
    my @residue_order = $self->_residue_order;

    my @aacid_order = ();

    foreach my $res (@residue_order) {
        
        my %atom_names = %{ $self->patch->resid_index->{$res} };
        
        my $index = [ values %atom_names ]->[0];
        my $resName = $self->patch->atom_array->[$index]->resName();

        my $onelc = eval { three2one_lc($resName) };
        $onelc = 'X' if ! $onelc;
        push(@aacid_order, $onelc);
    }   

    my $string = join( '', @aacid_order );
    
    my $rev_string
        = shift (@aacid_order) . reverse ( join( '', @aacid_order) );
    
    if (  abs ( scalar @{$posz_contacts} - scalar @{$negz_contacts} )
              < $threshold ) {
        push( @return, ($string, $rev_string) );
        
    }elsif ( scalar @{$posz_contacts} < scalar @{$negz_contacts} ) {
        push( @return, $string );
    }
    else {
        push( @return, $rev_string );
    }
    
    $self->reset_parent_coords(\%orig_coords);
    return @return;
}

sub _surface_atoms {
    my $self = shift;

    my @surface_atoms = ();
    my @errors = ();
    
    if ( ! $self->parent->has_read_ASA ) {
        my $ret = $self->parent->read_ASA();

        if ( ref $ret eq 'local::error' ) {
            push(@errors, $ret);
        }
        elsif ( ref $ret eq 'ARRAY' ){
            foreach ( @{$ret} ) {
                if ( ref $_ eq 'local::error' ){
                    push(@errors, $_);
                }
            }
        }
    }
    
    foreach my $atom ( @{ $self->parent->atom_array } ) {

        if ( ! $atom->has_ASAm() ){
            my $message
                = "atom " . $atom->serial() . " does not have ASAm assigned";
            
            my $error = local::error->new( message => $message,
                                           type => 'atom_no_ASAm',
                                           data => { atom => $atom, }
                                       );
            push(@errors, $error);
            next;
        }
        
        if ( $atom->ASAm > 0 ) {
            push(@surface_atoms, $atom);
        }
    }
    return ( \@errors, \@surface_atoms );
}

# Returns number of atoms contacting either side of patch
# 
sub _surface_sides {
    my ( $self, $ps_atom_h, $all_atom_h ) = @_;

    my @p_atom_serial = map { $_->serial  } values %{ $ps_atom_h };
    
    my %limit = ( xmin => 0, xmax => 0, ymin => 0, ymax => 0 );

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
        
        if ( ! $atom->has_radius  || ! defined $atom->radius ) {
            my $message
                =  "non-patch atom " . $atom->serial
                 . " has no radius set";

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

            if ( ! $ps_atom->has_radius || ! defined $atom->radius) {
                my $message
                    =  "patch atom " . $atom->serial
                        . " has no radius set";
                
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

sub reset_parent_coords {
    my $self = shift;
    my $orig_coords = shift;

    #return 1; #TEST
    
    foreach my $atom ( @{ $self->parent->atom_array() } ) {
        foreach my $coord ( qw ( x y z ) ) {
            $atom->$coord( $orig_coords->{$atom->serial}->{$coord} );
        }
    }
}

__PACKAGE__->meta->make_immutable;

1;
__END__
1
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
