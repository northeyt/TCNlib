package automatic_patches;

use Moose;
use Moose::Util::TypeConstraints;

use types;

use Carp;
use pdb::pdb;
use pdb::makepatch;

use write2tmp;

use TCNPerlVars;

# Subtypes

subtype 'PatchType',
    as 'Str',
    where { $_ =~ m{ \A (?: contact|normal ) \s* \z }xms },
    message { "$_ is not a valid patch type" };

subtype 'ValidPDBObject',
    as 'Ref',
    where { ref $_ eq 'pdb' || ref $_ eq 'chain' },
    message { "$_ is not a valid pdb object (pdb or chain)" };

# Import vars
my $pdbprep = $TCNPerlVars::pdbprep;
my $pdbext  = $TCNPerlVars::pdbext;
my $pdbdir  = $TCNPerlVars::pdbdir;

my $xmasprep = $TCNPerlVars::xmasprep;
my $xmasext  = $TCNPerlVars::xmasext;
my $xmasdir  = $TCNPerlVars::xmasdir;

my $makepatch = $TCNPerlVars::makepatch;
my $xmas2pdb  = $TCNPerlVars::xmas2pdb;
my $radii_file = $TCNPerlVars::radii_file;

my $pdb2xmas = $TCNPerlVars::pdb2xmas;

my $tmpdir = $TCNPerlVars::tmpdir . '/automatic_patches' ;

# Attributes

has 'radius' => (
    is => 'rw',
    isa => 'Int',
    required => 1,
);

has 'patch_type' => (
    is => 'rw',
    isa => 'PatchType',
    required => 1,
);

has 'pdb_object' => (
    is => 'rw',
    isa => 'ValidPDBObject',
    lazy => 1,
    builder => '_build_pdb_object',
);

has 'pdb_code' => (
    is => 'rw',
    isa => 'Str',
    required => 1
);

has 'chain_id' => (
    is => 'rw',
    isa => 'ValidChar',
    predicate => 'has_chain_id'
);

# Minimum absolute accessibility an atom must have if it is to be considered
# surface
has surf_min => (
    is => 'rw',
    isa => 'Num',
    default => 0,
);

for my $name ( 'pdb', 'xmas' ) {
    my $att_name = $name . '_file';
    my $builder = "_build_$name" . "_fname";
    
    has $att_name => ( is => 'rw',
                       isa => 'FileReadable',
                       builder => $builder,
                       lazy => 1,
                   );
    
}

# Methods

# Build an automatic_patches object straight from a pdb object
around BUILDARGS => sub {
    my $orig = shift;
    my $class = shift;

    my %arg = @_;

    if ( exists $arg{pdb_object} ) {
        my $pdb_obj = $arg{pdb_object};
        
        $arg{pdb_code} = $pdb_obj->pdb_code;
        $arg{chain_id} = $pdb_obj->chain_id if ref $pdb_obj eq 'chain';
        
        foreach my $type ( 'pdb', 'xmas' ) {
            my $attribute = $type . '_file';
            my $predicate = 'has_' . $attribute;
            if ($pdb_obj->$predicate) {
                $arg{$attribute} = $pdb_obj->$attribute;
            }
        }
    }
    return $class->$orig(%arg);
};

sub _build_pdb_fname {
    my $self = shift;

    if (@_) {
        return $_[0];
    }
    
    my $pdb_code = $self->pdb_code();
    my $fname = $pdbprep . lc $pdb_code . $pdbext;
    croak "no file found for $pdb_code in $pdbdir"
        if ! -e $fname;
    return $fname;
    
}
sub _build_xmas_fname {
    my $self = shift;

     if (@_) {
        return $_[0];
    }
    
    my $pdb_code = $self->pdb_code();
    my $fname = $xmasprep . lc $pdb_code . $xmasext;

    # Attempt to build xmas file from pdb file if xmas file not found
    if ( ! -e $fname) {
        my $pdb_file = $self->pdb_file();
        my $cmd = "$pdb2xmas -q -f $pdb_file";
        my @output = `$cmd`;
        croak "pdb2xmas produced no output. Command run:\n$cmd"
            if ! @output;
        
        my $tmp
            = write2tmp->new( data => [ @output ],
                               SUFFIX => '.xmas',
                          );

    $fname = $tmp->file_name;
    }
    return $fname;
}

sub _build_pdb_object {
    my $self = shift;
    
    my $class = $self->has_chain_id ? 'chain' : 'pdb';
    
    my %pdb_arg
        = ( pdb_code => $self->pdb_code,
            pdb_file => $self->pdb_file(),
            xmas_file => $self->xmas_file(),
            hydrogen_cleanup => 1,
            altLoc_cleanup => 1,
            solvent_cleanup => 1,
        );
    
    $pdb_arg{chain_id} = $self->chain_id() if $class eq 'chain';
        
    return $class->new(%pdb_arg);    
}

sub get_patches {
    my $self = shift;

    my $pdb_code = lc $self->pdb_code();
    
    # Get xmas file

    my $xmas_file = $self->xmas_file();
    croak "Could not read $xmas_file" if ! -r $xmas_file;
    
    my $class = $self->has_chain_id ? 'chain' : 'pdb';
    my $form = $class eq 'chain' ? 'monomer' : 'multimer' ; 
    

    my $pdb_obj = $self->pdb_object();
    
    if ( ! $pdb_obj->has_read_ASA() ) {
        
        $pdb_obj->read_ASA();
       
        my @ASA_read_err = ();
        
        # Read ASA values for pdb object, check for errors
        foreach my $ret ( $pdb_obj->read_ASA() ) {
            if (ref $ret eq 'local::error'){
                if ( $ret->type() eq 'ASA_read' ) {
                    push(@ASA_read_err, $ret);
                }
                else {
                    croak "Unrecognised error type '" . $ret->error()
                    . "' returned by read_ASA";
                }
            }
        }
    }
    
    # Create tmp pdb file with modified atom lines ala xmas2pdb output
    my $ASA_type = $form eq 'monomer' ? 'ASAm' : 'ASAc' ;
    my $predicate = 'has_' . $ASA_type;
    
    my %swap = ( occupancy => 'radius', tempFactor => $ASA_type );

    my @ATOM_lines = ();
    
    foreach my $atom ( @{ $pdb_obj->atom_array() } ) {

        # Avoid printing atoms to file that do not have ASA or are labelled
        # solvent
        next if ! $atom->$predicate || $atom->is_solvent();

        
        push( @ATOM_lines, $atom->stringify( {%swap} ) );
        if ( $atom->is_terminal() ) {
            push( @ATOM_lines, $atom->stringify_ter() );
        }
    }
    
    my $tmp = write2tmp->new( suffix => '.pdb',
                              data => [ @ATOM_lines ],
                          );
    
    my $tmp_file_name = $tmp->file_name();
    
    my @patches = ();

    my($pc_errors, $patch_centres) = $pdb_obj->patch_centres();
    
    foreach my $atom_index ( @{$patch_centres} ) {

        my $cent_atom = $pdb_obj->atom_array->[$atom_index];

        my %mkp_arg
            = ( makepatch_file => $makepatch,
                patch_type     => $self->patch_type,
                radius         => $self->radius,
                pdb_file       => $tmp_file_name,
                pdb_code       => $pdb_code,
                central_atom   => $cent_atom,
                surf_min       => $self->surf_min,
            );

        my $mkp_obj = makepatch->new(%mkp_arg);
    
        my $return = do {
            local $@;
            my $ret;
            eval { $ret = patch->new($mkp_obj); 1 };
            $ret ? $ret : $@;
        };
        
        push(@patches, $return);
    }
    return  @patches;
}


__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

automatic_patches - Perl extension for blah blah blah

=head1 SYNOPSIS

   use automatic_patches;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for automatic_patches, 

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
