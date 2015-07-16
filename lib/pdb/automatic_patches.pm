package automatic_patches;
use Moose;
use Moose::Util::TypeConstraints;

use types;

use Carp;
use pdb::pdb;
use pdb::makepatch;
use pdb::pdb2xmas;
use pdb::pdbFunctions;
use pdb::multiChain;
use GLOBAL qw(&rm_trail);

use Parallel::ForkManager;

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

subtype 'ArrayRefOfValidPDBObjects',
    as 'ArrayRef[ValidPDBObject]';

coerce 'ArrayRefOfValidPDBObjects',
    from 'ValidPDBObject',
    via { [$_] };

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
    isa => 'ArrayRefOfValidPDBObjects',
    lazy => 1,
    coerce => 1,
    builder => '_build_pdb_object',
);

has 'pdb_code' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
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

has 'ASA_type' => (
    is => 'rw',
    isa => 'Str',
    builder => '_build_ASA_type',
    lazy => 1,
);

for my $name (qw(pdb xmas)) {
    my $att_name = $name . '_file';
    my $builder = "_build_$name" . "_fname";
    
    has $att_name => (is => 'rw',
                      isa => 'FileReadable',
                      builder => $builder,
                      lazy => 1);
}

has 'build_patches_from_parent' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

has 'forkFlag' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

has 'numForks' => (
    is => 'rw',
    isa => 'Int',
    default => `nproc` - 1,
);


# Methods

# Build an automatic_patches object straight from a pdb object
around BUILDARGS => sub {
    my $orig = shift;
    my $class = shift;

    my %arg = @_;

    if ( exists $arg{pdb_object} ) {

        my $pdb_obj
            = ref $arg{pdb_object} eq 'ARRAY' ? $arg{pdb_object}->[0]
            : $arg{pdb_object};

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

sub _build_ASA_type {
    my $self = shift;

    # Assign ASA type based on the pdb objects assigned
    # If array of chains, assume ASAb
    # Otherwise, use object class of first element to set ASA type

    if (@{$self->pdb_object()} > 1) {
        my $multipleChainArrayTest = 1;
        foreach my $pdb_obj (@{$self->pdb_object()}) {
            if (ref $pdb_obj ne 'chain') {
                $multipleChainArrayTest = 0;
                last;
            }
        }
        return 'ASAb' if $multipleChainArrayTest;
    }
    
    my $type
        = ref $self->pdb_object->[0] eq 'pdb' ? 'ASAc'
            : 'ASAb';
        
    return $type;
}

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

        my $pdb2xmas = pdb::pdb2xmas->new(pdb_file => $pdb_file);
        
        my @output = $pdb2xmas->output();
        
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
    my %arg  = @_;
    
    my $pdb_code = lc $self->pdb_code();

    my $class = $self->has_chain_id ? 'chain' : 'pdb';
    my $form = $class eq 'chain' ? 'monomer' : 'multimer' ; 
    
    my $pdb_obj_aref = $self->pdb_object();

    # Use pdb::multiChain::readASAb to read ASAs if ASA type = ASAb
    # and ASA have not yet been read
    if ($self->ASA_type() eq 'ASAb' && ! $self->pdb_object->[0]->has_read_ASA) {
        pdb::multiChain::readASAb($self->pdb_object());
    }
    
    foreach my $pdb_obj (@{$pdb_obj_aref}) {
        if (! $pdb_obj->has_read_ASA()) {
            
            my @ASA_read_err = ();
        
            # Read ASA values for pdb object, check for errors
            foreach my $ret ($pdb_obj->read_ASA()) {
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
    }
        
    # Create tmp pdb file with modified atom lines ala xmas2pdb output
    my $ASA_type = $self->ASA_type();
    my $predicate = 'has_' . $ASA_type;
    
    my %swap = (occupancy => 'radius', tempFactor => $ASA_type);

    my @ATOM_lines = ();

    my $atomAref = pdb::pdbFunctions::generateAtomAref(@{$self->pdb_object});

    foreach my $atom (@{$atomAref}) {
        
        # Avoid printing atoms to file that do not have ASA or are labelled
        # solvent
        next if ! $atom->$predicate || $atom->is_solvent();
        
        # Get atom string with occupancy replaced with radius, and tempFactor
        # replaced with ASA
        push(@ATOM_lines, $atom->stringify(\%swap));
        if ($atom->is_terminal()) {
            push(@ATOM_lines, $atom->stringify_ter());
        }
    }
    
    my $tmp_file_name
        = write2tmp->new(suffix => '.pdb', data => \@ATOM_lines)->file_name();
            
    my $all_pc_errors = [];
    my $all_patch_centres = [];

    if (exists $arg{patch_centres} && @{$arg{patch_centres}}) {
        $all_patch_centres = $arg{patch_centres};
    }
    else {
        $all_patch_centres = [$self->_get_patch_centres()];
    }
    
    my $patchAref = $self->forkMakePatch($all_patch_centres, $tmp_file_name);

    return @{$patchAref};
}

sub _get_patch_centres {
    my $self = shift;
    my @all_patch_centres = ();
    foreach my $pdb_obj (@{$self->pdb_object()}) {
        my($pc_errors, $patch_centres)
            = $pdb_obj->patch_centres(type => $self->ASA_type());
        
        push(@all_patch_centres, @{$patch_centres});
    }
    return @all_patch_centres;
}

sub forkMakePatch {
    my $self = shift;
    
    my $patchCentreAref = shift;
    my $pdbDataFile = shift;
    
    my $forkFlag = $self->forkFlag();
    
    my %mPatchArgs = (makepatch_file => $makepatch,
                      patch_type     => $self->patch_type,
                      radius         => $self->radius,
                      pdb_file       => $pdbDataFile,
                      pdb_code       => lc $self->pdb_code,
                      surf_min       => $self->surf_min,
                  );
   
    if ($self->build_patches_from_parent()) {
        $mPatchArgs{pdb_object} = $self->pdb_object();
        $mPatchArgs{new_atoms} = 0;
    }

    # If parent pdb is actually a ref to an array of chains, then create a
    # multi-chain atom hash. Otherwise, take atom hash of single parent pdb
    my $atomSerialHref
        = ref $self->pdb_object eq 'ARRAY' ?
            pdb::multiChain::multiChainAtomSerialHref($self->pdb_object())
              : $self->pdb_object->atom_serial_hash();

    my $pm = $forkFlag ? Parallel::ForkManager->new($self->numForks) : 0;

    my $pProcID = $$;
    
    for (my $i = 0 ; $i < @{$patchCentreAref} ; ++$i) {
        
        my $pid;

        if ($pm) {
            $pid = $pm->start and next;
        }

        my $centAtom = $patchCentreAref->[$i];
        
        $mPatchArgs{central_atom} = $centAtom;
        
        my $mkpObj = makepatch->new(%mPatchArgs);

        my @atomSerials = ();
        
        foreach my $atomLine (@{$mkpObj->output()}) {
            # Get atom serial
            my $atomSerial = rm_trail(substr($atomLine, 6, 5));
            push(@atomSerials, $atomSerial);
        }

        my $fName = "/tmp/$pProcID.$i";
        open(my $OUT, ">", $fName)
            or die "Cannot open out file $fName, $!";
        
        print {$OUT} "@atomSerials\n";
             
        $pm->finish if $pm;
    }

    $pm->wait_all_children if $pm;

    my @patches = ();
    
    foreach (my $i = 0 ; $i < @{$patchCentreAref} ; ++$i){

        my $inFile = "/tmp/$pProcID.$i";
        open(my $IN, "<", $inFile) or die $!;
        
        my $serialLine = <$IN>;
        chomp $serialLine;
        my $atomSerialAref = [split(" ", $serialLine)];
        
        my $centralAtom = $patchCentreAref->[$i];

        # Get corresponding atom object for patch atom serial
        my @atoms = map {$atomSerialHref->{$_}} @{$atomSerialAref};
        
        my %patchArgs = (pdb_code => $self->pdb_code(),
                         parent_pdb => $self->pdb_object,
                         central_atom => $centralAtom,
                         atom_array => \@atoms);
        
        my $patch = patch->new(%patchArgs);
        
        push(@patches, $patch);
        close $IN;
        unlink $inFile;
    }    
    return \@patches;
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
