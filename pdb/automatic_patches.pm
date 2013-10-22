package automatic_patches;

use Moose;
use Moose::Util::TypeConstraints;

use types;

use Carp;
use pdb::pdb;
use pdb::xmas2pdb;
use pdb::makepatch;

use write2tmp;

use TCNPerlVars;

# Subtypes

subtype 'PatchType',
    as 'Str',
    where { $_ =~ m{ \A (?: contact|normal ) \s* \z }xms },
    message { "$_ is not a valid patch type" };

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

has 'pdb_code' => (
    is => 'rw',
    isa => 'ValidPDB',
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

sub get_patches {
    my $self = shift;

    my $pdb_code = lc $self->pdb_code();
    
    # Get files
    my $pdb_file  = $self->pdb_file();
    croak "Could not read $pdb_file" if ! -r $pdb_file;

    my $xmas_file = $self->xmas_file();
    croak "Could not read $xmas_file" if ! -r $xmas_file;
    
    my $class = $self->has_chain_id ? 'chain' : 'pdb';

    my %pdb_arg
        = ( pdb_code => $pdb_code,
            pdb_file => $pdb_file,
            xmas_file => $xmas_file,
        );

    my $form;
    
    if ($class eq 'chain') {
        $pdb_arg{chain_id} = uc $self->chain_id();
        $form = 'monomer';
    }
    else {
        $form = 'multimer';
    }
    
    my $pdb_obj = $class->new(%pdb_arg);
    
    my %x2p_arg
        = ( xmas2pdb_file => $xmas2pdb,
            radii_file    => $radii_file,
            xmas_file     => $xmas_file,
            form          => $form,
        );
    
    # xmas2pdb
    my $x2p_obj = xmas2pdb->new(%x2p_arg);

    # Read ASA values for pdb object
    $pdb_obj->read_ASA($x2p_obj);

    # Parse x2p output to avoid making patches with non-chain res
    # if object is chain

    my $tmp_pdb_file;
    
    if ($class eq 'chain') {

        my @ATOM_lines = ();
        
        foreach my $line ( @{ $x2p_obj->output } ) {
            if (   $line =~ /^ATOM/ 
                && substr($line, 21, 1) eq $pdb_obj->chain_id ){
                push( @ATOM_lines, $line);
            }
        }

        croak   "No ATOM lines parsed with chain id " . $pdb_obj->chain_id
            if ! @ATOM_lines;

        my $tmp_file = write2tmp->new( data => [@ATOM_lines], suffix => '.pdb' );
        
        $tmp_pdb_file  = $tmp_file->file_name;
    }
    else {
        $tmp_pdb_file = $x2p_obj->output_file();
    }
    
    my @patches = ();

    foreach my $atom_index ( $pdb_obj->patch_centres() ) {

        my $cent_atom = $pdb_obj->atom_array->[$atom_index];

        my %mkp_arg
            = ( makepatch_file => $makepatch,
                patch_type     => $self->patch_type,
                radius         => $self->radius,
                pdb_file       => $tmp_pdb_file,
                pdb_code       => $pdb_code,
                central_atom   => $cent_atom,
                surf_min       => $self->surf_min,
            );

        my $mkp_obj = makepatch->new(%mkp_arg);
        my $patch   = patch->new($mkp_obj);
        push(@patches, $patch);
        
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
