package pdb;

=head1 NAME

pdb - A class for creating pdb objects from pdb files.

=cut

=Head SYNOPSIS

use pdb;
$pdbObject = pdb->new(pdb_code => '1adjs',
                         pdb_file => '1adjs.pdb');

# Or ...
$pdbObject = pdb->new(pdb_code => '1djs'); # Grab file automatically!

=cut


=head1 DESCRIPTION

pdb is one of the classes I have used for the majority of my PhD work. It is a
class whose primary function is to parse information from PDB files in order
to create atom objects, then supply a range of attributes and methods for
accessing those atoms. pdb is the parent class of the chain and patch objects

=cut

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;
use TCNUtil::local::error;

use TCNUtil::GLOBAL qw(&rm_trail &three2one_lc &is_int);
use TCNUtil::VectorCalcs qw(rotate2pc
                          get_rms_difference_of_points_from_bestfit_plane
                     );

use Carp;

use TryCatch;
use Storable;
use Data::Dumper;

use Math::VectorReal qw(:all);
use Math::MatrixReal;
use Math::Trig;

use pdb::get_files;
use pdb::solv;
use pdb::ss;
use pdb::secstrCalculator;
use pdb::hbondFinder;
use pdb::RadiusFinder;

=head1 Methods

=over 12

=cut

### Attributes #################################################################
################################################################################

=item C<pdb_code>

Returns pdb code of object

=cut

has 'pdb_code' => (
    isa => 'Str',
    is  => 'rw',
    builder => '_build_pdb_code',
    lazy => 1,
);

=item C<pdb_file>

Returns pdb filename of object. This can be automatically assigned by
pdb::get_files if pdb_code is supplied by user during object construction.

=cut

=item C<pdb_data>

Ref to array of lines from object's pdb file.

=cut

# pdb_file, pdb_data, xmas_file, xmas_data atributes. If any are not supplied,
# builder methods will attempt to find files using pdb::get_files

has 'pdb_file' => (
    isa => 'FileReadable',
    is  => 'rw',
    predicate => 'has_pdb_file',
    builder => '_get_pbd_file',
    lazy => 1,
);

has 'pdb_data' => (
    isa => 'ArrayRef',
    is => 'rw',
    builder => '_build_pdb_data',
    predicate => 'has_pdb_data',
    lazy => 1,
);

=item C<remark_hash>

Hash of remark lines parsed from pdb data
Hash is of form:
    remarkNumber => Ref to array of scalar refs to pdb_data lines

=cut

has 'remark_hash' => (
    isa => 'HashRef',
    is => 'ro',
    builder => '_build_remark_hash',
    lazy => 1,
);

=item C<atom_array>

Ref to array of atom objects parsed from PDB data, or assigned by user

=cut

has 'atom_array' => (
    isa => 'ArrayRef[atom]',
    is  => 'rw',
    lazy => 1,
    builder => '_parse_atoms',
);

=item C<terminal_atom_array>

Ref to array of atoms labelled is_terminal()

=cut

has 'terminal_atom_array' => (
    isa => 'ArrayRef[atom]',
    is  => 'rw',
    builder => '_build_terminal_atom_array',
    lazy => 1,
);

=item C<atom_index>

Ref to hash of atoms, indexed as so:
    chainID -> resSeq -> atom_name -> atom

e.g.
    $firstCAlpha = $pdb->atom_index->{A}->{1}->{CA};

=cut

has 'atom_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_atom_index',
);

=item C<resid_index>

Ref to hash of atoms, indexed as so:
    resid -> atom_name -> atom

Where resid = chainID.resSeq (e.g. A.12)

=cut

has 'resid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_resid_index',
);

=item C<missing_resid_index>

Index in form of resid -> {} (i.e. ref to empty hash)
Where each resid has been identified as missing from structure,
e.g. information parsed from pdb file remarks

=cut

has 'missing_resid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_missing_resid_index',    
);

=item C<pdbresid_index>

This index is the same as resid_index, except that resids are prepended with
the object's pdb code, e.g. 4houA.213

=cut

has 'pdbresid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_pdbresid_index',
);

=item C<atom_serial_hash>

Hash of form atomSerial => atom, e.g.
    $firstAtom = $pdb->atom_serial_hash->{1};

=cut

has 'atom_serial_hash' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_atom_serial_hash',
);

=item C<multi_resName_resid>

Hash of resids with multiple resNames with form resid => 1, e.g.
    @multiNameResids = keys %{$pdb->multi_resName_resids()};

=cut

has 'multi_resName_resids' => (
    isa => 'HashRef',
    is => 'rw',
    lazy => 1,
    builder => '_build_multi_resName_resids',
);

=item C<altLoc_cleanup(BOOL)>

Flag that if set to TRUE, then atom altLocs will be processed so that only the
altLoc for a given residue atom with highest occupancy will be added to the
object's atom_array

DEFAULT = TRUE

=cut

has 'altLoc_cleanup' => (
    isa => 'Bool',
    is  => 'rw',
    default => 1,
);

=item C<hydrogen_cleanup(BOOL)>

Flag that if set to TRUE, then no hydrogen atoms will be included in
object's atom_array.

DEFAULT = 0

=cut

has 'hydrogen_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

=item C<het_atom_cleanup(BOOL)>

Flag that if set to TRUE, then not HETATM atoms will be included in the object's
atom_array.

=cut

has 'het_atom_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

=item C<solvent_cleanup>

Flag that if set to TRUE, then not solvent atoms will be included in the
object's atom_array.

=cut

# Avoid parsing solvent atoms from pdb data
has 'solvent_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

# Avoid parsing solvent atoms from pdb data
has 'solvent_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);



# Used by methods to indicate that the object has had
# accessible surface area (ASA) measures taken in some form
has 'has_read_ASA' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

=item C<experimental_method>

Experimental method of structure determination, parsed from pdb data.

=cut

has 'experimental_method' => (
    is => 'rw',
    isa => enum([ ('X-RAY DIFFRACTION', 'NMR', 'ELECTRON MICROSCOPY') ]),
    predicate => 'has_experimental_method',
    lazy => 1,
    builder => '_build_experimental_method',
);

=item C<resolution>

Resolution of structure, parsed from pdb data.

=cut

has 'resolution' => (
    is => 'rw',
    isa => 'Num',
    predicate => 'has_resolution',
    lazy => 1,
    builder => '_build_resolution',
);

=item C<r_value>

R-value of strucure, parsed from pdb data.

=cut

has 'r_value' => (
    is => 'rw',
    isa => 'Num',
    predicate => 'has_r_value',
    lazy => 1,
    builder => '_build_r_value',
);


has 'missing_residues' => (
    is => 'rw',
    isa => 'HashRef',
    builder => '_parse_remark465',
    lazy => 1,
);


=item C<resid2RelASAHref>

Where RelASA = Relative solvent accessible surface area
Hash of form resid -> rASA, normally assigned using read_ASA method.

=cut

has 'resid2RelASAHref' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
    lazy => 1,
);

=item C<resid2ModResAref>

This method returns a ref to a hash of form
resID => [modifiedResidueName, standardResidueName]
e.g. A.281 => [MSE, MET]
This is used by get_sequence to convert modified residues names to standard
residue names.

=cut

has 'resid2ModResAref' => (
    is => 'rw',
    isa => 'HashRef',
    builder => '_parseMODRESLines',
    lazy => 1,
);

=item C<resNameMod2StdHref>

This method returns a ref a  hash of form modifiedResidueName => StandardResidueName
e.g. MSE => MET

Therefore this hash does not specify the positions of modified residues, simply
the modified to standard residue name mapping.

This hash is derived from the resid2ModResAref by default. It is useful for when
some but not all modified residue ids are listed in the PDB header. Its main
purpose is to be used by get_sequence to attempt to map modified resid residue
names to standard residue names when they have not been listed in MODRES.

It could also be defined by the user before using get_sequence if no PDB header
data is present.

=cut

has 'resNameMod2StdHref' => (
    isa => 'HashRef',
    is  => 'rw',
    lazy => 1,
    builder => '_build_resNameMod2StdHref',
);

has 'threelc2hydrophobicValueHref' => (
    isa => 'HashRef',
    is => 'rw',
    lazy => 1,
    builder => '_buildHydroPhoHrefFromFile',
);

has 'threelc2hydrophobicValueFile' => (
    isa => 'FileReadable',
    is => 'rw',
    lazy => 1,
    default => $TCNPerlVars::hydroPhoValueFile,
);

has 'resID2secStructHref' => (
    isa => 'HashRef',
    is  => 'rw',
    lazy => 1,
    builder => '_build_resID2secStructHref',
);

has 'radiusFinder' => (
    is      => 'rw',
    isa     => 'pdb::RadiusFinder',
    default => sub {pdb::RadiusFinder->new()},
);

has 'atomRadiiAreAssigned' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

### Attribute Builder Methods ##################################################
################################################################################

sub _build_pdb_code {
    my $self = shift;

    if ($self->can("has_parent_pdb") && $self->has_parent_pdb) {
        return $self->parent_pdb->pdb_code();
    }
    else {
        $self->_buildPDBCodeFromData();
    }
}

sub _get_pdb_file {
    my $self = shift;

    my $pdbCode = $self->pdb_code();
    my $fName = eval {pdb::get_files->new(pdb_code => $pdbCode)->pdb_file()};

    if (! $fName) {
        croak "No pdb file specified: "
            . "Attempted to get pdb file from pdb code '$pdbCode', $@";
    }
    
    return $fName;
}

sub _build_pdb_data {
    my $self = shift;
    return $self->_build_data_from_file('pdb');
}

sub _build_data_from_file {
    my($self, $att) = @_;

    my $att_file = $att . '_file';

    my $file = eval {$self->$att_file};
    
    croak "Cannot get file data from $att file - $@"
        if ! $file;

    open(my $fh, '<', $file) or die "Cannot open file '$file', $!";

    my @array = <$fh>;

    croak "No lines found in file $file" if ! @array;
   
    return \@array;
}


# This is the central subroutine that deals with parsing atom lines from pdb
# data in order to create an array of atom objects.
sub _parse_atoms {
    my $self = shift;

    # Get raw ATOM lines from pdb data if ref to array of atoms has not been
    # passed
    my @ATOM_lines = (ref $_[0]) eq 'ARRAY' ? @{$_[0]}
        : $self->_parse_ATOM_lines();
    
    my @atoms = ();

    # Used to cache altLocs
    my %altLoc = ();
    
    # Used to cache terminal labelled atoms
    my %ter = ();
        
    foreach my $line (@ATOM_lines) {

        # Process TER lines to label atoms as terminal
        if ($line =~ /^TER/) {
            # Label previous atom as terminal
            my $prevAtom = $atoms[-1];
            $prevAtom->is_terminal(1);
        
            $ter{$prevAtom->chainID} = []
                if ! exists $ter{$prevAtom->chainID};
            
            push(@{$ter{$prevAtom->chainID}}, $prevAtom->serial());
            next;
        }

        # Parse atom line to create atom object
        my $atom = atom->new(ATOM_line => $line);
        
        # Skip this atom if type is flagged to avoid
        next if ($self->hydrogen_cleanup && $atom->element eq 'H'
                     || $self->het_atom_cleanup && $atom->is_het_atom);

        # Cache altLoc atom if altLoc_clean has been set
        if ($self->altLoc_cleanup && $atom->has_altLoc) {

            my $string = $atom->name . $atom->resName . $atom->chainID
                . $atom->resSeq;

            if (exists $altLoc{$string}) {
                push(@{$altLoc{$string}}, $atom);
            }
            else {
                $altLoc{$string} = [$atom];
                
                # If this is first location of altLocs, add this to atom hash
                # (this makes is possible to assign atoms with altlocs as
                # terminal)
                push(@atoms, $atom);
            }
        }
        else {   
            push(@atoms, $atom);   
        }
    }

    # Process all altLocs if altLoc clean has been set 
    $self->_clean_altLocs(\%altLoc) if $self->altLoc_cleanup;
             
    # Label solvent atoms
    $self->_determineSolventAtoms(\@atoms);
    
    # Avoid including solvent atoms if solventClean has been set 
    if ($self->solvent_cleanup) {
        my @nonsolvent = grep {! $_->is_solvent()} @atoms;
        
        @atoms = @nonsolvent;
    }
    
    croak "No atom objects were created" if ! @atoms;
    
    return \@atoms;
}

sub _clean_altLocs {
    my $self = shift;
    my $altLocHref = shift;
    
    foreach my $altLocAref (values %{$altLocHref}) {
        # Sort altlocs by occupancy
        my @sorted
            = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                    map { [ $_, $_->occupancy ] }
                        @{$altLocAref};

        my $top = $sorted[0];
        
        # Assign top occupant's x,y,z and serial to first element of array
        # (as this is the member that is contained within atom list)
        foreach my $attr (qw(x y z serial)) {
            $altLocAref->[0]->$attr($top->$attr());
        }
            
        # Clear altLoc
        $altLocAref->[0]->clear_altLoc();
    }        
}

sub _determineSolventAtoms {
    my $self = shift;
    my $atomAref = shift;
    
    # Order all atoms by chain, then serial
    my %chain = ();
    foreach my $atom (@{$atomAref}) {
        my $chainID = $atom->chainID();
        
        if (! exists $chain{$chainID}) {
            $chain{$chainID} = {$atom->serial => $atom};
        }
        else {
            $chain{$chainID}->{$atom->serial} = $atom;
        }
    }    
    
    foreach my $chainID (sort keys %chain) {
        my @chain_atoms = ();
        
        foreach my $atom (sort {$a <=> $b} keys %{$chain{$chainID}}) {
            push(@chain_atoms, $chain{$chainID}->{$atom});
        }
        
        # Find terminal atoms in order to properly label solvent atoms
        $self->_determineChainSolventAtoms(@chain_atoms);
    }
}

sub _determineChainSolventAtoms {
    my $self = shift;
    my @chain_atoms = @_;
    
    my @chain_terminal_index = ();
    
    for my $i (0 .. @chain_atoms - 1) {
        my $atom = $chain_atoms[$i];
        my $serial = $atom->serial();
        if (! defined $serial) {
            croak "Undefined $serial!";
        }
        if ($atom->is_terminal()){
            push(@chain_terminal_index, $i);
        }
    }
    
    if(@chain_terminal_index == 1) {
        # If there is one terminal index, then that signals
        # end of non-solvent chain
        my $start = $chain_terminal_index[0];
        
        # Label all atoms after chain terminal as solvent
        for my $i ($start + 1 .. @chain_atoms - 1) {
            $chain_atoms[$i]->is_solvent(1);
            
        }
    }
    elsif (@chain_terminal_index > 1) {
        # Determine which terminal signals end of chain and which signals
        # end of solvent
        my @sortedIndices = sort {$a <=> $b} @chain_terminal_index;
        
        my $index = _determineNonSolvTerIndex(\@chain_atoms, \@sortedIndices);
                        
        for (my $i = 0 ; $i >= @chain_atoms ; $i++) {
            # Label those atoms not in range of returned index and
            # returned index - 1, as solvent
            unless ($i > $chain_terminal_index[$i - 1]
                        && $i < $chain_terminal_index[$i]) {
                $chain_atoms[$i]->is_solvent(1);
            }
        }
    }
}

# Determines the terminal atom that signals the end of the non-solvent chain
# segment
sub _determineNonSolvTerIndex {
    my($chain_atoms, $chain_terminal_index) = @_;

    # Read through each array intersection and determine if any are all
    # HETATM
    
    my @ordered_indices =  @{ $chain_terminal_index };

    # -1 allows +1 to used in foreach range below to avoid including
    # previous terminal atom in proceeding range
    my $prev_index = -1;

    my %hash = ();

    foreach my $index ( @ordered_indices ) {
        my $nonhetflag = 0;
        foreach my $i ( $prev_index + 1 ..$index ) {
            if ( ! $chain_atoms->[$i]->is_het_atom ) {
                $nonhetflag = 1;
            }
        }
        $hash{$index} = $nonhetflag;
        $prev_index = $index;
    }
    
    my $hatom_range_count = 0;
    
    foreach my $index (keys %hash) {
        ++$hatom_range_count if $hash{$index};
    } 
    
    croak "More than one range containing non-HETATM atoms found for chain"
        if $hatom_range_count > 1;
    
    croak "No ranges found containing non-HETATM atoms for chain"
        if $hatom_range_count < 1;

    # Return the index for the atom that signals the end of non-solvent segment
    # of chain
    for my $i ( 0 .. @ordered_indices - 1 ) {
        return $i if $hash{ $ordered_indices[$i] };
    }
}

# Also parses HETATM and TER lines
sub _parse_ATOM_lines {
    
    my $self = shift;

    my %arg = @_;

    my @array = @{ $self->pdb_data };  
 
    my @ATOM_lines = ();
    
    foreach my $line (@array) {
        if ($line =~ /^(?:ATOM|HETATM|TER)/) {
            push(@ATOM_lines, $line);
        }
    }
        
    croak "No ATOM, HETATM or TER lines parsed from pdb data"
        if ! @ATOM_lines;


    if ((exists $arg{all} && $arg{all}) || ! exists $arg{chain_id}) {
        return @ATOM_lines;
    }
    else {
        my @chain_ATOM_lines = ();
        foreach my $line (@ATOM_lines) {
            if (substr($line, 21, 1) eq $arg{chain_id}){
                push(@chain_ATOM_lines, $line);
            }
        }
        croak "No ATOM lines were parsed with chainId " . $arg{chain_id}
            . " for pdb " . $self->pdb_code()
                if ! @chain_ATOM_lines;
        
        return @chain_ATOM_lines;
    }
}

sub _parse_ter {
    my($ter_line) = @_;

    my $serial  = rm_trail( substr( $ter_line, 6, 5 ) );
    croak "Could not parse serial from TER line $ter_line"
        if ! defined $serial;
    
    my $chainID = rm_trail( substr( $ter_line, 21, 1 ) );
    croak "Could not parse chainID from TER line $ter_line"
        if ! defined $chainID;

    return($serial, $chainID);
}

sub _build_multi_resName_resids {
    my $self = shift;

    my %multiResNameResids = ();

    foreach my $resid (keys %{$self->resid_index()}) {
        my @residAtoms = values %{$self->resid_index->{$resid}};

        my %testUnique = map {$_->resName() => $_} @residAtoms;
        
        $multiResNameResids{$resid} = 1 if scalar keys %testUnique > 1;
    }

    return \%multiResNameResids;
}

sub _build_resolution {
    my $self      = shift;
    my $remarkStr = join("", map {${$_}} @{$self->remark_hash()->{2}});
    my ($resValue)  = $remarkStr =~ /RESOLUTION\. \s+ ([0-9.]+)/gxms;
    croak "No resolution was parsed from pdb header remarks!" if ! $resValue;
    return $resValue;
}

sub _build_r_value {
    my $self      = shift;
    my $remarkStr = join("", map {${$_}} @{$self->remark_hash()->{3}});
    my ($rValue)
        = $remarkStr =~ /R \s VALUE \s+ \(WORKING \s (?: \+ \s TEST \s)* SET\) \s : \s
                         ([0-9.]+)/gxms;
    croak "No r-value was parsed from pdb header remarks!" if ! $rValue;
    return $rValue;
}

sub _build_experimental_method {
    my $self = shift;
    my ($expDataLine) = grep {/^EXPDTA/} @{$self->pdb_data()};
    my ($method) = $expDataLine =~ / \A EXPDTA \s+ (\S+.*?) \s* \z /gxms;
    croak "No experimental method was parsed from pdb header!" if ! $method;
    return $method;
}

sub _build_atom_serial_hash {
    my $self = shift;

    my %atomSerialHash = ();

    foreach my $atom (@{$self->atom_array()}) {
        $atomSerialHash{$atom->serial} = $atom;
    }

    return \%atomSerialHash;
}

sub _build_remark_hash {
    my $self = shift;
    
    my %remarks = ();
    
    for (my $i = 0 ; $i < @{$self->pdb_data()} ; ++$i){
        if ($self->pdb_data()->[$i] =~ /^REMARK \s* (\d+)/xms) {
            if (exists $remarks{$1}) {
                push(@{$remarks{$1}}, \$self->pdb_data()->[$i])
            }
            else {
                $remarks{$1} = [\$self->pdb_data->[$i]];
            }
        }
    }
    return \%remarks;
}

sub _build_atom_index {
    my $self = shift;

    my %hash = ();

    foreach my $atom (@{$self->atom_array()}) {
   
        my $chainID = $atom->chainID();
        my $resSeq  = $atom->resSeq();
        $resSeq .= $atom->iCode() if $atom->has_iCode();
        
        my $name    = $atom->name();

        if ( ! exists $hash{$chainID} )  {
            $hash{$chainID} = {};
        };

        if ( ! exists $hash{$chainID}->{$resSeq} ) {
            $hash{$chainID}->{$resSeq} = {};
        };

        $hash{$chainID}->{$resSeq}->{$name} = $atom;
    }
    return \%hash;    
}

sub _build_pdbresid_index {
    my $self = shift;

    my %pdbresid_index = ();
    
    foreach my $resid (keys %{$self->resid_index()}) {
        my $pdbresid = lc($self->pdb_code()) . $resid;
        my $value = $self->resid_index->{$resid};

        $pdbresid_index{$pdbresid} = $value;
    }
    return \%pdbresid_index;
}

sub _build_resid_index {
    my $self = shift;

    my %hash = ();

    foreach my $chain (keys %{ $self->atom_index } ) {
        my %chain_h = %{ $self->atom_index->{$chain} };
        
        foreach my $resSeq (keys %chain_h) {
            my $resid = "$chain.$resSeq";
            $hash{$resid} = $chain_h{$resSeq};
        }
    }

    croak "pdb: " . $self->pdb_code()
        . " - nothing indexed while attempting to index by resid"
            if ! %hash;

    return \%hash;
}

sub _build_missing_resid_index {
    my $self = shift;

    my %missing_resids = ();
    foreach my $chain (keys %{$self->missing_residues()}) {
        foreach my $resSeq (keys %{$self->missing_residues->{$chain}}){
            my $resid = "$chain.$resSeq";
            $missing_resids{$resid} = {};
        }
    }
    return \%missing_resids;
}

sub _build_terminal_atom_array {
    my $self = shift;

    my @terminalAtomArray = grep {$_->is_terminal()} @{$self->atom_array()};

    return \@terminalAtomArray;
}

sub _parse_remark465 {
    my $self = shift;
    
    # Return ref to empty hash if remark 465 d.n.e
    # (means that there are no missing residues) 
    return {}
        if ! exists $self->remark_hash->{465};

    my %resSeqs = ();

    my $headerFlag = 0;
    
    foreach my $lineRef (@{$self->remark_hash()->{465}}) {

        if (${$lineRef} =~ /REMARK 465   M RES C SSSEQI/){
            $headerFlag = 1;
            next;
        }
        
        # Skip lines until past header
        next if ! $headerFlag;
        
        my $line = ${$lineRef};

        my $resType = substr($line, 15, 3);
        my $chainID = substr($line, 19, 1);
        my $resSeq  = substr($line, 21, 6);

        $resSeq = rm_trail($resSeq);
        
        if (! exists $resSeqs{$chainID}) {
            $resSeqs{$chainID} = {$resSeq => $resType};
        }
        else {
            $resSeqs{$chainID}->{$resSeq} = $resType;
        }
    }
    
    return \%resSeqs;
}

sub _buildPDBCodeFromData {
    my $self = shift;

    my $headerLine;

    foreach my $line (@{$self->pdb_data}) {
        if ($line =~ /^HEADER/){
            $headerLine = $line;
            next;
        }
    }
    croak "No header line parsed from pdb data!"
        if ! defined $headerLine;

    return lc(substr($headerLine, 62, 4));
}

# This method returns a ref to a hash of form
# resID => [modifiedResidueName, standardResidueName]
# e.g. A.281 => [MSE, MET]
# This is used as a builder for attribute resid2ModResAref.
sub _parseMODRESLines {
    my $self = shift;

    # Parse MODRES lines from PDB data
    my @MODRESLines = grep {/^MODRES/} @{$self->pdb_data()};

    my %resid2ModResAref = ();
    
    foreach my $line (@MODRESLines) {
        # Use indices from PDB data format specfication to parse modified
        # residue information from line.
        # Use rm_trail to remove trailing whitespace
        my $chainID    = substr($line, 16, 1);
        my $resSeq     = rm_trail(substr($line, 18, 4));
        my $iCode      = rm_trail(substr($line, 22, 1));
        my $modResName = rm_trail(substr($line, 12, 3));
        my $stdResName = rm_trail(substr($line, 24, 3));

        my $resID = "$chainID.$resSeq";

        $resID .= $iCode if $iCode;
        
        $resid2ModResAref{$resID} = [$modResName, $stdResName];
    }
    
    return \%resid2ModResAref;
}

sub _build_resNameMod2StdHref {
    my $self = shift;

    my %modRes2StdResNames = ();

    # Each value of resid2ModResAref is a ref to an array of form
    # modifiedResidueName => standardResidueName
    # So loop through values and assign modRes => stdRes to hash
    foreach my $modAndStdAref (values %{$self->resid2ModResAref}) {
        my ($modRes, $stdRes) = @{$modAndStdAref};
        $modRes2StdResNames{$modRes} = $stdRes;
    }
    return \%modRes2StdResNames;
}

sub _buildHydroPhoHrefFromFile {
    my $self = shift;

    my $inFile = $self->threelc2hydrophobicValueFile();

    open(my $IN, "<", $inFile) or die "Cannot open file $inFile, $!";

    my %threelc2hydrophovalue = ();
    
    while (my $line = <$IN>) {
        next if $line =~ /^#/;

        my ($code, $value) = $line =~ /(\S+)\s+([0-9-.]+)/g;

        $threelc2hydrophovalue{uc $code} = $value;
    }
    return \%threelc2hydrophovalue;
}

sub _build_resID2secStructHref {
    my $self = shift;

    my %resID2secStruct
        = pdb::secstrCalculator->new(input => $self)->getResID2secStructHref();
    
    my %parsedResID2secStruct = ();
    while (my ($resID, $secStruct) = each %resID2secStruct) {
        if (exists $self->resid_index->{$resID}){
            $parsedResID2secStruct{$resID} = $secStruct;
        }
    }
    return \%parsedResID2secStruct;
}

### BUILD Method ###############################################################
################################################################################

# check if pdb is multi-model
sub BUILD {
    my $self = shift;

    return if ! $self->has_pdb_file();
     
    if ( $self->pdb_data && $self->_multi_model) {
        my $message
            = "pdb is a multi-model record: currently cannot handle these";
        
        my $error = local::error->new( message => $message,
                                       type => 'multi_model_pdb',
                                       data => {
                                           pdb_file => $self->pdb_data },
                                   );

        croak $error;
    }
}

sub _multi_model {
    my $self = shift;

     croak "Cannot determine multi-model status of pdb - no pdb data"
         if ! $self->has_pdb_data;

    if ( grep { m{ \A MODEL \s* \d+ \s* \z }xms } @{ $self->pdb_data() } ) {
        return 1;
    }
    else {
        return 0;
    }
};

### Methods ####################################################################
################################################################################

=item C<get_sequence(%args)>

Where %arg = (chain_id => ... , return_type => ... , include_missing => BOOL)

chain_id = Chain identifier of chain you wish to know the sequence of OR aref of
           chains identifiers.

return_type = 1 OR 3. (DEFAULT = 1)
    If 1, array of 1-letter AA codes is returned.
    If 3, array of 3-letter AA codes is returned.

including_missing = BOOL. (DEFAULT = TRUE)
    If TRUE, missing residues (normally parsed from pdb header data) are
    included in sequence.

std = BOOL. (DEFAULT = TRUE)
    If TRUE, modified residues are included using their standard name.
    e.g. MSE (selonmethionine) will be included as MET (methionine).

If a single chain identifier is passed, this method returns an array of 1 or
3-letter amino acid codes that represent the sequence of the chain specified.
e.g.

    @sequence = $pdb->get_sequence(chain_id => A, return_type => 3, include_missing => 1, std => 0);

If more than one chain identifier is passed, then a href is returned of form
chainID => sequenceAref.

If no chain identifiers are passed, then all chain sequences will be returned
(within a href as described above).

=cut

sub get_sequence {
    my $self = shift;

    # Currently haven't dealt with how to process multi-resName residues, so for
    # now throw an exception
    croak "pdb: " . $self->pdb_code()
        . " can't get_sequence for pdb containing multi-resName residues: "
        . Dumper $self->multi_resName_resids()
            if %{ $self->multi_resName_resids };
    
    my %arg = %{$self->_check_get_sequence_args(@_)};
    
    my %chainID2seqAref
        = map {$_ => $self->_get_chain_sequence($_, %arg)} @{$arg{chain_id}};

    if (keys %chainID2seqAref == 1) {
        # Get only key in the hash
        my $key = [keys %chainID2seqAref]->[0];
        return @{$chainID2seqAref{$key}};
    }
    else {
        return \%chainID2seqAref;
    }
}

sub _get_chain_sequence {
    my $self     = shift;
    my $chain_id = shift;
    my %arg      = @_;
    
    my %resSeq2ChainSeq
        = $self->map_resSeq2chainSeq(chain_id => $chain_id,
                                     include_missing => $arg{include_missing});
    
    # Sort resSeq by chainSeq to ensure corresponding residues are ordered
    # correctly
    my @sortedResSeqs
        = sort {$resSeq2ChainSeq{$a} <=> $resSeq2ChainSeq{$b}}
            keys %resSeq2ChainSeq;
    
    my @residues = ();
    
    # This hash is used to map modified residue names to standard residue names,
    # if arg std has been set to TRUE
    my %modRes2StdRes = %{$self->resNameMod2StdHref};
    
    # Find a residue name for each resSeq
    foreach my $resSeq (@sortedResSeqs){
        
        my $resName;
        
        # Value will eval to false if the chainID-resSeq is not present
        # (i.e. residue is missing)
        my $atom
            = eval {[values %{$self->atom_index->{$chain_id}->{$resSeq}}]->[0]};

        if($atom) {
            # Skip if atom from residue is solvent
            next if $atom->is_solvent();
            
            $resName = $atom->resName();
        }
        elsif ($arg{include_missing}) {
            # No atoms exists for this residue as it is missing, so use
            # missing_residues to find residue name
            $resName = $self->missing_residues->{$chain_id}->{$resSeq};
        }

        if ($arg{std} && exists $modRes2StdRes{$resName}) {
            # Map modified residue to standard residue name, if std arg
            # has been set to TRUE
            $resName = $modRes2StdRes{$resName};
        }
        
        push(@residues, $resName);
    }

    # Translate residue names to 1-lc if user has set return_type => 1
    if ($arg{return_type} == 1) {
        for my $i (0 .. @residues - 1) {
            my $onelc = eval { three2one_lc($residues[$i]) };
            if ($@) {
                # If three2one_lc has returned error then residue is a
                # non-standard AA.
                $onelc = 'X';
            }
            $residues[$i] = $onelc;
        }
    }
    return \@residues;
}

# Checks get_sequence args to make sure that a chain ID has been specified
# and that return_type value is valid (i.e. 1 or 3).
# Also supplies default values for other args. Defaults:
# return_type => 1, include_missing => 1, std => 1
sub _check_get_sequence_args {
    my $self = shift;
    my %arg = @_;

    $self->_process_get_sequence_chain_id_arg(\%arg);
    
    $arg{return_type} = 1 if ! exists $arg{return_type};
    $arg{include_missing} = 1 if ! exists $arg{include_missing};
    $arg{std} = 1 if ! exists $arg{std};

    croak "get_sequence: return_type '$arg{return_type}' is invalid."
        . " return_type must be 1 or 3."
            if $arg{return_type} != 1 && $arg{return_type} != 3;
    
    return \%arg;
}

sub _process_get_sequence_chain_id_arg {
    my $self    = shift;
    my $argHref = shift;

    # chain_id arg must either be an aref of chain_ids or a single chain_id
    # string which will be packaged into an aref. An exception will be thrown
    # if it anything else.
        
    if (! exists $argHref->{chain_id}) {
        # Set arg to aref of all chain ids present in pdb
        $argHref->{chain_id} = [keys %{$self->atom_index}];
    }
    elsif (ref $argHref->{chain_id} eq "") {
        # Arg is a string (hopefully corresponding to a chain ID!)
        # Package into an aref
        $argHref->{chain_id} = [$argHref->{chain_id}];
    }
    elsif (ref $argHref->{chain_id} ne 'ARRAY') {
        # Throw exception because arg is not valid form
        croak "Chain id arg " . $argHref->{chain_id} . "is not valid";
    }
}

=item C<getFASTAStr(%args)>

Where %args = (chain_id => "A", header => $headerStr, includeMissing => BOOL,)

chain_id = Chain identifier of chain you wish to obtain a FASTA String for.

header = header for FASTA String.
If no header is supplied, header if formed from pdb_code and supplied chain_id.

includeMissing = BOOL. (DEFAULT = 1)
If TRUE, residues missing in structure but present in pdb data are included in
FASTA sequence.

std = BOOL. (DEFAULT = TRUE)
If TRUE, modified residues are included using their standard name.
e.g. MSE (selonmethionine) will be included as MET (methionine).

This method returns a string containing a FASTA-formatted sequence for the given
chain. The option is given to include residues missing from the structure
(default is false) e.g.

    $pdb->getFASTAStr(chain_id => "A", header => "4houA_MUTANT", includeMissing => 1)

=cut

sub getFASTAStr {
    my $self = shift;
    my %arg = @_;

    # Throw error if no chain id was passed
    croak "No chain id was passed!" if ! exists $arg{chain_id};

    # Set default arg values if not supplied
    $arg{includeMissing} = 1 if ! exists $arg{includeMissing};
    $arg{std} = 1 if ! exists $arg{std};

    # Create default header if none was passed
    $arg{header} =  ">" . $self->pdb_code() . $arg{chain_id} . "\n"
        if ! exists $arg{header};

    # Get sequence in array of 1-letter codes
    my @seq = $self->get_sequence(chain_id => $arg{chain_id},
                                  return_type => 1,
                                  include_missing => $arg{includeMissing},
                                  std => $arg{std});
    
    my $seqStr = join("", @seq) . "\n";
   
    return $arg{header} . $seqStr;
}

=item C<seq_range_atoms(START, END)>

This method returns an array of atoms for residues within the given index range.

Range acts like perl indices i.e 0 = first residue, -1 = last residue, e.g.

    @first50ResidueAtoms = $pdb->seq_range_atoms(0, 49);

=cut

sub seq_range_atoms {
    my $self = shift;
    my ($start, $end) = @_;

    # Validate and inputs
    foreach my $input ($start, $end) {
        croak "seq_range_atoms: '$input' is not an int!"
            if ! is_int($input);
    }
    
    croak "seq_range_atoms: invalid range! $start - $end"
        if $start > $end && $end != -1;

    # Get splice of sorted_atom_arrays, then flatten arrays

    my @atom_arrays = @{$self->sorted_atom_arrays()};
    if ($end == -1) {
        $end =  $#atom_arrays;
    }
    my @atom_arrays_splice = @atom_arrays[$start .. $end];
    my @atoms = ();
    
    foreach my $atomARef (@atom_arrays_splice){
        push(@atoms, @{$atomARef});
    }
    return @atoms;
}

=item C<sorted_atom_arrays>

This method returns a ref to an array of refs, where each ref is an array of
residue atoms. Residue atom array refs are sorted by resid, e.g.

    $sortedAtomArraysAref = $pdb->sorted_atom_arrays();
    @firstResidueAtoms = @{$sortedAtomArraysAref->[0]};

=cut

sub sorted_atom_arrays {
    my $self = shift;

    my @atom_arrays = ();
    
    foreach my $chain_id (sort {$a cmp $b} keys %{$self->atom_index()}) {
        my $resSeqHref = $self->atom_index->{$chain_id};
        foreach my $resSeq (sort {pdb::pdbFunctions::compare_resSeqs($a,$b)}
                                keys % {$resSeqHref}){
            my @atoms = values %{$resSeqHref->{$resSeq}};

            my @sorted_atoms = sort_atoms(@atoms);
            push(@atom_arrays, \@sorted_atoms);
        }
    }
    return \@atom_arrays;
}

=item C<sort_atoms>

Given an array of atoms, returns an array of those atoms sorted by serial

=cut

sub sort_atoms {
    my @atoms = @_;
    
    @atoms = sort {$a->serial() <=> $b->serial()} @atoms;

    return @atoms;
}

=item C<read_ASA>

read_ASA assigns ASA values to the object's atoms and also assigns
the object's resid2RelASAHref attribute

Note that atom have three ASA attributes: ASAc, ASAm and ASAb

ASAc = Complex. This is the ASA value of an atom when the full structure is
       considered.

ASAm = Monomer. This is the ASA value of an atom when only the chain structure
       is considered, i.e. atoms that may be buried within complex interfaces
       are accessible to solvent in monomer form.

ASAb = SuBset. This is the ASA value of an atom when a subset of chains is
       considered. This is useful for cases where you want to treat multiple
       chains as one entity, e.g. a multi-chain antigen.

For a pdb object, atoms will be assigned ASAc values and relASA values will be
taken in context of full complex.

=cut

sub read_ASA {

    my $self = shift;
    my @errors = ();

    # solv must also be run in order to get relative ASAs
    my $solv = pdb::solv->new(input => $self);
    my $atomSerial2ASAHref = $solv->getOutput();
    
    my $ASAType = ref $self eq 'pdb' ? 'ASAc' : 'ASAb' ;
    
    # Use solv per atom output
    foreach my $atom (@{$self->atom_array()}) {
        next if $atom->is_solvent()
                || ! $atomSerial2ASAHref->{$atom->serial()};
        
        $atom->$ASAType($atomSerial2ASAHref->{$atom->serial()});
    }    
    $self->resid2RelASAHref($solv->resid2RelASAHref());    
    $self->has_read_ASA(1);
    
    return \@errors;
}

=item C<patch_centres(%args)>

where %args = (threshold => NUM, type => ... )

threshold = minimum rASA value of residue to be counted as a patch centre
            residue. DEFAULT = 25

type = ASA type (ASAm, ASAc or ASAb)

patch_centres returns a ref to an array of patch centre atoms.
Each patch centre atom is the most solvent accessible atom of a patch centre
residue. Each patch centre residue has an  rASA value over the given threshold.

patch_centres also returns a ref to an array of errors. It is important to check
this array.

Example:
    ($errorAref, $pCentreAref) = $pdb->patch_centres(threshold => 20, type => 'ASAc');
    if (@{$errorAref}) {
        croak "There were errors when patch_centres ran!" . Dumper $errorAref;
    }
    else {
        print "Here are my patch centres", @{$pCentreAref};
    }


=cut

sub patch_centres {
    my $self = shift;
    
    my %arg = @_;

    my @errors = ();
 
    if (! %{$self->resid2RelASAHref()}) {
        my $message
            = "read_ASA must be run before patches_centres "
                . "to supply a resid2RelASAHref";
        my $error = local::error->new(
            message => $message,
            type => 'no resid2RelASAHref',
            data => {},
        );
        push(@errors, $error);

        return \@errors;
    }
    
    my $threshold = exists $arg{threshold} ? $arg{threshold} : 25;

    # Get ASA type to decide on patch central atom from each path centre
    # If type has been supplied, use supplied type - else, use ASAm if self is a
    # chain, otherwise use ASAc
    my $type
        = exists $arg{type} ? $arg{type}
        : ref $self eq 'pdb' ? 'ASAc'
        : 'ASAb';
    
    my @central_atoms = ();
    
    foreach my $chain ( keys %{ $self->atom_index } ) {
        my %chain_h = %{ $self->atom_index->{$chain} };
        
        foreach my $resSeq ( keys %chain_h ) {
            
            my %atom_h = %{ $chain_h{$resSeq} };
                       
            my $CA_flag = 0;
            
            foreach my $atom_name (keys %atom_h) {
                $CA_flag = 1 if $atom_name eq 'CA';
                
                my $atom = $atom_h{$atom_name};

                next if $atom->is_solvent();                
            }

            # Patch centre residues must have a CA atom to be run through
            # makepatch
            next if ! $CA_flag;

            # Test to see if residue rASA is over threshold
            my $resid = $chain . "." . $resSeq;
            my $relASA =  $self->resid2RelASAHref->{$resid}->{allAtoms};

            next if $relASA < $threshold;

            my $central_atom = $self->highestASA($resid, $type);
            push(@central_atoms, $central_atom);
        }
    }  
    return(\@errors, \@central_atoms);
}

sub _is_patch_centre {

    my $self = shift;
    my $threshold = shift;
    my $attribute = shift;
    my @atoms = @_;

    my @nonHydAtoms = grep {$_->element() ne 'H'} @atoms;

    foreach my $atom (@nonHydAtoms) {
        my $predicate = "has_$attribute";
        
            if( ! $atom->$predicate() ){
                                
                my $message
                    =  "Could not determine patch centre status because "
                      . "atom " . $atom->serial
                      . " has no $attribute value";
                
                my $error = local::error->new(
                    message => $message,
                    type => "no_atom_$attribute" . "_set",
                    data => { atom => $atom, },
                );

                return $error;
            }        
    }
    
    @nonHydAtoms
        = sort {$b->$attribute() <=> $a->$attribute()} @nonHydAtoms;

    my $total = 0;

    map { $total += $_->$attribute() } @nonHydAtoms;

    if ($total >= $threshold) {
        return $nonHydAtoms[0];
    }

    return -1;
}

=item C<highestASA(RESID, ASATYPE)>

RESID = residue identifier (e.g. A.213)
ASATYPE = ASA value type (ASAc, ASAm or ASAb)

This method returns the atom from the specified residue with the highest
ASA value, of type specified.

e.g.

    $highestASAAtom = $pdb->highestASA("A.213", "ASAc");

=cut

sub highestASA {
    my $self  = shift;
    my $resid = shift or croak "highestASA must be passed a resid";

    my $ASA_type
        =  $_[0] ? $_[0] :
           ref $self eq 'pdb'   ? 'ASAc'
         : ref $self eq 'chain' ? 'ASAm'
         : '' ;

    croak "highestASA: something went wrong assigning ASA type"
        if ! $ASA_type;
    
    croak "resid '$resid' was not found in resid index"
        if ! exists $self->resid_index->{$resid};

    my @atoms = values %{ $self->resid_index->{$resid} };

    # If element is defined, use to avoid hydrogen
    @atoms = grep {! defined $_->element || $_->element() ne 'H'} @atoms;
    
    foreach my $atom (@atoms) {        
        my $predicate = "has_$ASA_type";
        croak 'atom ' . $atom->serial() . " has no $ASA_type value\n"
            . $atom . "\n"
                if ! $atom->$predicate();
    }

    my @sorted = sort { $b->$ASA_type <=> $a->$ASA_type  } @atoms;

    my $top_ASA_atom = $sorted[0];
    
    return $top_ASA_atom;
}

=item C<map_resSeq2chainSeq(chain_id => CHAIN_ID, include_missing => BOOL)>

CHAIN_ID = chain identifier of chosen chain.
include_missing = include missing residues. DEFAULT = 1

Maps resSeq numbers to chainSeq count numbers (equivalent to pdbcount num
in pdbsws.)
Returns hash with form: resSeq => chainSeq

=cut

sub map_resSeq2chainSeq {
    my $self = shift;

    my %args = @_;

    my $chain_id = $args{chain_id};
    my $include_missing = exists $args{include_missing} ? $args{include_missing} : 1;
        
    croak "pdb: " . $self->pdb_code() . " no residues found for chain "
        . " $chain_id" if ! exists $self->atom_index->{$chain_id};
    
    my $chainSeq = 0;
    my %return_map = ();

    my $prevResID = "";

    my @orderedResSeqs = ();

    # Get resSeqs as ordered in atom array
    foreach my $atom (@{$self->atom_array()}) {
        if ($atom->resid() ne $prevResID && $atom->chainID() eq $chain_id
                && ! $atom->is_solvent()) {
            push(@orderedResSeqs, [split(/\./, $atom->resid())]->[1]);
            $prevResID = $atom->resid();
        }
    }

    my %rS2Ord;
    @rS2Ord{@orderedResSeqs} = (1 .. scalar @orderedResSeqs);

    if ($include_missing) {
        # Get missing resSeqs from chain
        my @missingResSeqs = keys %{$self->missing_residues()->{$chain_id}};

        my @sortedResSeqs
            = sort {pdb::pdbFunctions::compare_resSeqs($a, $b, \%rS2Ord)}
                (@orderedResSeqs, @missingResSeqs);
        
        @rS2Ord{@sortedResSeqs} = (1 .. scalar @sortedResSeqs);

    }
    return %rS2Ord;
}

=item C<map_chainSeq2resSeq(CHAIN_ID)>

CHAIN_ID = chain identifier of chosen chain.

Returns hash chainSeq => resSeq, i.e. the reverse of map_resSeq2chainSeq

=cut
    
sub map_chainSeq2resSeq {
    my $self = shift;
    
    my $chain_id
        = shift or croak "map_resSeq2chainSeq must be passed a chain id";
    
    croak "pdb: " . $self->pdb_code() . " no residues found for chain "
        . " $chain_id" if ! exists $self->atom_index->{$chain_id};
    
    my %resSeq2chainSeq = $self->map_resSeq2chainSeq(chain_id => $chain_id);
    
    my %return_map
        = map { $resSeq2chainSeq{$_} => $_ } keys %resSeq2chainSeq;
    
    return %return_map;
}

=item C<create_chains(CHAIN_ID_ARRAY)>    

CHAIN_ID_ARRAY = array of chain identifiers

Method to create chain objects from pdb object. Returns array of chains.
An array of chain ids can be passed to the method; if this is the case then
the returned chains will be in the order specified in the passed array.

If CHAIN_ID_ARRAY contains a hyphen, then all chains present in the pdb
that do have chain ids corresponding to the chain ids following the hyphen
will be returned. e.g,

    ('-', 'A')

will return all chains except A.

=cut
    
sub create_chains {
    my $self = shift;
    my @passed_chain_ids = @_;
    my @return_chain_ids = ();
    my @except_chain_ids = ();
    
    # If chain ids have been passed, check that all chain ids are found in
    # this pdb
    if (@passed_chain_ids) {
        my %chain_ids = map { $_ => 1 } @{$self->get_chain_ids()};
        while (@passed_chain_ids) {
            my $chain_id = shift @passed_chain_ids;
            if ($chain_id eq '-') {
                @except_chain_ids = @passed_chain_ids;
                last;
            }
            else {
                croak "Passed chain_id $chain_id was not found in pdb!"
                    if ! exists $chain_ids{$chain_id};
                push(@return_chain_ids, $chain_id);
            }
        }
    }
        
    my %atoms = map { $_ => [] } @{$self->get_chain_ids()};

    # Hash atoms by chain id
    foreach my $atom (@{$self->atom_array()}) {
        push(@{$atoms{$atom->chainID()}}, $atom);        
    }

    my %chains = ();

    # Create chains
    foreach my $chain_id (keys %atoms) {
        my $chain = chain->new(pdb_code => $self->pdb_code(),
                               chain_id => $chain_id,
                               atom_array => $atoms{$chain_id},
                               pdb_data => $self->pdb_data());

        $chains{$chain->chain_id()} = $chain;
    }
        
    my @return_chains = ();
    
    # If chain_ids have been specified, return chains in specified order
    if (@return_chain_ids) {
        foreach my $chain_id (@return_chain_ids) {
            push(@return_chains, $chains{$chain_id});
        }
    }
    if (@except_chain_ids) {
        foreach my $chain_id (keys %chains) {
            push(@return_chains, $chains{$chain_id})
                if ! grep {$_ eq $chain_id} @except_chain_ids;
        }
    }
    return @return_chains ? @return_chains : values %chains;
}

=item C<get_chain_ids>

Returns ref to array containing chainIDs found in pdb

=cut

sub get_chain_ids {
    my $self = shift;
    my @chain_ids = keys (%{$self->atom_index()});
    return \@chain_ids;
}

=item C<getResIDs>

Returns an array of all the resIDs of the chain
    include_missing : include resIDs of missing residues. DEFAULT = TRUE

=cut

sub getResIDs {
    my $self = shift;
    my %opt  = (@_);

    my @resIDs = keys %{$self->resid_index}; 
    my $includeMissing
        = exists $opt{include_missing} && ! $opt{include_missing} ? 0
        : 1;
        
    if (! $includeMissing) {
        @resIDs = grep {! exists $self->missing_resid_index->{$_}} @resIDs;
    }
    return @resIDs;
}

=item C<getAbPairs(CHAIN_ARRAY)>

Determines antibody pairs found within pdb.
INPUT:
This fuction can optionally be passed an array of pdb chains.
e.g.
    $pdb->getAbPairs()
    $pdb->getAbPairs(@chains)
    $pdb->getAbPairs($pdb->create_chains())
OUTPUT:
     Returns three refs to arrays of:
      1. Antibody Pairs:  [ [$heavyChain, $lightChain, $numContacts], ... ]
      2. Unpaired Chains: [ $chainA, $chainB, ... ]
      3. scFv Chains:     [ $chainScFv, ... ]

=cut

sub getAbPairs {
    my $self = shift;
    my @chains = @_;

    # Create chains if none have been passed
    @chains = $self->create_chains() if ! @chains;

    my %chainTypes = (Heavy => [], Light => [],
                      antigen => [], scFv => []);

    # Hash to keep track of which chains have been paired
    my %paired = ();

    # Hash chains by chain type
    _hashChains(\@chains, \%chainTypes, \%paired);

    # Get inter-chain contacts
    my $cContacts = pdb::chaincontacts->new(input => $self);
    my $cResult = $cContacts->getOutput();
    
    my @heavyLightCombs
        = _heavyLightCombinations($chainTypes{Heavy}, $chainTypes{Light},
                                  $cResult);
    
    # Sort heavy-light chain combinations by number of inter-chain contacts
    @heavyLightCombs = sort { $b->[2] <=> $a->[2] } @heavyLightCombs;

    my @finalCombinations = _getFinalCombinations(\@heavyLightCombs, \%paired);

    my @unpaired = ();

    # Get those chains that were not paired
    foreach my $chain ( @{$chainTypes{Heavy}}, @{$chainTypes{Light}} ) {
        if (! $paired{$chain->chain_id()}) {
            push(@unpaired, $chain);
        }
    }
    
    # Return refs to arrays of:
    #  final combinations, non-paired ab chains and  scFv chains
    return(\@finalCombinations, \@unpaired, $chainTypes{scFv});
}

# This subroutine processes an array of heavy-light chain combinations to find
# the set of pairs where each chain is only paired once and the sum of all
# contacts is the highest possible
sub _getFinalCombinations {
    my($combAref, $pairedHref) = @_;

    my @finalCombinations = ();
    
    foreach my $combination (@{$combAref}) {
        my $heavy = $combination->[0]->chain_id();
        my $light = $combination->[1]->chain_id();
        # Check that neither heavy nor light has been paired already
        unless ($pairedHref->{$heavy} || $pairedHref->{$light}) {
            # Set chains to paired
            $pairedHref->{$heavy} = 1;
            $pairedHref->{$light} = 1;
            
            push(@finalCombinations, $combination);
        }
    }
    return @finalCombinations;
}

# This subroutine hashes chain by abVariable type, into two hashes, in
# preparation for pairing.
# Input: refs to array of chains, hash for chain types,
# hash for heavy/light chains (to track pairing)
sub _hashChains {
    my($chainsAref, $chainTypesHref, $pairedHref) = @_;
    
    foreach my $chain (@{$chainsAref}) {
        my $chainType;
        my $ret = eval {$chainType = $chain->isAbVariable(); 1};
        
        if (! $ret) {
            print {*STDERR} "ERROR!: " . Dumper $@;
            croak $@;
        }
            
        if (! $chainType) {
            $chainType = 'antigen';
        }
        elsif ($chainType eq 'Heavy' || $chainType eq 'Light') {
            $pairedHref->{$chain->chain_id()} = 0;
        }

        # Set to is_ab_variable to ab variable chain type, unless chain is
        # antigen.
        $chain->is_ab_variable($chainType eq 'antigen' ? 0 : $chainType);
        
        push(@{$chainTypesHref->{$chainType}}, $chain);
    }
}

# Given references to arrays of heavy and light chains and a chaincontacts
# result for the pdb, returns an array of arrays with form:
# ( [heavyChain, lightChain, numContacts], ... )
sub _heavyLightCombinations {
    my($heavyAref, $lightAref, $cResult) = @_;

    my @heavyLightCombinations = ();
    
    # Get all combinations of heavy and light chains
    foreach my $heavy (@{$heavyAref}) {
        foreach my $light (@{$lightAref}) {
            # Find number of contacts between heavy and light chain
            my $residAref = $cResult->chain2chainContacts([$heavy], [$light]);
            my $numContacts = scalar @{$residAref};

            my $combinationAref = [$heavy, $light, $numContacts];
          
            push(@heavyLightCombinations, $combinationAref);
        }
    }
    return @heavyLightCombinations;
}

=item C<rotate2PCAs(ATOM_SUBSET)>

This method rotates atoms of the pdb so that the PC1 and PC2 of the atoms
lay on the x and y axes respectively. An array of atoms or resids can be
passed; in this case, the PC1 and PC2 of this subset of atoms (or atoms from
resids) can be used. 

=cut

sub rotate2PCAs {
    my $self = shift;
    my @centering = @_;

    my @centralAtoms = $self->_validateRotate2PCAsArgs(@centering);

    # Centre pdb using all atoms if no centre atoms have been passed
    if (! @centralAtoms) {
        @centralAtoms = @{$self->atom_array()};
    }
    
    # Find centre point of centralAtoms and move this point to origin
    my($cx, $cy, $cz) = pdb::pdbFunctions::findAtomMean(@centralAtoms);

    # Translate all so that centre point is at origin (0, 0, 0)
    foreach my $atom (@{$self->atom_array()}) {
        $atom->x($atom->x() - $cx);
        $atom->y($atom->y() - $cy);
        $atom->z($atom->z() - $cz); 
    }

    # Get rot matrix for transforming x and y unit vectors to
    # centralAtoms PC1 and PC2
    my @vectors = map { vector($_->x, $_->y, $_->z) }  @centralAtoms;
    my $rotationMatrix = rotate2pc(@vectors);

    # Rotate all points using matrix
    pdb::pdbFunctions::rotateAtoms($rotationMatrix, $self->atom_array());

    return 1;
}

# This sub obtains an array of atoms from rotate2PCAs args
sub _validateRotate2PCAsArgs {
    my $self = shift;
    
    my @args = @_;
    my @atoms = ();
    
    foreach my $ele (@args) {
        if (ref $ele eq 'atom') {
            push(@atoms, $ele);
        }
        # Is element a resid?
        elsif (exists $self->resid_index->{$ele}) {
            push(@atoms, values %{$self->resid_index->{$ele}});
        }
        else {
            croak "Invalid arg '$ele' passed to rotate2PCAs: "
                . "args must be atoms or valid resids";
        }
    }
    return @atoms;
}

=item C<rotate2Face>

This method checks the number of +z and -z atom co-ordinates and flips the pdb
180degrees around the z axis so that the z axis points through the "body" of
the pdb; i.e. flips the pdb so that you are NOT looking through the body of
the pdb, if you are looking down the z axis

=cut

sub rotate2Face {
    my $self = shift;

    my $zPlusCount = 0;
    my $zMinusCount = 0;

    # Get counts of atoms that have positive or minus z co-ordinates
    foreach my $atom (@{$self->atom_array()}) {
        abs $atom->z() == $atom->z() ? ++$zPlusCount : ++$zMinusCount;
    }

    if ($zPlusCount > $zMinusCount) {
        # Rotate pdb 180 degrees around Z axis

        # Get rotation matrix
        my $RM = Math::MatrixReal->new_from_rows(
            [ [cos(pi), -sin(pi), 0],
              [sin(pi), cos(pi), 0],
              [0, 0, 1] ]
        );

        # Rotate all atoms
        pdb::pdbFunctions::rotateAtoms($RM, $self->atom_array());
    }
    return 1;
}

=item C<clearEpitopeLabels>

This method sets all atom->is_epitope() labels of pdb atoms to 0

=cut

sub clearEpitopeLabels {
    my $self = shift;

    map { $_->is_epitope(0) } @{$self->atom_array()};
}

=item C<storeInFile>

This method saves the pdb in a binary file. A filename can be passed,
otherwise the pdb code is used. Returns file name of binary file
e.g
    $pdb->store();
    $storedPDBFname = $pdb->store();
    $pdb->store("myStoredPDB.obj");

=cut

sub storeInFile {
    my $self = shift;
    my ($fName) = @_;

    $fName = $self->pdb_code() . ".obj" if ! $fName;

    store $self, $fName;

    return $fName;
}

=item C<squaredDistance(ATOM1, ATOM2)>

This method returns the squared distance between the two atoms passed.
Only use this method if speed is not a concern.

=cut

sub squaredDistance {
    my $self  = shift;
    my $atom1 = shift;
    my $atom2 = shift;

    croak "squaredDistance must be passed two atoms" if ! ($atom1 && $atom2);

    my @deltas = map {($atom1->$_ - $atom2->$_)**2} qw(x y z);

    my $sum = 0;
    $sum += $_ foreach @deltas;

    return $sum;
}

sub calcAverageHydrophobicity {
    my $self = shift;

    my $threelc2valueHref = $self->threelc2hydrophobicValueHref;

    return $self->calcAveragePropensity($threelc2valueHref);
}

sub calcAveragePropensity {
    my $self = shift;
    my $threelc2valueHref = shift;
    
    # Get residues using get_sequence
    my @getSequenceRet = $self->get_sequence(return_type => 3);
    
        my @sequence = ();
    
    if (ref $getSequenceRet[0] eq 'HASH') {
        # Flatten sequence arrays (referenced as values in returned hash)
        # from each chain
        @sequence = map {@{$_}} values %{$getSequenceRet[0]};
    }
    else {
        # Sequence has been returned as array
            @sequence = @getSequenceRet;
        }
    
    my $total;
        $total += $threelc2valueHref->{$_} foreach @sequence;
    my $avg = $total / scalar @sequence;
    
    return $avg;
}

sub planarity {
    my $self = shift;
    my @points = map {vector($_->x, $_->y, $_->z)} @{$self->atom_array()}; 
    return get_rms_difference_of_points_from_bestfit_plane(@points);
}

# pp = protein-protein
sub labelppHbondedAtoms {
    my $self = shift;

    my @ppHbs
        = pdb::hbondFinder->new(input => $self, type => 'pp')->getHbonds();

    foreach my $hb (@ppHbs) {
        my $dSerial = $hb->donorSerial();
        my $aSerial = $hb->acceptorSerial();
        if (exists $self->atom_serial_hash->{$dSerial}
                && exists $self->atom_serial_hash->{$aSerial}) {
            
            my $dAtom = $self->atom_serial_hash->{$dSerial};
            my $aAtom = $self->atom_serial_hash->{$aSerial};
            $dAtom->HbAcceptor($aAtom);
            $aAtom->HbDonor($dAtom);            
        }
    }
}

sub labelSSbondedAtoms {
    my $self = shift;

    my @ssBonds = pdb::ssFinder->new(input => $self)->getssArray($self);

    foreach my $ssBond (@ssBonds) {
        my ($aSerial, $bSerial) = @{$ssBond->atomSerialPairAref};
        
        if (exists $self->atom_serial_hash->{$aSerial}
                && exists $self->atom_serial_hash->{$bSerial}) {
            my $aAtom = $self->atom_serial_hash->{$aSerial};
            my $bAtom = $self->atom_serial_hash->{$bSerial};

            $aAtom->SSbond($bAtom);
            $bAtom->SSbond($aAtom);
        }
    }
}


=item C<getResNames(resIDS)>

Given an array of residue IDs (e.g, A.101, B.102), this method
returns the corresponding residue names (e.g, HIS, GLY)

=cut

sub getResNames {
    my $self = shift;
    my @resIDs = @_;

    # For each resID, get first atom from resid_index then get that atom's
    # name
    return map { [values %{$self->resid_index->{$_}}]->[0]->resName() } @resIDs;
}

sub readAtomRadii {
    my $self = shift;
    map {$_->readRadius($self->radiusFinder)} @{$self->atom_array};
    $self->atomRadiiAreAssigned(1);
}

__PACKAGE__->meta->make_immutable;

################################ END OF PACKAGE ################################
################################################################################

package chain;


=head1 CLASS

pdb - A class for creating chain objects from pdb files.

=cut

=head1 SYNOPSIS

use pdb::chain;
$chainObject = chain->new(pdb_code => '1adjs',
                          chain_id => 'A',
                          pdb_file => '1adjs.pdb');

# Or ...
$pdbObject = pdb->new(pdb_code => '1djs', chain_id =>'A'); # Grab file automatically!

=cut


=head1 DESCRIPTION

The chain class is used for performing functions on single chains, rather than
an entire pdb, e.g.
     is a chain an antibody variable chain?
     what residues form an interface between this chain and another chain?

chain is an extension of pdb, so all pdb methods can be used with chain objects,
e.g. atom_array, resid_index, etc.

=cut


use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;
use Carp;
use TryCatch;

use pdb::pdbsws;
use pdb::idabchain;
use pdb::kabatnum;
use pdb::chaincontacts;

use TCNUtil::write2tmp;

extends 'pdb';

=head1 Methods

=over 12

=cut

### Attributes #################################################################
################################################################################

=item C<chain_id>

Chain identifier of chain

=cut

has 'chain_id' => (
    isa => 'ValidChar',
    is => 'rw',
    required => 1,
);

has 'pdbID' => (
    is => 'ro',
    lazy => 1,
    default => sub {$_[0]->pdb_code() . $_[0]->chain_id()},
);


=item C<accession_codes>

Returns a ref to an array of accession codes that the chain object is assigned
to, via pdbsws.

=cut

has 'accession_codes' => (
    isa => 'ArrayRef[ValidAC]',
    is => 'ro',
    lazy => 1,
    builder => '_get_acs',
);

=item C<chain_length>

returns length of chain sequence

=cut

has 'chain_length' => (
    isa => 'Int',
    is => 'ro',
    lazy => 1,
    builder => '_build_chain_length',
);

=item C<is_ab_variable>

Returns TRUE if chain can be identified as an antibody variable region chain.

=cut

has 'is_ab_variable' => (
    isa => 'Str',
    is => 'rw',
    lazy => 1,
    builder => 'isAbVariable',
);

=item C<is_het_chain>

Returns TRUE if chain consists entirely of HETATMS.

=cut

has 'is_het_chain' => (
    isa => 'Bool',
    is => 'ro',
    lazy => 1,
    builder => '_build_is_het_chain',
);

=item C<is_nt_chain>

Returns TRUE if chain consists entirely of nucleotide bases.

=cut

has 'is_nt_chain' => (
    isa => 'Bool',
    is => 'ro',
    lazy => 1,
    builder => '_build_is_nt_chain',
);

=item C<cluster_id>

This is assigned when cdhit is run with chain as input.

=cut

has 'cluster_id' => (
    is => 'rw',
    isa => 'Num',
);

### Attribute Builder Methods ##################################################
################################################################################

# Returns number of non-solvent residues found in chain
sub _build_chain_length {
    my $self = shift;

    my %resid_h = %{ $self->resid_index };

    my $count = 0;

    foreach my $resid (keys %resid_h) {
        my @atoms = values %{ $resid_h{$resid} };
        
        ++$count if ! $atoms[0]->is_solvent();
    }
    return $count;
}

sub _get_acs {
    my $self = shift;

    my $pdbsws = pdb::pdbsws->new();

    my $pdbid = $self->pdb_code . $self->chain_id;

    return [ $pdbsws->get_ac($pdbid) ];
}

sub _build_is_het_chain {
    my $self = shift;

    foreach my $atom (@{$self->atom_array()}) {
        return 0 if ! $atom->is_het_atom;
    }

    return 1;
}

sub _build_is_nt_chain {
    my $self = shift;

    my %nts = map {$_ => 1} qw(A D G C U);

    my $flag = 1;
    
    foreach my $atom (@{$self->atom_array()}) {
        if (! exists $nts{$atom->resName()}) {
            $flag = 0;
            last;
        }
    }

    return $flag;
}


### Around modifiers ###########################################################
################################################################################

around '_parse_remark465' => sub {
    my $orig = shift;
    my $self = shift;

    my $chain_id = $self->chain_id();

    my $origOutput =  $self->$orig;
    if (exists $origOutput->{$chain_id}) {
        my $chainOnly = {$chain_id => $origOutput->{$chain_id}};
        return $chainOnly;
    }
    else {
        return {};
    }
    
};

# Automatically set arg chain_id
around [qw(get_sequence getFASTAStr map_resSeq2chainSeq _parse_ATOM_lines)] => sub {

    my $orig = shift;
    my $self = shift;

    my %arg = @_;

    $arg{chain_id} = $self->chain_id() if ! exists $arg{chain_id};
    
    return $self->$orig(%arg);
    
};

# Automatically send chain id
around [qw(map_chainSeq2resSeq)] => sub {
    my $orig = shift;
    my $self = shift;

    my @arg = ($self->chain_id());
    
    return $self->$orig(@arg);    
};

### Methods ####################################################################
################################################################################

# This method checks if the chain is an antibody variable chain, or antigen
# Is the chain is antibody variable, the type is returned i.e. Light or Heavy
sub isAbVariable {
    my $self = shift;

    my $idabchain = pdb::idabchain->new(input => $self);
    my $chainType = "";
    
    try {
        $chainType = $idabchain->chainIs();
    }
    catch ($err where {ref $_ eq 'local::error'}) {
        if ($err->type() eq 'AllHETATMInputFile'){
            # Chain conists of HETATMs only, thus chain is not AbVariable
            $chainType = 'Antigen';
        }
        elsif ($err->type() eq 'NoOutputForChainID'){
            if ($self->is_het_chain()) {
                # Chain consists of HETATMS only, thus chain is not AbVar
                $chainType = 'Antigen';
            }
            else {
                croak $err;
            }
        }
        elsif ($err->message() =~ /Can.*'t split sequence/) {
            if ($self->is_nt_chain()) {
                # Chain consists of nucleotide atoms only, thus chain is not
                # AbVar
                $chainType = 'Antigen';
            }
        }
        else {
            croak $err;
        }
    };
    
    if ($chainType eq 'Antigen') {
        return 0;
    }
    else {
        return $chainType;
    }
}

=item C<kabatSequence>

This method runs kabatnum using pdb::kabatnum and assigns CDR and kabatSeq
attributes to atoms of the chain

=cut

sub kabatSequence {
    my $self = shift;

    my $kabatnum = pdb::kabatnum->new(input => $self);

    $kabatnum->sequenceChain();
}

=item C<determineEpitope(ABCHAINS, DIST1, DIST2)>

ABCHAINS = ref to an array of antibody chains.
DIST1 = primary contact distance paramter. DEFAULT = 4 (angstrom)
DIST2 = secondary contact distance parameter DEFAULT = 4 (angstrom)

This method determines the atoms that are within a given distance threshold
(default 4A) of the CDRs of the given antibody chains. The antibody chain
CDR atoms must be labelled as such (i.e. $atom->is_CDR() == 1).

Non-CDR contacting residues can be included by passing an additional, larger
distance parameter. This distance parameter is used to find antigen residues
close to the CDRs that are also in contact with any antibody residue. The
larger the second distance threshold, the larger the epitope.
  e.g. $chain->determineEpitope(\abChains, 4, 8);

If DIST1 = DIST2 then no non-CDR contacting peripheral residues will be
included.

=cut

sub determineEpitope {
    my($self, $abChainsAref, $normDistance, $exDistance) = @_;

    # Default distances if not supplied by user.
    $normDistance = 4 if ! $normDistance;
    $exDistance = 4 if ! $exDistance;
    
    ## STEP 1: CDR Atoms - Ag Chain Contacts, EXTENDED Distance 
    # Get array of CDR atoms from antibody chains
    my @CDRAtoms = ();
    foreach my $chain (@{$abChainsAref}) {
        push(@CDRAtoms, $chain->getCDRAtoms());
    }

    croak "No CDR atoms were found in passed antibody chains!"
        if ! @CDRAtoms;
    
    # Create array of CDR atoms and antigen chain atoms
    my @CDRAndAgAtoms = (@CDRAtoms, @{$self->atom_array()});

    # Test for contacts between CDRs and antigen
    my $chainContacts = pdb::chaincontacts->new(input => \@CDRAndAgAtoms,
                                                threshold => $exDistance);
    my $exContactResult = $chainContacts->getOutput();

    # STEP 2: All Ab Atoms - Ag Chain Contacts, NORMAL Distance
    
    # Create array of all atoms (all ab + antigen)
    my $allAtomsAref
        = pdb::pdbFunctions::generateAtomAref($self, @{$abChainsAref});

    # Test for contacts between whole ab chains and antigen, NORMAL distance
    $chainContacts->input($allAtomsAref);
    $chainContacts->threshold($normDistance);
    
    my $allcContactResult = $chainContacts->getOutput();

    # Get array of resids that are contacting ab chains and in proximity to CDRs
    my $eptiopeResidAref = _determineEpitopeResids($abChainsAref, $self,
                                                   $exContactResult,
                                                   $allcContactResult);

    # STEP 3: Filter these residues so that only those residues that are also
    # labelled interface are included
    my @interfaceResidues = $self->getInterfaceResidues($abChainsAref);

    my %epitopeResids = map {$_ => 1} @{$eptiopeResidAref};

    my @epitopeAndInterface
        = grep {exists $epitopeResids{$_}} @interfaceResidues;
    
    # Label atoms from resids as epitope
    $self->labelEpitopeAtoms(@epitopeAndInterface);
}

# This subroutine determines the residues that are part of the epitope, from the
# two chaincontact::result objects that are generated within determineEpitope
# INPUT: AbChainsAref, AgChain, ExtendedContactResult, AllContactResult
sub _determineEpitopeResids {
    my($abChainsAref, $AgChain, $eCR, $aCR) = @_;

    my @residHrefs = ();
    
    foreach my $result ($eCR, $aCR) {
        my $residAref = $result->chain2chainContacts($abChainsAref, [$AgChain]);

        my %resids = map { $_ => 1 } @{$residAref};
        
        push(@residHrefs, \%resids);
    }

    # Keep any resid that is within EXTENDED distance of the CDRs AND normal
    # distance of any ab residue
    my @epitopeResids
        = grep ($residHrefs[1]->{$_}, keys %{$residHrefs[0]});

    return \@epitopeResids;
}

=item C<labelEpitopeAtoms(RESIDS)>

RESIDS = array of resids (example resid = A.123)

Labels atoms of chain as epitope according to the array of resids passed to it.

=cut

sub labelEpitopeAtoms {
    my $self = shift;
    my @resids = @_;

    foreach my $resid (@resids) {

        foreach my $atom (values %{$self->resid_index->{$resid}}) {
            $atom->is_epitope(1);
        }
    }
}

=item C<labelInterfaceAtoms(RESIDS)>

RESIDS = array of resids (example resid = A.123)

Labels atoms of chain as epitope according to the array of resids passed to it.

=cut

sub labelInterfaceAtoms {
    my $self = shift;
    my @resids = @_;

    foreach my $resid (@resids) {
        foreach my $atom (values %{$self->resid_index->{$resid}}) {
            $atom->is_interface(1);
        }
    }
}


=item C<getResSeqs(include_mising => BOOL)>

Returns an array of all the resSeqs of the chain

Opts:
    include_missing : include resSeqs of missing residues. DEFAULT = TRUE

=cut

sub getResSeqs {
    my $self = shift;
    my %opt  = (@_);

    my @resSeqs = keys %{$self->atom_index->{$self->chain_id}};

    my $includeMissing
        = exists $opt{include_missing} && ! $opt{include_missing} ? 0
            : 1;
    
    if (! $includeMissing) {
        my @resIDs = map {$self->chain_id() . ".$_"} @resSeqs;
        my @presentResIDs
            = grep {! exists $self->missing_resid_index->{$_}} @resIDs;
        @resSeqs = map { [split(/\./, $_)]->[1] } @presentResIDs;
    }
    return @resSeqs;
}

=item C<getEpitopeResSeqs>

Returns array of residue resSeqs, where each residue has at least one atom
labelled as epitope

=cut

sub getEpitopeResSeqs {
    my $self = shift;

    my %epitopeResSeqs = ();
    
    foreach my $atom (@{$self->atom_array()}) {
        $epitopeResSeqs{$atom->resSeq()} = 1
            if $atom->is_epitope();
    }
    return keys %epitopeResSeqs;
}

=item C<getInterfaceResidues(CHAINS or IN COMPLEX RESIDUE ASAs)>

CHAINS = ref to array of chains that the chain forms an interface with.
IN COMPLEX RESIDUE ASAs = ref to hash of form ResID => relASA
 e.g A.139 => 0.19
relASA measures should be for the resiudes in their complexed state.
This method returns the resids of residues that form the interface between
the chain and the set of chains supplied. The difference between the residue's
rASA values in and out of complex must be greater or equal to 10 to be defined
as interface.

This method can be passed either an aref of chains that are complexed with
the chain or a href of ResID => relASA values, where relASA values has been
calculated in the context of a complex of interest.

=cut

sub getInterfaceResidues {
    my $self = shift;

    my $complexResid2RelASAHref;
    
    if (ref $_[0] eq 'ARRAY') {
        my $otherChainsAref = shift;

        $self->read_ASA() if ! $self->has_read_ASA();

        my $atomAref
            = pdb::pdbFunctions::generateAtomAref($self, @{$otherChainsAref});
    
        # Calculate ASAb values of self + otherChains
        my $solv = pdb::solv->new(input => $atomAref);
        my $atomSerialHref = $solv->getOutput();
        $complexResid2RelASAHref = $solv->resid2RelASAHref();
    }
    elsif (ref $_[0] eq 'HASH') {
        $complexResid2RelASAHref = shift;
    }
    else {
        croak "getInterfaceResidues must be passed an aref of chains or an href"
            . " of resID => inComplexRelASA";
    }

    $self->read_ASA if ! %{$self->resid2RelASAHref};
    
    my @interfaceResidues = ();

    foreach my $resid (keys %{$self->resid_index()}) {

        # Skip if this a solvent residue
        next if [values %{$self->resid_index->{$resid}}]->[0]->is_solvent();
        
        my $complexRelASA   = $complexResid2RelASAHref->{$resid}->{allAtoms};
        my $monomericRelASA = $self->resid2RelASAHref->{$resid}->{allAtoms};
        
        confess "Undefined in complex relASA for value for "
            . $self->pdb_code() . $self->chain_id() . " $resid!\n"
                if ! defined $complexRelASA;

        confess "Undefined relASA for value for "
            . $self->pdb_code() . $self->chain_id() . " $resid!\n"
                if ! defined $monomericRelASA;
                
        my $ASAdiff = $monomericRelASA - $complexRelASA;
        
        push(@interfaceResidues, $resid) if $ASAdiff >= 10;
    }
    return @interfaceResidues;
}

=item C<getCDRAtoms>

This method returns an array of chain atoms labelled as CDR
i.e. those atoms where $atom->is_CDR() is TRUE

=cut

sub getCDRAtoms {
    my $self = shift;

    my @CDRAtoms = ();

    foreach my $atom (@{$self->atom_array()}) {
        if ($atom->is_CDR()) {
            push(@CDRAtoms, $atom);
        }
    }
    return @CDRAtoms;
}

=item C<isInContact(CHAINS, CONTACT_MIN, DIST_THRESHOLD)>

CHAINS = ref to array of chains to consider contact with.
CONTACT_MIN = Minimum number of residue-residue contacts that defines a
              chain-chain contact. DEFAULT = 1
DIST_THRESHOLD = atom-atom distance that defines contact. DEFAULT = 4 (ang)

This method checks if the chain is in contact with the other chains passed.
A tolerance can be used to set the minimum of residue-residue contacts that
must exist for the chains to be considered contacting. Default = 1
A threshold can also be set that specifies that maximum atom-atom distance
that will quality residues as touching.
  e.g. $antigenChain->isInContact([$heavyChain, $lightChain], 5, 3)
returns 1 if chain is in contact with any of the other chains passed. 

=cut

sub isInContact {
    my $self = shift;
    my($chainAref, $contactMin, $threshold) = @_;

    # Set default contactMins and threshold if not passed
    $contactMin = 1 if ! $contactMin;
    $threshold = 4 if ! $threshold;
    
    # Create array of atoms from all chains
    my @allAtoms = ();
    
    foreach my $ele ($self, @{$chainAref}) {
        try {
            my $chain = $ele;
            push(@allAtoms, @{$chain->atom_array()});
        }
        catch ($err) {
            # If user has passed a ref to an array containing
            # atoms, rather than chains, push atom onto array
            if (ref $ele eq 'atom') {
                my $atom = $ele;
                push(@allAtoms, $atom);
            }
            else {
                croak $err;
            }
        };
    }
    
    my $cContacts = pdb::chaincontacts->new(input => \@allAtoms,
                                            threshold => $threshold);
    my $contResults = $cContacts->getOutput();
    
    my $contactAref = $contResults->chain2chainContacts([$self], $chainAref);
    
    if (scalar @{$contactAref} >= $contactMin) {
        return 1;
    }
    else {
        return 0;
    }
}

=item C<processAlnStr(%args)>

where %args = (alnStr => ALIGNMENT_STR, includeMissing => BOOL)

ALIGNMENT_STR = Alignment string, e.g. = "AKLWRTPNM----QRCV",
                where "-" = gap in the alignment

This method processes an alignment string in order to label atoms with alnSeq.
If includeMissing is TRUE, residues missing in the structure
will be included in the alignment
example input:
 $chain->processAlnStr(alnStr => $myStr, includeMissing => 1);

=cut

sub processAlnStr {
    my $self = shift;
    my %args = @_;
    my $alnStr = $args{alnStr};

    my $includeMissing
        = exists $args{includeMissing} && $args{includeMissing} ? 1 : 0;

    my %resids = ();
    
    if ($includeMissing) {
        %resids = (%{$self->resid_index()}, %{$self->missing_resid_index()});
    }
    else {
        %resids = %{$self->resid_index()};
    }
    
    # Create array of residue atoms, ordered by resSeq 
    my @residue_atoms
        = map { [values %{$self->resid_index->{$_}}] }
            sort {pdb::pdbFunctions::compare_resids($a, $b)}
                keys %resids;

    my $chainSeqPos = 0;
    
    for (my $i = 0 ; $i < length($alnStr) ; ++$i) {
        if (substr($alnStr, $i, 1) ne '-') {
            # Label atoms of this residue
            # $i + 1, as alnSeq should start from 1
            map {$_->alnSeq($i+1)} @{$residue_atoms[$chainSeqPos]};
            ++$chainSeqPos;
        }
    }
}

=item C<labelAtomsWithClusterSeq>

Label chain atoms with clusterSeq.
clusterSeq = chain cluster id + atom alnSeq.

This can be used to find atoms from residues that are aligned - the cluster id
assures that residues from the same alignment are considered.

=cut

sub labelAtomsWithClusterSeq {
    my $self = shift;
    my $cluster_id = $self->cluster_id();

    # clusterSeq is formed from chain cluster id and atom alnSeq
    foreach my $atom (@{$self->atom_array()}) {
        if (! defined $atom->alnSeq()) {
            # If atom is hetatm, assume that atom is solvent (and therefore
            # should not have an aln or cluster seq)
            if ($atom->is_het_atom()) {
                next;
            }
            else {
                use Data::Dumper;
                print "Atom has no alnSeq: $atom\n";
                print $self->get_sequence(include_missing => 1,
                                          return_type => 1),  "\n";
                print Dumper $self;
                exit;
            }
        }
        else {
            $atom->clusterSeq($cluster_id . "." . $atom->alnSeq());
            
        }
    }
}

# This method allows the creation of chains from the same pdb file as the
# calling chain. The chains will be created by parsing atoms from the chain's
# pdb data
sub createOtherChains {
    my $self = shift;

    # Parse all atoms from pdb data
    my $atomLines = [$self->_parse_ATOM_lines(all => 1)];
    my @atoms = @{$self->_parse_atoms($atomLines)};

    my %chainID2Atoms = ();
    
    # Hash atoms by chain id
    foreach my $atom (@atoms) {
        push(@{$chainID2Atoms{$atom->chainID()}}, $atom);        
    }

    my @chains = ();
    
    foreach my $chainID (keys %chainID2Atoms) {
        next if $chainID eq $self->chain_id();

        my $chain = chain->new(pdb_code => $self->pdb_code(),
                               chain_id => $chainID,
                               atom_array => $chainID2Atoms{$chainID});

        push(@chains, $chain);
    }

    return @chains;
}

__PACKAGE__->meta->make_immutable;

################################ END OF PACKAGE ################################
################################################################################

package patch;

=head1 NAME

patch - A class for creating patch objects, i.e. makepatch output.

=cut

=head1 SYNOPSIS

use pdb::patch;

# Either build directly from makepatch object ...
$patchObject = patch->new($makePatch);

# Or supply a parent pdb or chain object and build from a patch summary line
$patchObject = patch->new(summary => "<patch B.160> B:153 B:154 B:155 B:160",
                          parent_pdb => $pdb4houB);

=cut


=head1 DESCRIPTION

patch extends the pdb class. The main function is to deal with output from
makepatch, a program that creates overlapping surface patches from an input
set of atoms (see makepatch.pm)

=cut

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;
use Carp;
use TCNUtil::GLOBAL qw(&rm_trail);

extends 'pdb';
use overload '""' => 'stringify';

=head1 Methods

=over 12

=cut

### Attributes #################################################################
################################################################################

=item C<central_atom>

Central atom of patch. This is specified during the patch creation process.
If a patch has been built from a summary line, then the central atom of the
patch is guessed to be the C-alpha atom of the patch centre residue.

=cut

has central_atom => (
    is => 'rw',
    isa => 'atom',
    required => 1,
);

=item C<summary>

Summary line of patch. E.g.
     <patch B.160> B:153 B:154 B:155 B:156 B:160
Where B.160 is the resid of the patch centre residue.

A patch can be built from a patch summary, if a parent_pdb object is supplied.

=cut

has summary => (
    is => 'rw',
    isa => 'Str',
    lazy => 1,
    builder => '_build_summary',
);

=item C<parent_pdb>

parent_pdb object. If this has been supplied during the patch creation process
via makepatch, then patch atoms will be shared with the parent pdb. That means
that if you make changes to the patch atoms (e.g. change their co-ordinates)
then the parent pdb will also be affected!

=cut

has 'parent_pdb' => (
    is => 'rw',
    predicate => 'has_parent_pdb',
);



has 'ASAb_hash' => (
    is => 'ro',
    isa => 'HashRef',
    lazy => 1,
    builder => '_build_ASAb_hash',
);

foreach my $ASAtype (qw(ASAc ASAm ASAb)) {
    my $attr = 'total_' . $ASAtype; 
    my $builder = '_build_ASA';
    has $attr => (
        isa => 'Num',
        is => 'ro',
        lazy => 1,
        builder => $builder,
    );
}

has 'is_epitope' => (
    isa => 'Bool',
    is => 'rw',
    predicate => 'has_epitope_label',
    clearer => 'clear_epitope_label',
);

has 'id' => (
    is => 'ro',
    isa => 'Str',
    lazy => 1,
    builder => '_build_patch_id',
);

has 'is_multi_chain' => (
    is => 'ro',
    isa => 'Bool',
    lazy => 1,
    builder => '_build_is_multi_chain',
);

### Attribute Builder Methods ##################################################
################################################################################

around 'get_sequence' => sub {
    my $orig = shift;
    my $self = shift;

    my %arg = @_;

    # Ensure that no missing residues that may be present in patch pdb data
    # are included in sequence of patch
    $arg{include_missing} = 0;

    return $self->$orig(%arg);
};

sub _build_patch_id {
    my $self = shift;

    return $self->pdb_code . $self->summary();
}

sub _build_ASA {
    my $self = shift;

    # Looking at calling attribute to determine ASA type to sum
    my $toBeBuilt = (caller(1))[3];
    # example toBeBuilt = patch::total_ASAb
    my $ASAtype = [split("_", $toBeBuilt)]->[-1];
    
    my $totalASA = 0;
    map {$totalASA += $_->$ASAtype()} @{$self->atom_array()};
    
    return $totalASA;
}

sub _build_ASAb_hash {
    my $self = shift;

    my %ASAmHash = ();
    
    foreach my $residAtomHref (values %{$self->resid_index()}) {
        my $totalASA = 0;
        my $residAtomAref = [values %{$residAtomHref}];
        foreach my $atom (@{$residAtomAref}) {
            croak "ASAb not defined for atom:\n$atom"
                if ! defined $atom->ASAb();
        }
        map { $totalASA += $_->ASAb() } @{$residAtomAref};
        my $key
            = $residAtomAref->[0]->alnSeq() . $residAtomAref->[0]->resName();

        $ASAmHash{$key} = $totalASA;
    }
    return \%ASAmHash;
}

sub _build_summary {
    my $self = shift;

    my @chainIDAndResSeqs = ();

    foreach my $chain (keys %{$self->atom_index}) {
        foreach my $resSeq (keys %{$self->atom_index->{$chain}}) {
            push(@chainIDAndResSeqs, [$chain, $resSeq]);
        }
    }

    @chainIDAndResSeqs
        = sort { $a->[0] cmp $b->[0] ||
                     pdb::pdbFunctions::compare_resSeqs($a->[1],$b->[1]) }
            @chainIDAndResSeqs;

    my $cA = $self->central_atom();

    my $resSeqListStr
        = join(" ", map {$_->[0] . ":" . $_->[1]} @chainIDAndResSeqs);
    
    my $summary = "<patch " . $cA->chainID() . "." . $cA->resSeq()
        . "> $resSeqListStr\n";

    return $summary;
}

sub _build_is_multi_chain {
    my $self = shift;

    my %chainIDs = map {$_->chainID() => 1} @{$self->atom_array()};

    return keys %chainIDs > 1 ? 1 : 0;
}

### BUILD Method ###############################################################
################################################################################

sub BUILD {
    my $self = shift;

    my @atom_array = ();

    if ($self->has_parent_pdb()) {
        
        # If parent pdb is actually a ref to an array of chains, then create a
        # multi-chain atom hash. Otherwise, take atom hash of single parent pdb
        my $atomSerialHref
            = ref $self->parent_pdb eq 'ARRAY' ?
                pdb::multiChain::multiChainAtomSerialHref($self->parent_pdb())
              : $self->parent_pdb->atom_serial_hash();

        
        # Replace atoms with corresponding atoms from parent pdb
        # (including central atom)
        foreach my $atom (@{$self->atom_array()}) {
            my $parentAtom = $atomSerialHref->{$atom->serial()};
            if ($atom->serial eq $self->central_atom->serial) {
                $self->central_atom($parentAtom);
            }
            push(@atom_array, $parentAtom);
        }
        $self->atom_array(\@atom_array);
    }
}

### Around Modifiers ###########################################################
################################################################################

# Allow a makepatch object to be passed directly to new method
# Or, if summary line and parent pdb object has been passed,
# parse summary line
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;

    my %arg = int (scalar @_ / 2) == (scalar @_ / 2) ? @_ : ();
    
    if (ref $_[0] eq 'makepatch') {
        my $makepatch = $_[0];
        
        if ( ref $makepatch->output()->[0] eq 'local::error' ) {
            croak $makepatch->output()->[0];
        }

        my %arg
            = ( central_atom => $makepatch->central_atom,
                pdb_data     => $makepatch->output,
                pdb_code     => $makepatch->pdb_code,
            );

        
        if ((! $makepatch->new_atoms) && $makepatch->has_pdb_object) {
            # Patch will contain atoms from parent pdb, rather than new atoms
            $arg{parent_pdb} = $makepatch->pdb_object;
        }
        $class->$orig(%arg);
    }
    elsif ($arg{summary} && $arg{parent_pdb}) {
        
        # Build from summary and parent pdb
        my @resids = parseSummaryLine($arg{summary});
        
        my @atoms = map {values %{$arg{parent_pdb}->resid_index->{$_}}} @resids;
        
        $arg{atom_array} = \@atoms;
        
        if ($arg{central_atom}) {
            # To be consistent, ensure central atom is from parent pdb
            $arg{central_atom}
                = $arg{parent_pdb}->atom_index->{$arg{central_atom}->serial()};
        }
        else {
            # If no central atom has been supplied, assign central residue
            # C-alpha
            
            foreach my $atom (@atoms) {
                next if $atom->name() ne 'CA';

                # First C-alpha should be from central residue
                croak "Central residue has no C-alpha atom!"
                    if $atom->resid() ne $resids[0]; # 1st ele is central resid

                $arg{central_atom} = $atom;
                last;
            }
        }
        $class->$orig(%arg);
    }
    else {
        $class->$orig(@_);
    }
};

around '_parse_ATOM_lines' => sub {

    my $orig = shift;
    my $self = shift;
      
    my @patch_ATOM_lines = ();
    
    foreach my $line ( $self->$orig() ) {
        if ( $line =~ /^(?:ATOM|HETATM)/ && substr($line, 60, 6) =~ /1\.00/ ){
            push(@patch_ATOM_lines, $line);
        }
    }
    
    croak "No ATOM lines were parsed for patch"
        if ! @patch_ATOM_lines;

    return @patch_ATOM_lines;
};

### Methods ####################################################################
################################################################################

sub parseSummaryLine {
    my $summaryLine = shift;

    # example summary line : <patch G.409> G:335 G:397 G:398 G:407 G:408 G:409
    # parse all resids
    my @resids = $summaryLine =~ /(\w+[\.:]-*\w+)/g;
    
    # change any : separators to .
    map {s/:/./} @resids;
    
    my $centralResid = shift @resids;

    # Remove repeat of central residue
    @resids = grep {$_ ne $centralResid} @resids;
    
    return ($centralResid, @resids);
}

sub stringify {
    my $self = shift;
    return $self->pdb_code . ":" . $self->central_atom->chainID()
        . "." . $self->central_atom->resSeq();
}

__PACKAGE__->meta->make_immutable;

################################################################################
################################################################################

package atom;

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;
use Math::Trig;
use TCNUtil::GLOBAL qw(&rm_trail);
use TryCatch;

use Carp;

### Attributes #################################################################
################################################################################

has 'ATOM_line' => (
    isa => 'Str',
    is  => 'rw',
);

has [qw(name resName element charge resSeq kabatSeq chothiaSeq ichothiaSeq
        alnSeq clusterSeq) ]
    => (is => 'rw', isa => 'Str');

foreach my $name ('altLoc', 'chainID', 'iCode') {
    my $predicate = 'has_'   . $name;
    my $clearer   = 'clear_' . $name;
    
    has $name => (is => 'rw',
                  isa => 'Str',
                  predicate => $predicate,
                  clearer => $clearer); 
}

has ['serial'] => (is => 'rw', => isa => 'Int');

foreach my $name (qw(radius ASAm ASAc ASAb x y z occupancy tempFactor)) {
    my $predicate = 'has_' . $name;
    
    has $name => (is => 'rw', isa => 'Num', predicate => $predicate);
}

has 'resid' => (
    is => 'ro',
    isa => 'Str',
    lazy => 1,
    builder => '_get_resid',
);

my @labels = qw(is_het_atom is_terminal is_solvent
                is_CDR is_epitope is_interface is_modified
                is_from_modRes);

foreach my $label (@labels) {
    has $label => (
        isa => 'Bool',
        is => 'rw',
        default => 0,
    );
}

foreach (qw(HbDonor HbAcceptor SSbond)) {
    has $_ => (
        isa => 'atom',
        is  => 'rw',
        predicate => "has_" . $_,
    );
}

### Attribute Builder Methods ##################################################
################################################################################

sub _get_resid {
    my $self = shift;
    my $resid = join(".", ($self->chainID, $self->resSeq));

    if ($self->has_iCode()) {
        $resid .= $self->iCode();
    }
    
    return $resid;
}

### BUILD Method ###############################################################
################################################################################

sub BUILD {
    
    my $self = shift;

    my $ATOM_line = $self->ATOM_line();

    return if ! $ATOM_line;

    # Avoid substr complaining if there are missing columns at end of string
    $ATOM_line = pack ( "A81", $ATOM_line );
    my %record
        = ( ATOM => rm_trail(substr($ATOM_line, 0, 6)),
            serial =>  rm_trail(substr($ATOM_line, 6, 5)),
            name => rm_trail(substr($ATOM_line, 12, 4)),
            altLoc => rm_trail(substr($ATOM_line, 16, 1)),
            resName => rm_trail(substr($ATOM_line, 17, 3)),
            chainID => rm_trail(substr($ATOM_line, 21, 1)),
            resSeq => rm_trail(substr($ATOM_line, 22, 4)),
            iCode => rm_trail(substr($ATOM_line, 26, 1)),
            x => rm_trail(substr($ATOM_line, 30, 8)),
            y => rm_trail(substr($ATOM_line, 38, 8)),
            z => rm_trail(substr($ATOM_line, 46, 8)),
            occupancy => rm_trail(substr($ATOM_line, 54, 6)),
            tempFactor => rm_trail(substr($ATOM_line, 60, 6)),
            element => rm_trail(substr($ATOM_line, 76, 2 )),
            charge => rm_trail(substr($ATOM_line, 78, 2)),
        );

    croak "It looks like you're trying to parse a non-ATOM or HETATM line: "
        . "$ATOM_line"
            if ! ($record{ATOM} eq 'ATOM' || $record{ATOM} eq 'HETATM') ;

    if ( $record{ATOM} eq 'HETATM' ) {
        $self->is_het_atom(1);
    }
    
    delete $record{ATOM};
    
    foreach my $value (keys %record) {
        next if $record{$value} eq '';
        $self->$value( $record{$value} );
    }
}

### Methods ####################################################################
################################################################################

use overload '""' => \&stringify, fallback => 1;

sub stringify {
    my $self = shift;

    my @arg = ();
    
    foreach my $arg (@_) {
        if ( defined $arg && $arg =~ /\S/ ) {
            push(@arg, $arg);
        }
    }
    
    my %hash = ();
    
    if (@arg) {
        croak "stringify must be passed a hashref: this hash must contain "
            . "column replacements. e.g. { tempFactor => ASAm } to replace "
            . "tempFactor with ASAm"
                if ref $_[0] ne 'HASH';

        %hash = %{ $_[0] };
    }
    
    my @ordered_attr = qw(serial name altLoc resName chainID
                          resSeq iCode x y z occupancy
                          tempFactor element charge);
    
    # Make column replacements
    for my $i ( 0 .. @ordered_attr - 1 ) {
        if ( exists $hash{ $ordered_attr[$i] } ) {
            my $predicate = 'has_' . $ordered_attr[$i];
            
            $ordered_attr[$i]
                = $self->$predicate ? $hash{ $ordered_attr[$i] } : ' ' ;
        }
    }
    
    @ordered_attr = map { $self->$_ } @ordered_attr;
    
    # Process 'name' to match pdb formatting by padding with whitespace
    $ordered_attr[1]
        =  length $ordered_attr[1] == 0 ? ' ' x 4
         : length $ordered_attr[1] == 1 ? ' ' . $ordered_attr[1] . ' ' x 2
         : length $ordered_attr[1] == 2 ? ' ' . $ordered_attr[1] . ' '
         : length $ordered_attr[1] == 3 ? ' ' . $ordered_attr[1] 
         : $ordered_attr[1];

    unshift (@ordered_attr, ( $self->is_het_atom ? 'HETATM' : 'ATOM' ) );

    for my $i ( 0 .. @ordered_attr) {
        if ( ! defined $ordered_attr[$i] ) {
            $ordered_attr[$i] = '';
        }
    }

    my $string
        = sprintf(  '%-6.6s%5.5s' . ' ' . '%s%1.1s%3.3s %1.1s%4.4s%1.1s   '
                   .'%8.3f' x 3 .  '%6.2f' x 2 . ' ' x 10 . '%2.2s'
                   . '%2.2s' ,
                    @ordered_attr );

    $string .= "\n";
    
    return $string;
}

sub stringify_ter {
    my $self = shift;

    # Temporarily increment serial
    $self->serial( $self->serial() + 1);
    
    my $string = $self->stringify();

    # Deincrement serial
    $self->serial( $self->serial() - 1 );
    
    # Clip string to TER line length
    $string = pack( "A27", $string );

    # Modify start of string
    substr($string, 0, 6) = 'TER   ';

    # Add newlines
    $string = $string . "\n\n";
    
    return $string;
}

# This method tests to see if any of the passed atoms are in contact with self
# atom. Contact is determined using the distance threshold supplied
# (Default = 4A)
# Input: ref to array of atoms, distance threshold (optional)
# e.g. $atom->anyContacting($atomsAref, 4)
sub anyContacting {
    my ($self, $atomsAref, $distance) = @_;

    $distance = 4 if ! $distance;
    
    foreach my $atom (@{$atomsAref}) {
        if ($self->contacting($atom)) {
            return 1;
        }
    }
    return 0;
}

# This method tests to see if two atoms are contacting eachother,
# based on a contact distance threshold (default = 4A)
sub contacting {
    my ($self, $atom2test, $distance) = @_;

    $distance = 4 if ! $distance;
    my $distSq = $distance ** 2;
    
    my $ax = $self->x();
    my $ay = $self->y();
    my $az = $self->z();

    my $bx = $atom2test->x();
    my $by = $atom2test->y();
    my $bz = $atom2test->z();

    my $d = (($ax - $bx)**2) + (($ay - $by)**2) + (($az - $bz)**2);

    if ($d <= $distSq) {
        return 1;
    }
    else {
        return 0;
    }
}

sub readRadius {
    my $self         = shift;
    my $radiusFinder = shift;
    my $radius;
    
    try {
        $radius = $radiusFinder->findRadiusOfAtomFromNames($self->resName, $self->name);
    }
    catch {
        try {
            $radius = $radiusFinder->findRadiusOfAtomFromElement($self->element);
        }
        catch {
            $radius = 1.80;
        }
    };
    
    $self->radius($radius);
}

__PACKAGE__->meta->make_immutable;


################################################################################
################################################################################

1;


__END__



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
 
