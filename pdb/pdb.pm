package pdb;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use local::error;

use GLOBAL qw(&rm_trail &three2one_lc &is_int);

use Carp;

use TryCatch;
use Storable;

use Math::VectorReal qw(:all);
use Math::MatrixReal;
use Math::Trig;

use pdb::xmas2pdb;
use pdb::getresol;
use pdb::rotate2pc;
use pdb::get_files;
use pdb::asurf64;

# Subtypes

### Attributes

has 'pdb_code' => (
    isa => 'Str',
    is  => 'rw',
    builder => '_build_pdb_code',
    lazy => 1,
);

sub _build_pdb_code {
    my $self = shift;

    if ($self->has_parent_pdb()) {
        return $self->parent_pdb->pdb_code();
    }
    else {
        croak "No pdb code has been set!";
    }
}


for my $name ( 'pdb', 'xmas' ) {
    
    my $att_file = $name . '_file';
    my $att_data   = $name . '_data'  ;
    
    has $att_file => (
        isa => 'FileReadable',
        coerce => 1,
        # See types.pm for FileReadable coercion from ArrayRef[Str]
        is => 'rw',
        predicate => 'has_' . $att_file,
        builder => "_get_$att_file",
        lazy => 1,
    );

    has $att_data => (
        isa => 'ArrayRef',
        is => 'rw',
        lazy => 1,
        builder => '_build_' . $att_data,
        predicate => 'has_' . $att_data,
    );
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

sub _get_xmas_file {
    my $self = shift;

    my $pdbCode = $self->pdb_code();
    my $fName = eval {pdb::get_files->new(pdb_code => $pdbCode)->xmas_file()};

    if (! $fName) {
        croak "No xmas file specified: "
            . "Attempted to get xmas file from pdb code '$pdbCode', $@";
    }
    
    return $fName;
}

sub _build_pdb_data {
    my $self = shift;
    return $self->_build_data('pdb');
}
sub _build_xmas_data {
    my $self = shift;
    return $self->_build_data('xmas');
}


sub _build_data {
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

has 'remark_hash' => (
    isa => 'HashRef',
    is => 'ro',
    builder => '_build_remark_hash',
    lazy => 1,
);

has 'atom_array' => (
    isa => 'ArrayRef[atom]',
    is  => 'rw',
    lazy => 1,
    builder => '_parse_atoms',
);

# Index is form of chain -> resSeq -> atom_name -> atom
has 'atom_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_atom_index',
);

has 'terminal_atom_index' => (
    isa => 'ArrayRef[atom]',
    is  => 'rw',
    default => sub { [] },
);

# Index in form of resid -> atom_name -> atom
# Resid = chainResSeq
has 'resid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_resid_index',
);

# Index in form of resid -> {}
# Where each resid has been identified as missing from structure
# (see missing_residues)
has 'missing_resid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_missing_resid_index',    
);

# This index is the same as resid_index, except that resids are prepended with
# the object's pdb code
has 'pdbresid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_pdbresid_index',
);

has 'atom_serial_hash' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_atom_serial_hash',
);

has 'multi_resName_resid' => (
    isa => 'HashRef',
    is => 'rw',
    default => sub { {} },
);

# Selects one altLoc atom for a given residue atom according to highest
# occupancy
has 'altLoc_cleanup' => (
    isa => 'Bool',
    is  => 'rw',
    default => 1,
);

# Avoid parsing hydrogen atoms from pdb data
has 'hydrogen_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

has 'het_atom_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

has 'has_read_ASA' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

has 'solvent_cleanup' => (
    isa => 'Bool',
    is => 'rw',
    default => 0,
);

has 'experimental_method' => (
    is => 'rw',
    isa => 'Str',
    predicate => 'has_experimental_method',
    lazy => 1,
    builder => '_build_experimental_method',
);

has 'resolution' => (
    is => 'rw',
    isa => 'Num',
    predicate => 'has_resolution',
    lazy => 1,
    builder => '_build_resolution',
);

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

has 'resid2RelASAHref' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
    lazy => 1,
);


# Consume antigen role
with 'pdb::antigen';

### Methods

# Attempt to build experimental_method, resolution and r_factor from getresol
# object. Also check if pdb is multi-model
sub BUILD {
    my $self = shift;

    return if ! $self->has_pdb_file();

    my ($expMethod, $resolution, $rValue) = eval {$self->_run_getresol()};

    $self->experimental_method($expMethod) if $expMethod;
    $self->resolution($resolution) if $resolution;
    $self->r_value($rValue) if $rValue;
     
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

sub _build_experimental_method {
    my $self = shift;

    my($expMethod, $resolution, $rValue) = $self->_run_getresol();

    # Set the other two values
    $self->resolution($resolution);
    $self->r_value($rValue);
    
    return $expMethod;
}

sub _build_resolution {
    my $self = shift;

    my($expMethod, $resolution, $rValue) = $self->_run_getresol();

    # Set the other two values
    $self->experimental_method($expMethod);
    $self->r_value($rValue);
    
    return $resolution;
}

sub _build_r_value {
    my $self = shift;

    my($expMethod, $resolution, $rValue) = $self->_run_getresol();

    # Set the other two values
    $self->experimental_method($expMethod);
    $self->resolution($resolution);
   
    return $rValue;
}

sub _run_getresol {

    my $self = shift;
    
    my $getresol = pdb::getresol->new(pdb_file => $self->pdb_file);

    my $expMethod;
    my $resolution;
    my $rValue;

    my $ret = $getresol->run();
    
    if (ref $ret  ne 'local::error'){
        $expMethod  = $getresol->experimental_method();
        $resolution = $getresol->resolution();
        $rValue     = $getresol->r_value();
    }
    else {
        croak $ret;
    }

    return ($expMethod, $resolution, $rValue);
}


sub _build_atom_serial_hash {
    my $self = shift;

    my %atomSerialHash = ();

    foreach my $atom (@{$self->atom_array()}) {
        $atomSerialHash{$atom->serial} = $atom;
    }

    return \%atomSerialHash;
}

# Creates hash of form:
# remarkNumber => Ref to array of scalar refs to pdb_data lines
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
        my $resSeq  = substr($line, 21, 5);

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

sub _parse_atoms {
    my $self = shift;
    
    my @ATOM_lines = $self->_parse_ATOM_lines();
    
    my @atoms = ();
    
    my $aL_clean = $self->altLoc_cleanup();
    my %altLoc = ();
    
    my $h_clean = $self->hydrogen_cleanup();
    my $HETATM_clean = $self->het_atom_cleanup();
    my $solvent_clean = $self->solvent_cleanup();
    
    my %ter = ();
    
    # Used to find multi-residue resids
    my %test_unique = ();
    
    foreach my $line (@ATOM_lines) {
        
        if ($line =~ /^TER/) {
            # Label previous atom as terminal
            my $prevAtom = $atoms[-1];
            $prevAtom->is_terminal(1);
        
            $ter{$prevAtom->chainID} = []
                if ! exists $ter{$prevAtom->chainID};
            
            push(@{$ter{$prevAtom->chainID}}, $prevAtom->serial());
            next;
        }
        
        my $atom = atom->new(ATOM_line => $line);

        my $chain   = $atom->chainID();
        my $resSeq  = $atom->resSeq();
        my $resName = $atom->resName();

        if (! exists $test_unique{$chain}) {
            $test_unique{$chain} = {};
        }
        
        if (! exists $test_unique{$chain}->{$resSeq}) {
            $test_unique{$chain}->{$resSeq} = $resName;
        }
        else {
            my $recorded_resName = $test_unique{$chain}->{$resSeq};
            if ($recorded_resName ne $resName){
               
                my %multi = %{$self->multi_resName_resid};
                $multi{ $atom->_get_resid } = 1;
                
                $self->multi_resName_resid({%multi});
            }
        }
        
        next if ($h_clean && $atom->element eq 'H'
                      || $HETATM_clean && $atom->is_het_atom);
        
        if ($aL_clean && $atom->has_altLoc) {

            my $string
                =  $atom->name . $atom->resName . $atom->chainID
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

    if ($aL_clean) {
        foreach my $arr (values %altLoc) {

            # Sort altlocs by occupancy
            my @sorted
                = map { $_->[0] }
                    sort { $b->[1] <=> $a->[1] }
                        map { [ $_, $_->occupancy ] }
                            @{$arr};

            my $top = $sorted[0];
            
            # Assign top occupant's x,y,z and serial to first element of array
            # (as this is the member that is contained within atom list)
            foreach my $attr (qw(x y z serial)) {
                $arr->[0]->$attr($top->$attr());
            }
            
            # Clear altLoc
            $arr->[0]->clear_altLoc();
        }        
    }


    # Order all atoms by chain, then serial
    my %chain = ();
    foreach my $atom (@atoms) {
        my $chainID = $atom->chainID();

        if (! exists $chain{$chainID}) {
            $chain{$chainID} = {$atom->serial => $atom};
        }
        else {
            $chain{$chainID}->{$atom->serial} = $atom;
        }
    }    
    
    my @sorted_atoms = ();
    my @term_array   = ();
    
    foreach my $chainID (sort keys %chain) {
        my @chain_atoms = ();
        
        foreach my $atom ( sort {$a <=> $b} keys %{$chain{$chainID}}) {
            push(@chain_atoms, $chain{$chainID}->{$atom});
        }

        my @chain_terminal_index = ();
        
        # Label terminal atoms
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
 
        # Determine solvent atoms
        if( @chain_terminal_index == 1 ) {
            my $start = $chain_terminal_index[0];
            # Label all atoms after chain terminal as solvent
            for my $i ( $start + 1 .. @chain_atoms - 1 ) {
                $chain_atoms[$i]->is_solvent(1);

            }
        }
        elsif ( @chain_terminal_index > 1 ) {
            # Determine which terminal signals end of chain and which signals
            # end of solvent
            
            my $index
                = _determine_solvent(\@chain_atoms,
                                    [ sort { $a <=> $b }
                                          @chain_terminal_index ] );

            
            for ( my $i = 0 ; $i >= @chain_atoms ; $i++ ) {
                # Label those atoms not in range of returned index and
                # returned index - 1, as solvent
                unless ($i > $chain_terminal_index[$i - 1]
                            && $i < $chain_terminal_index[$i]) {
                    $chain_atoms[$i]->is_solvent(1);
                }
            }
        }

        if ($solvent_clean) {
            my @nonsolvent = ();
            foreach my $atom (@chain_atoms) {
                push(@nonsolvent, $atom) if ! $atom->is_solvent();
            }
            @chain_atoms = @nonsolvent;
        }
        
        push(@sorted_atoms, @chain_atoms);
    }

    $self->terminal_atom_index( [ @term_array ] );

    croak "No atom objects were created" if ! @sorted_atoms;

    return [ @sorted_atoms ];
}

# Determines the terminal atom that signals the end of the non-solvent chain
# segment
sub _determine_solvent {
    my( $chain_atoms, $chain_terminal_index ) = @_;

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
    
    my @array = @{ $self->pdb_data };  
 
    my @ATOM_lines = ();

    foreach my $line (@array) {
        if ($line =~ /^(?:ATOM|HETATM|TER) /) {
            push(@ATOM_lines, $line);
        }
    }
    
    croak "No ATOM, HETATM or TER lines parsed from pdb data"
        if ! @ATOM_lines;

    return @ATOM_lines;
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


sub get_sequence {
    my $self = shift;

    croak "pdb: " . $self->pdb_code()
        . " can't get_sequence for pdb containing multi-resName residues: "
        . Dumper $self->multi_resName_resid()
            if %{ $self->multi_resName_resid };
    
    my $USAGE
        = 'get_sequence(chain_id => (chain_id), return_type => ( 1 | 3 ), '
        . 'include_missing_res => (1|0, default = 0)';
    
    my %arg = @_;

    # Check args are okay
    if (   ! (   exists  $arg{return_type} && exists  $arg{chain_id} )
        || ! (   defined $arg{return_type} && defined $arg{chain_id} )
        || ! ( $arg{return_type} eq '1' || $arg{return_type} eq '3'  ) ) {
        croak $USAGE;
    }

    my $inc_missing = 0;
    
    if (exists $arg{include_missing_res}) {
        $inc_missing = $arg{include_missing_res};
    }
    
    my %chain_resSeq_index = ();

    # Filter out residues that are not from specified chain
    foreach my $chain (keys %{$self->atom_index()}) {
        next if $chain ne $arg{chain_id};
        foreach my $resSeq (keys %{$self->atom_index->{$chain}}){
            $chain_resSeq_index{$resSeq}
                = $self->atom_index->{$chain}->{$resSeq};
        }
        
    }

    if ($inc_missing) {
        # Get missing residues for specified chain and add them to chain hash
        if (exists $self->missing_residues->{$arg{chain_id}}) {

            my %missingResSeq = %{$self->missing_residues->{$arg{chain_id}}};
            
            foreach my $resSeq (keys %missingResSeq) {
                $chain_resSeq_index{$resSeq} = $missingResSeq{$resSeq};
            }
        }
    }
    
    #  Sort resSeqs
    my @sortedResSeqs = sort {pdb::pdbFunctions::compare_resSeqs($a, $b)}
        keys %chain_resSeq_index;


    my @residues = ();

    foreach my $resSeq (@sortedResSeqs){
        # If hashed value is a ref to a hash where values are atoms
        if (ref $chain_resSeq_index{$resSeq} eq 'HASH'){
            # If not solvent, get resName of first atom
            my $atom = [values %{$chain_resSeq_index{$resSeq}}]->[0];
            if (! $atom->is_solvent()) {
                my $resName = $atom->resName();
                push(@residues, $resName);
            }
        }
        else {
            # If hashed value is a string containing 3lc resname
            push(@residues, $chain_resSeq_index{$resSeq});
        }
    }
    
    if ($arg{return_type} == 1) {
        for my $i (0 .. @residues - 1) {
            my $onelc = eval { three2one_lc($residues[$i]) };
            if ($@) {
                $onelc = 'X';
            }
            $residues[$i] = $onelc;
        }
    }
    return @residues;
}

# This method returns a string containing a FASTA-formatted sequence
# for the given chain. The option is given to include residues missing from
# the structure (default is false)
# Example input:
# $chain->getFASTAStr(chain_id => "A", header => "myHeader", includeMissing => 1)
sub getFASTAStr {
    my $self = shift;
    my %args = @_;
    
    my $chainID = $args{chain_id} or croak "No chain id was passed!";
    my $includeMissing
        = exists $args{includeMissing} ? $args{includeMissing} : 0;
    
    my $header
        = exists $args{header} ? $args{header}
        : ">" . $self->pdb_code() . $chainID . "\n";
   
    
    my @seq = $self->get_sequence(chain_id => $chainID,
                                  return_type => 1,
                                  include_missing_res => $includeMissing);
    
    my $seqStr = join("", @seq) . "\n";
    $header = ">" . $self->pdb_code() . $chainID . "\n" if ! $header;
   
    return $header . $seqStr;
}

# This method returns an array of atoms for residues within
# the given range. Range acts like perl indices
# i.e 0 = first residue, -1 = last residue
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

# This method returns an array ref of atom arrays, sorted by
# chain then resSeq
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

# Given an array of atoms, sorts the atoms by serial and returns the array
sub sort_atoms {
    my @atoms = @_;
    
    @atoms = sort {$a->serial() <=> $b->serial()} @atoms;

    return @atoms;
}

sub read_ASA {

    my $self = shift;
    my $xmas2pdb;

    my @errors = ();

   # asurf64 must also be run in order to get relative ASAs
    my $asurf = pdb::asurf64->new(input => $self);
    my $atomSerial2ASARadHref = $asurf->getOutput();
    
    # Use xmas2pdb object if it has been passed, otherwise create one
    if (@_) {
        $xmas2pdb = shift;
        
        croak "read_ASA must be passed an xmas2pdb object"
            if ref $xmas2pdb ne 'xmas2pdb';
        
        @errors = $self->_parseXmas2PDBOutput($xmas2pdb);
    }
    else {

        my $ASAType = ref $self eq 'pdb' ? 'ASAc' : 'ASAm' ;
        
        # Use asurf64 per atom output
        foreach my $atom (@{$self->atom_array()}) {
            next if $atom->is_solvent()
                || ! $atomSerial2ASARadHref->{$atom->serial()};
            
            my ($ASA, $radius) = @{$atomSerial2ASARadHref->{$atom->serial()}};
            
            $atom->$ASAType($ASA);
            $atom->radius($radius);
        }
    }
    
    $self->resid2RelASAHref($asurf->resid2RelASAHref());    
    $self->has_read_ASA(1);
    
    return \@errors;
}

sub _parseXmas2PDBOutput {
    my $self = shift;
    my $xmas2pdb = shift;

    
    my $form = $xmas2pdb->form();

    my $attribute = 'ASA' . ( $form eq 'monomer' ? 'm' : 'c' ) ;
    
    my  %ASAs  = ();
    my  %radii = ();
    
    foreach my $line ( @{ $xmas2pdb->output() } ){
        next unless $line =~ /^(?:ATOM|HETATM) /;

        my $serial = rm_trail( substr($line, 6, 5) );
        my $radius = rm_trail( substr($line, 54, 6) );
        my $ASA    = rm_trail( substr($line, 60, 6) );

        $ASAs{$serial} = $ASA;
        $radii{$serial}= $radius;
    }

    croak "Nothing parsed from xmas2pdb output!" if ! ( %ASAs && %radii );

    my @errors = ();
    
    foreach my $atom ( @{ $self->atom_array() } ) {
        my $serial = $atom->serial();

        if ( ! exists $ASAs{$serial} ) {
            my $message =  "Nothing from xmas2pdb object for atom "
                          . $atom->serial();
            
            my $error
                = local::error->new( message => $message,
                                     type => 'ASA_read',
                                     data => { xmas2pdb => $xmas2pdb,
                                               atom     => $atom, },
                                 );
            push(@errors, $error);
            
            next;
        }   
        $atom->radius( $radii{$serial} );
        $atom->$attribute( $ASAs{$serial} );
    }
    return @errors;
}


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
        : ref $self eq 'chain' ? 'ASAm'
        : 'ASAc';
    
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

    # Remove any . separators from resid, e.g. A.133 -> A133
    #$resid =~ s/\.//;
    
    croak "resid '$resid' was not found in resid index"
        if ! exists $self->resid_index->{$resid};

    my @atoms = values %{ $self->resid_index->{$resid} };

    # Avoid hydrogen atoms
    @atoms = grep {$_->element() ne 'H'} @atoms;
    
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

# Returns hash chainSeq => resSeq, by reversing keys and values of
# map_resSeq2chainSeq
sub map_chainSeq2resSeq {
    my $self = shift;
    
    my $chain_id
        = shift or croak "map_resSeq2chainSeq must be passed a chain id";
    
    croak "pdb: " . $self->pdb_code() . " no residues found for chain "
        . " $chain_id" if ! exists $self->atom_index->{$chain_id};
    
    my %resSeq2chainSeq = $self->map_resSeq2chainSeq($chain_id);
    
    my %return_map
        = map { $resSeq2chainSeq{$_} => $_ } keys %resSeq2chainSeq;
    
    return %return_map;
}
    
# Maps resSeq numbers to chainSeq count numbers (equivalent to pdbcount num
# in pdbsws. Returns hash resSeq => chainSeq
sub map_resSeq2chainSeq {
    my $self = shift;

    my $chain_id
        = shift or croak "map_resSeq2chainSeq must be passed a chain id";

    croak "pdb: " . $self->pdb_code() . " no residues found for chain "
        . " $chain_id" if ! exists $self->atom_index->{$chain_id};
    
    my $chainSeq = 0;
    my %return_map = ();
    
    my $chain_id_h = $self->atom_index->{$chain_id};
    
    foreach my $resSeq ( sort { $a <=> $b } keys %{ $chain_id_h  } ) {

        my $atom = [ values %{ $chain_id_h->{$resSeq} } ]->[0];
        
        next if $atom->is_solvent();
        
        ++$chainSeq;

        $return_map{$resSeq} = $chainSeq;
    }
    
    return %return_map;
}

# Method to create chain objects from pdb object. Returns array of chains.
# An array of chain ids can be passed to the method; if this is the case then
# the returned chains will be in the order specified in the passed array
sub create_chains {
    my $self = shift;
    my @passed_chain_ids = @_;

    # If chain ids have been passed, check that all chain ids are found in
    # this pdb
    if (@passed_chain_ids) {
        my %chain_ids = map { $_ => 1 } @{$self->get_chain_ids()};
        
        foreach my $chain_id (@passed_chain_ids) {
            croak "Passed chain_id $chain_id was not found in pdb!"
                if ! exists $chain_ids{$chain_id};
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
                               atom_array => $atoms{$chain_id});

        $chains{$chain->chain_id()} = $chain;
    }

    my @return_chains = ();
    
    # If chain_ids have been specified, return chains in specified order
    if (@passed_chain_ids) {
        foreach my $chain_id (@passed_chain_ids) {
            push(@return_chains, $chains{$chain_id});
        }
    }
    else {
        @return_chains = values %chains;
    }
    return @return_chains;
}

# Returns arrayref containing chainIDs found in pdb
sub get_chain_ids {
    my $self = shift;

    my @chain_ids = keys (%{$self->atom_index()});

    return \@chain_ids;
}

# Determines antibody pairs found within pdb.
# INPUT:
# This fuction can optionally be passed an array of pdb chains.
# e.g.
#  $pdb->getAbPairs()
#  $pdb->getAbPairs(@chains)
#  $pdb->getAbPairs($pdb->create_chains())
# OUTPUT:
# Returns three refs to arrays of:
#  1. Antibody Pairs:  [ [$heavyChain, $lightChain, $numContacts], ... ]
#  2. Unpaired Chains: [ $chainA, $chainB, ... ]
#  3. scFv Chains:     [ $chainScFv, ... ]
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
        my $chainType = $chain->isAbVariable();
        if (! $chainType) {
            $chainType = 'antigen';
        }
        elsif ($chainType eq 'Heavy' || $chainType eq 'Light') {
            $pairedHref->{$chain->chain_id()} = 0;
        }
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

# This method rotates atoms of the pdb so that the PC1 and PC2 of the atoms
# lay on the x and y axes respectively. An array of atoms or resids can be
# passed; in this case, the PC1 and PC2 of this subset of atoms (or atoms from
# resids) can be used. 
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
    my $rotationMatrix = rotate2pc::rotate2pc(@vectors);

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

# This method checks the number of +z and -z atom co-ordinates and flips the pdb
# 180degrees around the z axis so that the z axis points through the "body" of
# the pdb; i.e. flips the pdb so that you are NOT looking through the body of
# the pdb, if you are looking down the z axis
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

# This method sets all is_epitope labels to 0
sub clearEpitopeLabels {
    my $self = shift;

    map { $_->is_epitope(0) } @{$self->atom_array()};
}

# This method saves the pdb in a binary file. A filename can be passed,
# otherwise the pdb code is used. Returns file name of binary file
# e.g
# $pdb->store();
# my $storedPDBFname = $pdb->store();
# $pdb->store("myStoredPDB.obj");
sub storeInFile {
    my $self = shift;
    my ($fName) = @_;

    $fName = $self->pdb_code() . ".obj" if ! $fName;

    store $self, $fName;

    return $fName;
}

__PACKAGE__->meta->make_immutable;

package chain;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use Carp;
use TryCatch;

use pdb::pdbsws;
use pdb::idabchain;
use pdb::kabatnum;
use pdb::chaincontacts;

use write2tmp;

extends 'pdb';

# Attributes

has 'chain_id' => (
    isa => 'ValidChar',
    is => 'rw',
    required => 1,
);

has 'accession_codes' => (
    isa => 'ArrayRef[ValidAC]',
    is => 'ro',
    lazy => 1,
    builder => '_get_acs',
);

has 'chain_length' => (
    isa => 'Int',
    is => 'ro',
    lazy => 1,
    builder => '_build_chain_length',
);

has 'is_ab_variable' => (
    isa => 'Str',
    is => 'rw',
    lazy => 1,
    builder => 'isAbVariable',
);

has 'is_het_chain' => (
    isa => 'Bool',
    is => 'ro',
    lazy => 1,
    builder => '_build_is_het_chain',
);

has 'cluster_id' => (
    is => 'rw',
    isa => 'Num',
);

# Methods

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

around '_parse_ATOM_lines' => sub {

    my $orig = shift;
    my $self = shift;
    
    my @chain_ATOM_lines = ();
    
    foreach my $line ( $self->$orig() ) {
        if ( substr($line, 21, 1) eq $self->chain_id ){
            push(@chain_ATOM_lines, $line);
        }
    }
    
    croak "No ATOM lines were parsed with chainId " . $self->chain_id()
        . " for pdb " . $self->pdb_code()
            if ! @chain_ATOM_lines;

    return @chain_ATOM_lines;
};

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

sub _is_het_chain {
    my $self = shift;

    foreach my $atom (@{$self->atom_array()}) {
        return 0 if ! $atom->is_het_atom;
    }

    return 1;
}


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
    };
    
    if ($chainType eq 'Antigen') {
        return 0;
    }
    else {
        return $chainType;
    }
}

# This method runs kabatnum using pdb::kabatnum
sub kabatSequence {
    my $self = shift;

    my $kabatnum = pdb::kabatnum->new(input => $self);

    $kabatnum->sequenceChain();
}

# This method determines the atoms that are within a given distance threshold
# (default 4A) of the CDRs of the given antibody chains. The antibody chain
# CDR atoms must be labelled as such (i.e. $atom->is_CDR() == 1).
#
# Non-CDR contacting residues can be included by passing an additional, larger
# distance parameter. This distance parameter is used to find antigen residues
# close to the CDRs that are also in contact with any antibody residue. The
# larger the second distance threshold, the larger the epitope.
# INPUT: ref to array of antibody chains, two distance thresholds (optional)
# e.g. $chain->determineEpitope(\abChains, 4, 8);
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

    # Filter these residues so that only those residues that are also labelled
    # interface are included
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

# Labels atoms of chain as epitope according to the array of resids passed to it
sub labelEpitopeAtoms {
    my $self = shift;
    my @resids = @_;

    foreach my $resid (@resids) {
        # Remove chain-resSeq "." separator if present
        #$resid =~ s/\.//;

        foreach my $atom (values %{$self->resid_index->{$resid}}) {
            $atom->is_epitope(1);
        }
    }
}

# Returns array of residue resSeqs, where each residue has at least one atom
# labelled as epitope
sub getEpitopeResSeqs {
    my $self = shift;

    my %epitopeResSeqs = ();
    
    foreach my $atom (@{$self->atom_array()}) {
        $epitopeResSeqs{$atom->resSeq()} = 1
            if $atom->is_epitope();
    }
    return keys %epitopeResSeqs;
}

sub getInterfaceResidues {
    my $self = shift;
    my $otherChainsAref = shift;

    $self->read_ASA() if ! $self->has_read_ASA();

    my $atomAref
        = pdb::pdbFunctions::generateAtomAref($self, @{$otherChainsAref});
    
    # Calculate ASAb values of self + otherChains
    my $asurf = pdb::asurf64->new(input => $atomAref);
    my $atomSerialHref = $asurf->getOutput();
    my $complexResid2RelASAHref = $asurf->resid2RelASAHref();

    my @interfaceResidues = ();
    
    foreach my $resid (keys %{$self->resid_index()}) {
        
        my $complexRelASA = $complexResid2RelASAHref->{$resid}->{allAtoms};
        my $molRelASA = $self->resid2RelASAHref->{$resid}->{allAtoms};

        unless (defined $complexRelASA && defined $molRelASA) {
            croak "Undefined value for " . $self->pdb_code()
                . $self->chain_id() . " $resid!\n";
        }
        
        my $ASAdiff = $molRelASA - $complexRelASA;
        
        push(@interfaceResidues, $resid) if $ASAdiff >= 10;
    }
    return @interfaceResidues;
}


# This method returns an array of chain atoms labelled as CDR
# i.e. those atoms where $atom->is_CDR() is TRUE
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

# This method checks if the chain is in contact with the other chains passed.
# A tolerance can be used to set the minimum of residue-residue contacts that
# must exist for the chains to be considered contacting. Default = 1
# A threshold can also be set that specifies that maximum atom-atom distance
# that will quality residues as touching.
# e.g. $antigenChain->isInContact([$heavyChain, $lightChain], 5, 3)
# returns 1 if chain is in contact with any of the other chains passed. 
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

# This method processes an alignment string in order to label atoms with alnSeq
# If includeMissing is specified as true, residues missing in the structure
# will be included in the alignment
# example input:
#  $chain->processAlnStr(alnStr => $myStr, includeMissing => 1);
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
                print $self->get_sequence(include_missing_res => 1,
                                          return_type => 1),  "\n";
                print Dumper $self;
                exit;
            }
        }
        else {
            $atom->clusterSeq($cluster_id . "." . $atom->alnSeq());
        }
    }
    #map {$_->clusterSeq($cluster_id . $_->alnSeq())} @{$self->atom_array()};
}


# Automatically set arg chain_id
around [qw(get_sequence getFASTAStr)] => sub {

    my $orig = shift;
    my $self = shift;

    my %arg = @_;

    $arg{chain_id} = $self->chain_id();
    
    return $self->$orig(%arg);
    
};

# Automatically send chain id
around [qw(map_chainSeq2resSeq map_resSeq2chainSeq)] => sub {
    my $orig = shift;
    my $self = shift;

    my @arg = ( $self->chain_id() );
    
    return $self->$orig(@arg);    
};

__PACKAGE__->meta->make_immutable;


package patch;
use Moose;
use Moose::Util::TypeConstraints;
use types;
use Carp;

use pdb::PatchOrder_v3 qw(&PatchOrder);
use GLOBAL qq(&rm_trail);

extends 'pdb';

# Subtypes

# Attributes

has central_atom => (
    is => 'rw',
    isa => 'atom',
    required => 1,
);

has summary => (
    is => 'rw',
    isa => 'Str',
    lazy => 1,
    builder => '_build_summary',
);

has 'porder' => (
    is => 'ro',
    isa => 'Str',
    builder => 'run_PatchOrder',
    lazy => 1,
);

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

# Methods

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
        foreach my $atom (@{$self->atom_array()}, $self->central_atom()) {
            push(@atom_array, $atomSerialHref->{$atom->serial()});
        }
        $self->atom_array(\@atom_array);
    }
}


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

sub parseSummaryLine {
    my $summaryLine = shift;

    # example summary line : <patch G.409> G:335 G:397 G:398 G:407 G:408 G:409
    # parse all resids
    my @resids = $summaryLine =~ /(\w+[\.:]\w+)/g;
    
    # change any : separators to .
    map {s/:/./} @resids;
    
    my $centralResid = shift @resids;

    # Remove repeat of central residue
    @resids = grep {$_ ne $centralResid} @resids;
    
    return ($centralResid, @resids);
}


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


sub run_PatchOrder {    
    my $self = shift;

    my @ATOM_lines = ();
    
    foreach my $atom ( @{ $self->atom_array } ) {
        croak "atom $atom does not have an ATOM_line assigned"
            if ! $atom->ATOM_line;
        push( @ATOM_lines, $atom->ATOM_line );
    }

    croak "No ATOM_lines were assigned to array from patch $self"
        if ! @ATOM_lines;

    my $v = '';
    
    my $pOrder = PatchOrder_v3::PatchOrder( $v, $self->summary, @ATOM_lines );
    
    return $pOrder;
}

__PACKAGE__->meta->make_immutable;


package atom;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use Math::Trig;
use GLOBAL qw(&rm_trail);

use Carp;

### Attributes #################################################################

has 'ATOM_line' => (
    isa => 'Str',
    is  => 'rw',
);

has [qw(name resName element charge resSeq kabatSeq chothiaSeq ichothiaSeq
        alnSeq clusterSeq) ]
    => ( is => 'rw', isa => 'Str' );

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

my @labels = qw(is_het_atom is_terminal is_solvent is_CDR is_epitope);

foreach my $label (@labels) {
    has $label => (
        isa => 'Bool',
        is => 'rw',
        default => 0,
    );
}

### Methods ####################################################################

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
    
    my @ordered_attr = ( 'serial', 'name', 'altLoc', 'resName', 'chainID',
                       'resSeq','iCode', 'x', 'y', 'z', 'occupancy',
                         'tempFactor', 'element', 'charge' );
     
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

# Methods
sub BUILD {
    
    my $self = shift;

    my $ATOM_line = $self->ATOM_line();

    return if ! $ATOM_line;

    # Avoid substr complaining if there are missing columns at end of string
    $ATOM_line = pack ( "A81", $ATOM_line );
    my %record
        = ( ATOM => rm_trail( substr($ATOM_line, 0, 6) ),
            serial =>  rm_trail( substr($ATOM_line, 6, 5) ),
            name => rm_trail( substr($ATOM_line, 12, 4) ),
            altLoc => rm_trail( substr($ATOM_line, 16, 1) ),
            resName => rm_trail( substr($ATOM_line, 17, 3) ),
            chainID => rm_trail( substr($ATOM_line, 21, 1) ),
            resSeq => rm_trail( substr($ATOM_line, 22, 4) ),
            iCode => rm_trail( substr($ATOM_line, 26, 1) ),
            x => rm_trail( substr($ATOM_line, 30, 8) ),
            y => rm_trail( substr($ATOM_line, 38, 8) ),
            z => rm_trail( substr($ATOM_line, 46, 8) ),
            occupancy => rm_trail( substr($ATOM_line, 54, 6) ),
            tempFactor => rm_trail( substr($ATOM_line, 60, 6) ),
            element => rm_trail( substr($ATOM_line, 76, 2 ) ),
            charge => rm_trail( substr($ATOM_line, 78, 2) ),
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

sub _get_resid {
    my $self = shift;
    my $resid = join(".", ($self->chainID, $self->resSeq));

    return $resid;
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


__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

pdb - Perl extension for blah blah blah

=head1 SYNOPSIS

   use pdb;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb, 

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
 
