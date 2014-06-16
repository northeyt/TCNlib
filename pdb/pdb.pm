package pdb;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use local::error;

use GLOBAL qw(&rm_trail &three2one_lc &is_int);

use Carp;
use Scalar::Util qw(looks_like_number);
use TryCatch;

use Math::VectorReal qw(:all);
use Math::MatrixReal;
use Math::Trig;

use pdb::xmas2pdb;
use pdb::getresol;
use pdb::rotate2pc;

# Subtypes

### Attributes

has 'pdb_code' => (
    isa => 'Str',
    is  => 'rw',
    default => '????', # Add a builder that attempts to determine a code
                       # from pdb_file, if has_pdb_file ?
);

for my $name ( 'pdb', 'xmas' ) {
    
    my $att_file = $name . '_file';
    my $att_data   = $name . '_data'  ;
    
    has $att_file => (
        isa => 'FileReadable',
        coerce => 1,
        # See types.pm for FileReadable coercion from ArrayRef[Str]
        is => 'rw',
        predicate => 'has_' . $att_file,
    );

    has $att_data => (
        isa => 'ArrayRef',
        is => 'rw',
        lazy => 1,
        builder => '_build_' . $att_data,
        predicate => 'has_' . $att_data,
    );
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
    
    my $predicate = "has_$att" . '_file';
    
    croak "Cannot get file data from $att file - no file specified"
        if ! $self->$predicate;

    my $file = $self->$att_file;
    
    open(my $fh, '<', $file) or die "Cannot open file '$file', $!";

    my @array = <$fh>;

    croak "No lines found in file $file" if ! @array;
   
    return \@array;
}

has 'atom_array' => (
    isa => 'ArrayRef[atom]',
    is  => 'ro',
    lazy => 1,
    builder => '_parse_atoms',
);

# Index is form of chain -> resSeq -> atom_name -> index for atom_array
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

# Index in form of resid -> atom_name -> index for atom_array
# Resid = chainResSeq
has 'resid_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_resid_index',
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
);

has 'resolution' => (
    is => 'rw',
    isa => 'Num',
    predicate => 'has_resolution',
);

has 'r_value' => (
    is => 'rw',
    isa => 'Num',
    predicate => 'has_r_value',
);


# Consume antigen role
with 'pdb::antigen';

### Methods

# Attempt to build experimental_method, resolution and r_factor from getresol
# object. Also check if pdb is multi-model
sub BUILD {
    my $self = shift;

    return if ! $self->has_pdb_file();

    my $getresol = pdb::getresol->new( pdb_file => $self->pdb_file );

    if ( ref $getresol->run() ne 'local::error' ){
        $self->experimental_method( $getresol->experimental_method() );
        $self->resolution( $getresol->resolution() );
        $self->r_value( $getresol->r_value() );
    }
 
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

    my $test_i = 0;
    
    foreach my $line (@ATOM_lines) {
        ++$test_i;
        
        if ( $line =~ /^TER/ ) {
            my($serial, $chainID) = _parse_ter($line);

            $ter{ $serial - 1 } = $chainID;
            next;
        }
        
        my $atom = atom->new( ATOM_line => $line );

        my $chain   = $atom->chainID();
        my $resSeq  = $atom->resSeq();
        my $resName = $atom->resName();

        if ( ! exists $test_unique{$chain} ) {
            $test_unique{$chain} = {};
        }
        
        if ( ! exists $test_unique{$chain}->{$resSeq} ) {
            $test_unique{$chain}->{$resSeq} = $resName;
        }
        else {
            my $recorded_resName = $test_unique{$chain}->{$resSeq};
            if ($recorded_resName ne $resName){
               
                my %multi = %{ $self->multi_resName_resid };
                $multi{ $atom->_get_resid } = 1;
                
                $self->multi_resName_resid( { %multi } );
            }
        }
        
        next if ( $h_clean && $atom->element eq 'H'
                      || $HETATM_clean && $atom->is_het_atom );
        
        if ( $aL_clean && $atom->has_altLoc ) {
            
            my $string
                =  $atom->name . $atom->resName . $atom->chainID
                 . $atom->resSeq;

            if ( exists $altLoc{$string} ) {
                push( @{ $altLoc{$string} }, $atom );
            }
            else {
                $altLoc{$string} = [ $atom ];
            }
            
        }
        else {   
            push(@atoms, $atom);   
        }
    }

    if ($aL_clean) {
        foreach my $arr ( values %altLoc ) {

            my @sorted
                = map { $_->[0] }
                    sort { $b->[1] <=> $a->[1] }
                        map { [ $_, $_->occupancy ] }
                            @{$arr};

            # Clear altLoc
            $sorted[0]->clear_altLoc();
            
            push( @atoms, $sorted[0] );
        }
    }

    # Order all atoms by chain, then serial
    my %chain = ();
    foreach my $atom (@atoms) {
        my $chainID = $atom->chainID();

        if ( ! exists $chain{$chainID} ) {
            $chain{$chainID} = { $atom->serial => $atom };
        }
        else {
            $chain{$chainID}->{ $atom->serial } = $atom;
        }
    }

    my @sorted_atoms = ();
    my @term_array   = ();
    
    foreach my $chainID (sort keys %chain) {
        my @chain_atoms = ();
        
        foreach my $atom ( sort { $a <=> $b } keys %{ $chain{$chainID} } ) {
            push( @chain_atoms, $chain{$chainID}->{$atom} );
        }

        my @chain_terminal_index = ();
        
        # Label terminal atoms
        for my $i ( 0 .. @chain_atoms - 1 ) {
            my $atom = $chain_atoms[$i];
            my $serial = $atom->serial();
            if (! defined $serial) {
                croak "Undefined $serial!";
            }
            if ( exists $ter{$serial} && $atom->has_chainID
                     && $atom->chainID eq $ter{$serial} ){
                $atom->is_terminal(1);
                
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

sub _build_resid_index {
    my $self = shift;

    my %hash = ();

    foreach my $chain (keys %{ $self->atom_index } ) {
        my %chain_h = %{ $self->atom_index->{$chain} };
        
        foreach my $resSeq (keys %chain_h) {
            my $resid = $chain . $resSeq;
            $hash{$resid} = $chain_h{$resSeq};
        }
    }

    croak "Nothing indexed while attempting to index by resid"
        if ! %hash;

    return \%hash;
}

sub get_sequence {
    my $self = shift;

    croak "pdb: " . $self->pdb_code()
        . " can't get_sequence for pdb containing multi-resName residues: "
        . Dumper $self->multi_resName_resid()
            if %{ $self->multi_resName_resid };
    
    my $USAGE
        = 'get_sequence( chain_id => (chain_id), return_type => ( 1 | 3 ) )';
    
    my %arg = @_;

    # Check args are okay
    if (   ! (   exists  $arg{return_type} && exists  $arg{chain_id} )
        || ! (   defined $arg{return_type} && defined $arg{chain_id} )
        || ! ( $arg{return_type} eq '1' || $arg{return_type} eq '3'  ) ) {
        croak $USAGE;
    }

    my %chain_resid_index = ();

    
    # Filter out residues that are not from specified chain
    foreach my $resid ( keys %{ $self->resid_index() } ){
        next if substr($resid, 0, 1) ne $arg{chain_id};
        $chain_resid_index{$resid} = $self->resid_index->{$resid};
    }
    
    #  Sort values of each resid of resid index by resSeq
    my @residue_atoms
        = map { [ values %{ $self->resid_index->{$_} } ]  }
            sort { substr($a, 1) <=> substr($b, 1) }
                keys %chain_resid_index;

    my @residues = ();
    
    # For each residue, get resName. Avoid solvent residues
    foreach my $atom_array (@residue_atoms) {

        my $atom = $atom_array->[0];
        next if $atom->is_solvent();
        
        my $res_name = $atom->resName();
        
        push(@residues, $res_name);
    }
    
    if ( $arg{return_type} == 1 ) {
        for my $i ( 0 .. @residues - 1 ) {
            my $onelc = eval { three2one_lc( $residues[$i] ) };
            if ($@) {
                $onelc = 'X';
            }
            $residues[$i] = $onelc;
        }
    }
    return @residues;
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
        foreach my $resSeq (sort {compare_resSeqs($a,$b)} keys % {$resSeqHref}){
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

            
# This functions compares two resSeqs to order them, like inbuilt cmp
# -1 is returned if $rsA should before $rsB
#  1 is returned is $rsA should come after $rsB
sub compare_resSeqs {
    my($rsA, $rsB) = @_;

    if (looks_like_number($rsA) && looks_like_number($rsB)) {
        # Both are resSeqs are numbers and can be sorted with a simple <=>
        return $rsA <=> $rsB;
    }
    else {
        my ($rsA_num) = $rsA =~ /(\d+)+/g;
        my ($rsB_num) = $rsB =~ /(\d+)+/g;
        
        my ($rsA_suffix) = $rsA =~ /([A-Z]+)$/g;
        my ($rsB_suffix) = $rsB =~ /([A-Z]+)$/g;

        # If resSeq does not have a suffix, set var to empty string
        # This will ensure that cmp orders a resSeq with no prefix
        # before those with prefixes
        foreach my $suffref (\$rsA_suffix, \$rsB_suffix) {
            ${$suffref} = "" if ! ${$suffref};
        }

        # Sort first by resSeq number then suffix
        return ($rsA_num <=> $rsB_num || $rsA_suffix cmp $rsB_suffix);
    }
}

sub read_ASA {

    my $self = shift;

    my $xmas2pdb;
    
    # Use xmas2pdb object if it has been passed, otherwise create one
    if (@_) {
        $xmas2pdb = shift;
       
        croak "read_ASA must be passed an xmas2pdb object"
            if ref $xmas2pdb ne 'xmas2pdb';
    }
    else {
        croak "read_ASA: pdb object must have an xmas file\n"
            if ! $self->has_xmas_file();

        my $form = ref $self eq 'pdb' ? 'multimer' : 'monomer' ;
        
        my %x2p_arg = ( xmas_file     => $self->xmas_file(),
                        form          => $form, );

        $xmas2pdb = xmas2pdb->new(%x2p_arg);
    }
    
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

    $self->has_read_ASA(1);
    
    return \@errors;
};

sub patch_centres {

    my $self = shift;
    
    my %arg = @_;

    my %opt = (
        ASA_threshold => [ 25, 'num' ],
        #makepatch     => [ '', 'required.object', 'makepatch' ],
    );

    # Get opts
    foreach my $value (keys %opt) {
        if ( ! exists $arg{$value} ) {
            croak "$value arg must be specified"
                  if $opt{$value}->[1] =~ /required/; 

            # Set to default
            $arg{$value} = $opt{$value}->[0];
        }
        elsif ( $opt{$value}->[1] =~ /object/ ){
            my $ref = ref $arg{$value};
            croak "$value arg must be a $ref object"
                if $ref ne $opt{$value}->[2];
        }
        elsif ( $opt{$value}->[1] =~ /int/ ) {
            croak "$value arg must be an int"
                if int $arg{$value} ne $arg{$value};
        }
           
    }

    my @central_atoms = ();

    my @errors = ();
    
    foreach my $chain_h ( keys %{ $self->atom_index } ) {
        my %chain_h = %{ $self->atom_index->{$chain_h} };
        
        foreach my $resSeq ( keys %chain_h ) {
            my %atom_h = %{ $chain_h{$resSeq} };

            my @atoms = ();

            my $CA_flag = 0;
            
            foreach my $atom_name (keys %atom_h) {

                $CA_flag = 1 if $atom_name eq 'CA';
                
                my $atom = $atom_h{$atom_name};

                next if $atom->is_solvent();
                
                push( @atoms, $atom );

            }

            # Patch centre residues must have a CA atom to be run through
            # makepatch
            
            next if ! $CA_flag;

            my $ret = $self->_is_patch_centre( $arg{ASA_threshold},
                                                       'ASAc',
                                                       @atoms );
            if ( ref $ret eq 'local::error' ) {
                my $message
                    =  "Could not determine if residue " . $resSeq
                      . " is valid patch center: " . $ret->message();

                my $error = local::error->new(
                    message => $message,
                    type => 'no_patch_centre_value_for_residue',
                    data => { residue => $resSeq,
                              chain_id => $chain_h,
                              atoms => [@atoms],
                          },
                    parent => $ret,
                );

                push(@errors, $error);
            }
            elsif ($ret ne "-1" ) {
                push( @central_atoms, $ret);
            }
        }
    }
    return(\@errors, \@central_atoms);
}

sub _is_patch_centre {

    my $self = shift;
    my $threshold = shift;
    my $attribute = shift;
    my @atoms = @_;

    foreach my $atom (@atoms) {
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
    
    @atoms
        = sort {$b->$attribute() <=> $a->$attribute()} @atoms;

    my $total = 0;

    map { $total += $_->$attribute() } @atoms;

    if ($total >= $threshold) {
        return $atoms[0];
    }

    return -1;
}

sub highestASA {
    my $self  = shift;
    my $resid = shift or croak "highestASA must be passed a resid";

    my $ASA_type
        =  ref $self eq 'pdb'   ? 'ASAc'
         : ref $self eq 'chain' ? 'ASAm'
         : '' ;

    croak "highestASA: something went wrong assigning ASA type"
        if ! $ASA_type;

    croak "resid '$resid' was not found in resid index"
        if ! exists $self->resid_index->{$resid};

    my @atoms = values %{ $self->resid_index->{$resid} };
    
    foreach my $atom (@atoms) {
        my $predicate = "has_$ASA_type";
        croak 'atom ' . $atom->serial() . " has no $ASA_type value"
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
        my $chain = chain->new(chain_id => $chain_id,
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

sub getAbPairs {
    my $self = shift;

    my @chains = $self->create_chains();

    my %chainTypes = (Heavy => [], Light => [],
                      antigen => [], scFv => []);

    # Hash to keep track of which chains have be paired
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
    foreach my $chain ( @{$chainTypes{heavy}}, @{$chainTypes{light}} ) {
        if ( ! $paired{$chain->chain_id()} ) {
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

# This method checks if the chain is an antibody variable chain, or antigen
# Is the chain is antibody variable, the type is returned i.e. Light or Heavy
sub isAbVariable {
    my $self = shift;

    my $idabchain = pdb::idabchain->new(input => $self);

    my $chainType = $idabchain->chainIs();
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
# CDR atoms must be labelled as such (i.e. $atom->is_CDR() == 1)
# INPUT: ref to array of antibody chains, distance threshold (optional)
# e.g. $chain->determineEpitope(\abChains, 4);
sub determineEpitope {
    my($self, $abChainsAref, $distance) = @_;

    # Default distance if not supplied by user
    $distance = 4 if ! $distance;

    # Get array of CDR atoms from antibody chains
    my @CDRAtoms = ();
    foreach my $chain (@{$abChainsAref}) {
        push(@CDRAtoms, $chain->getCDRAtoms());
    }

    # Create array of CDR atoms and antibody chain atoms
    my @allAtoms = (@CDRAtoms, @{$self->atom_array()});

    # Test for contacts
    my $chainContacts = pdb::chaincontacts->new(input => \@allAtoms);
    my $cContactResult = $chainContacts->getOutput();

    # Get resids of antigen residues that are in contact with ab chains
    my $residAref
        = $cContactResult->chain2chainContacts($abChainsAref, [$self]);

    # Label atoms from resids as epitope
    $self->labelEpitopeAtoms(@{$residAref});
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
sub determineEpitope2 {
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
    my $residAref = _determineEpitopeResids($abChainsAref, $self, $exContactResult,
                                      $allcContactResult);
    
    # Label atoms from resids as epitope
    $self->labelEpitopeAtoms(@{$residAref});
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
        = grep ($residHrefs[0]->{$_}, keys %{$residHrefs[1]});

    return \@epitopeResids;
}

# Labels atoms of chain as epitope according to the array of resids passed to it
sub labelEpitopeAtoms {
    my $self = shift;
    my @resids = @_;

    foreach my $resid (@resids) {
        # Remove chain-resSeq "." separator if present
        $resid =~ s/\.//;

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
# e.g. $antigenChain->isInContact([$heavyChain, $lightChain], 5)
# returns 1 if chain is in contact with any of the other chains passed. 
sub isInContact {
    my $self = shift;
    my($chainAref, $contactMin) = @_;

    # Set default is contactMin not passed
    $contactMin = 1 if ! $contactMin;

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

    print @allAtoms;
    
    my $cContacts = pdb::chaincontacts->new(input => \@allAtoms);
    my $contResults = $cContacts->getOutput();

    my $contactAref = $contResults->chain2chainContacts([$self], $chainAref);

    use Data::Dumper;
    print "DEBUG\n";
    print Dumper $contactAref;
    
    if (scalar @{$contactAref} >= $contactMin) {
        return 1;
    }
    else {
        return 0;
    }
}


### around MODIFIERS
# Modify _is_patch_centre to assess on monomer ASAm
around '_is_patch_centre' => sub {
    
    my $orig = shift;
    my $self = shift;
    
    my @arg = @_;
    
    $arg[1] = 'ASAm';
 
    return $self->$orig(@arg);
    
};

# Automatically set arg chain_id
around 'get_sequence' => sub {

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
);

has 'porder' => (
    is => 'ro',
    isa => 'Str',
    builder => 'run_PatchOrder',
    lazy => 1,
);

# Methods

# Allow a makepatch object to be passed directly to new method
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;

    if ( ref $_[0] eq 'makepatch' ) {
        my $makepatch = $_[0];
        
        if ( ref $makepatch->output()->[0] eq 'local::error' ) {
            croak $makepatch->output()->[0];
        }
        
        my %arg
            = ( central_atom => $makepatch->central_atom,
                pdb_data     => $makepatch->output,
                summary      => rm_trail( $makepatch->output->[-1] ),
                pdb_code     => $makepatch->pdb_code,
            );
        
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

has [ 'name', 'resName', 'element', 'charge', 'resSeq',
      'kabatSeq', 'chothiaSeq', 'ichothiaSeq' ]
    => ( is => 'rw', isa => 'Str' );

foreach my $name ( 'altLoc', 'chainID', 'iCode' ) {
    my $predicate = 'has_'   . $name;
    my $clearer   = 'clear_' . $name;
    
    has $name => ( is => 'rw',
                   isa => 'Str',
                   predicate => $predicate,
                   clearer => $clearer ); 
}

has [ 'serial'] => ( is => 'rw', => isa => 'Int' );

foreach my $name ( 'radius', 'ASAm', 'ASAc', 'x', 'y', 'z', 'occupancy',
                'tempFactor', ) {
    my $predicate = 'has_' . $name;
    
    has $name => ( is => 'rw', isa => 'Num', predicate => $predicate );
}

foreach my $name ('rASAm', 'rASAc') {
    my $predicate = 'has_' . $name;
    my $builder = '_build_' . $name;
    has $name => ( is => 'rw', isa => 'Num', predicate => $predicate,
                   builder => $builder, lazy => 1,); 
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

sub _build_rASAm {
    my $self = shift;
    croak "Atom has no ASAm value, cannot calculate rASAm"
        if ! $self->has_ASAm();

    return $self->_build_rASA($self->ASAm());
}

sub _build_rASAc {
    my $self = shift;
    croak "Atom has no ASAc value, cannot calculate rASAc"
        if ! $self->has_ASAc();

    return $self->_build_rASA($self->ASAc());
}

sub _build_rASA {
    my $self = shift;
    my $ASA  = shift;
    
    croak "Atom has no radius, cannot calculate a rASA"
        if ! $self->has_radius();

    my $r = $self->radius();
    
    return $ASA / (4 * pi * $r**2); 
}

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
        next if $record{$value} eq '' ;
        $self->$value( $record{$value} );
    }
}

sub _get_resid {
    my $self = shift;
    my $resid = $self->chainID . $self->resSeq;

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
 
