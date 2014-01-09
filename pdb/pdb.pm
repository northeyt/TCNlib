package pdb;

use FileHandle;
use Tie::File;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use local::error;

use GLOBAL qw(&rm_trail);

use Data::Dumper;
use Carp;

use pdb::xmas2pdb;
use pdb::getresol;

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
    
    my $predicate = $self . '->'. $att . '_file';
    
    croak "Cannot get file data from $att file - no file specified"
        if ! $predicate;

    my $file = $self->$att_file;
    
    open(my $fh, '<', $file) or die "Cannot open file $file, $!";

    my @array = ();

    tie @array, 'Tie::File', $file, mode => O_RDONLY
        or die "Cannot tie array to file $file, $!"; 

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
    isa => 'ArrayRef[Int]',
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
    
    for (my $i = 0 ; $i <  @{ $self->atom_array() } ; ++$i ) {
        
        my $atom = $self->atom_array->[$i];
        my $index = $i;
        
        my $chainID = $atom->chainID();
        my $resSeq  = $atom->resSeq();
        my $name    = $atom->name();

        if ( ! exists $hash{$chainID} )  {
            $hash{$chainID} = {};
        };

        if ( ! exists $hash{$chainID}->{$resSeq} ) {
            $hash{$chainID}->{$resSeq} = {};
        };

        $hash{$chainID}->{$resSeq}->{$name} = $index;
    }

    return { %hash };
    
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

    return {%hash};
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

            my @atom_indices = ();
            foreach my $atom_name (keys %atom_h) {

                my $index = $atom_h{$atom_name};
                my $atom  = $self->atom_array->[$index];

                next if $atom->is_solvent();
                
                push( @atom_indices, $index );

            }


            my $ret = $self->_is_patch_centre( $arg{ASA_threshold},
                                                       'ASAc',
                                                       @atom_indices );
            if ( ref $ret eq 'local::error' ) {
                my $message
                    =  "Could not determine if residue " . $resSeq
                      . " is valid patch center: " . $ret->message();

                my $error = local::error->new(
                    message => $message,
                    type => 'no_patch_centre_value_for_residue',
                    data => { residue => $resSeq,
                              chain_id => $chain_h,
                              atoms =>
                                  [ map { $self->atom_array->[$_] }
                                      @atom_indices ],
                          },
                    parent => $ret,
                );

                push(@errors, $error);
            }
            elsif ( $ret != -1 ) {
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
    my @indices = @_;

    foreach my $index (@indices) {
        my $predicate = "has_$attribute";
        
            if( ! $self->atom_array->[$index]->$predicate() ){
                my $message
                    =  "Could not determine patch centre status because "
                      . "atom " . $self->atom_array->[$index]->serial
                      . " has no $attribute value";
                
                my $error = local::error->new(
                    message => $message,
                    type => "no_atom_$attribute" . "_set",
                    data => { atom => $self->atom_array->[$index], },
                );

                return $error;
            }        
    }
    
    @indices
        = sort {    $self->atom_array->[$b]->$attribute()
                <=> $self->atom_array->[$a]->$attribute() } @indices;

    my $total = 0;

    map { $total += $self->atom_array->[$_]->$attribute() } @indices;

    if ($total >= $threshold) {
        return $indices[0];
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

    my @atoms
        = map { $self->atom_array->[$_] }
            values %{ $self->resid_index->{$resid} };
    
    foreach my $atom (@atoms) {
        my $predicate = "has_$ASA_type";
        croak 'atom ' . $atom->serial() . " has no $ASA_type value"
            if ! $atom->$predicate();
    }

    my @sorted = sort { $b->$ASA_type <=> $a->$ASA_type  } @atoms;

    my $top_ASA_atom = $sorted[0];

    return $top_ASA_atom;
}


__PACKAGE__->meta->make_immutable;


package chain;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use Carp;

use pdb::pdbsws;

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


# Methods

# Returns number of non-solvent residues found in chain
sub _build_chain_length {
    my $self = shift;

    my %resid_h = %{ $self->resid_index };

    my $count = 0;

    foreach my $resid (keys %resid_h) {
        my @atom_index = values %{ $resid_h{$resid} };
        
        ++$count if ! $self->atom_array->[ $atom_index[0] ]->is_solvent();
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
    
    croak "No ATOM lines were parsed with chainId " . $self->chain_id
        if ! @chain_ATOM_lines;

    return @chain_ATOM_lines;
};

# Modify _is_patch_centre to assess on monomer ASAm
around '_is_patch_centre' => sub {
    
    my $orig = shift;
    my $self = shift;
    
    my @arg = @_;
    
    $arg[1] = 'ASAm';
 
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

use Carp;

# Subtypes

subtype 'Character',
    as 'Str',
    where { $_ =~ /[^\s]{1}/ },
    message { "$_ is not a valid single-char string!" };

# Attributes

has 'ATOM_line' => (
    isa => 'Str',
    is  => 'rw',
);



has [ 'name', 'resName', 'element', 'charge' ]
    => ( is => 'rw', isa => 'Str' );

foreach my $name ( 'altLoc', 'chainID', 'iCode' ) {
    my $predicate = 'has_'   . $name;
    my $clearer   = 'clear_' . $name;
    
    has $name => ( is => 'rw',
                   isa => 'Character',
                   predicate => $predicate,
                   clearer => $clearer ); 
}

has [ 'serial', 'resSeq', ] => ( is => 'rw', => isa => 'Int' );

foreach my $name ( 'radius', 'ASAm', 'ASAc', 'x', 'y', 'z', 'occupancy',
                'tempFactor' ) {
    my $predicate = 'has_' . $name;
    
    has $name => ( is => 'rw', isa => 'Num', predicate => $predicate );
}

has 'resid' => (
    is => 'ro',
    isa => 'Str',
    lazy => 1,
    builder => '_get_resid',
);

my @labels = ( 'is_het_atom', 'is_terminal', 'is_solvent' );

foreach my $label (@labels) {
    has $label => (
        isa => 'Bool',
        is => 'rw',
        default => 0,
    );
}

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

    use GLOBAL qw(&rm_trail);
    
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
