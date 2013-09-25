package pdb;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use GLOBAL qw(&rm_trail);
use search_hash;

use Data::Dumper;
use Carp;

# Subtypes

### Attributes

has 'pdb_code' => (
    isa => 'Str',
    is  => 'rw',
    default => '????', # Add a builder that attempts to determine a code
                       # from pdb_file, if has_pdb_file ?
);

has 'pdb_file' => (
    isa => 'FileReadable',
    coerce => 1, # See types.pm for FileReadable coercion from ArrayRef[Str]
    is => 'rw',
    required => 1,
    predicate => 'has_pdb_file',
);

has 'xmas_file' => (
    isa => 'FileReadable',
    is => 'rw',
    predicate => 'has_xmas_file',
);


has 'atom_array' => (
    isa => 'ArrayRef[atom]',
    is  => 'ro',
    lazy => 1,
    builder => '_parse_atoms',
);

# Index is form of chain -> resSeq -> name -> index for atom_array
has 'atom_index' => (
    isa => 'HashRef',
    is => 'ro',
    lazy => 1,
    builder => '_build_atom_index',
);

# Consume antigen role
with 'pdb::antigen';

### Methods

sub _parse_atoms {
    my $self = shift;

    croak "pdb_file attribute must be assigned before atoms can be parsed"
        if ! $self->has_pdb_file();
    
    my @ATOM_lines = $self->_parse_ATOM_lines();
  
    my @atoms = ();
    
    foreach my $line (@ATOM_lines) {
        my $atom = atom->new( ATOM_line => $line );
        push(@atoms, $atom);
    }

    croak "No atom objects were created" if ! @atoms;

    return [ @atoms ];
}


sub _parse_ATOM_lines {
    
    my $self = shift;

    my $file = $self->pdb_file;
    
    open(my $fh, '<', $file)
        or die "Cannot open file $file to parse ATOM lines";
  
    my @ATOM_lines = ();

    while ( my $line = <$fh> ) {
        if ($line =~ /^ATOM /) {
            push(@ATOM_lines, $line);
        }
    }

    close $fh;
    
    croak "No ATOM lines parsed from pdb file $file"
        if ! @ATOM_lines;

    return @ATOM_lines;
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

sub read_ASA {

    my $self = shift;
    my $xmas2pdb = shift;
       
    croak "read_ASA must be passed an xmas2pdb object"
        if ref $xmas2pdb ne 'xmas2pdb';

    my $form = $xmas2pdb->form();

    my $attribute = 'ASA' . ( $form eq 'monomer' ? 'm' : 'c' ) ;

    my  %ASAs  = ();
    my  %radii = ();
    
    foreach my $line ( @{ $xmas2pdb->output() } ){
        next unless $line =~ /^ATOM /;

        my $serial = rm_trail( substr($line, 6, 5) );
        my $radius = rm_trail( substr($line, 54, 5) );
        my $ASA    = rm_trail( substr($line, 60, 6) );

        $ASAs{$serial} = $ASA;
        $radii{$serial}= $radius;
    }

    croak "Nothing parsed from xmas2pdb output!" if ! ( %ASAs && %radii );
    
    foreach my $atom ( @{ $self->atom_array() } ) {
        my $serial = $atom->serial();
        croak "Nothing parsed from xmas2pdbfor atom " . $atom->serial()
            if ! exists $ASAs{$serial};
        
        $atom->radius( $radii{$serial} );
        $atom->$attribute( $ASAs{$serial} );
    }
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
    
    foreach my $chain_h ( keys %{ $self->atom_index } ) {
        my %chain_h = %{ $self->atom_index->{$chain_h} };
        
        foreach my $resSeq ( keys %chain_h ) {
            my %atom_h = %{ $chain_h{$resSeq} };

            my @atom_indices = ();
            foreach my $atom_name (keys %atom_h) {

                my $index = $atom_h{$atom_name};
                
                push( @atom_indices, $index );

            }
            
            my $central_atom
                = $self->_is_patch_centre( $arg{ASA_threshold},
                                           'ASAc',
                                           @atom_indices );
            
            if ( $central_atom != -1 ) {
                push( @central_atoms, $central_atom);
            }
            
        }
    }

    @central_atoms ? return @central_atoms
                   : croak "No patch centres found!";
    
}

sub _is_patch_centre {
   
    my $self = shift;
    my $threshold = shift;
    my $attribute = shift;
    my @indices = @_;

    foreach my $index (@indices) {
        my $predicate = "has_$attribute";
        croak "$attribute has not been set"
            if ! $self->atom_array->[$index]->$predicate();
        
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

package chain;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use Carp;

extends 'pdb';

# Subtypes
subtype 'ValidId',
    as 'Str',
    where { $_ =~ /^[A-Z]$/i },
    message { "$_ is not a valid Chain Id" };

# Attributes

has 'chain_id' => (
    isa => 'ValidId',
    is => 'rw',
    required => 1,
);

# Methods
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
        my %arg
            = ( central_atom => $makepatch->central_atom,
                pdb_file     => $makepatch->output,
                summary      => rm_trail( $makepatch->output->[-1] )
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
        if ( substr($line, 60, 6) =~ /1\.00/ ){
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

package atom;

use Moose;
use Moose::Util::TypeConstraints;
use types;

use Carp;

# Subtypes

subtype 'Character',
    as 'Str',
    where { $_ =~ /\w{1}/ },
    message { "$_ is not a valid single-char string!" };

# Attributes

has 'ATOM_line' => (
    isa => 'Str',
    is  => 'rw',
);

has [ 'name', 'resName', 'charge' ] => ( is => 'rw', isa => 'Str' );

has [ 'altLoc', 'chainID', 'iCode', ] => ( is => 'rw', isa => 'Character' );

has [ 'serial', 'resSeq', ] => ( is => 'rw', => isa => 'Int' );

foreach my $name ( 'radius', 'ASAm', 'ASAc' ) {
    my $predicate = 'has_' . $name;
    
    has $name => ( is => 'rw', isa => 'Num', predicate => $predicate );
}

has [ 'x', 'y', 'z', 'occupancy', 'tempFactor' ]
    => ( is => 'rw', isa => 'Num');

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
            charge => rm_trail( substr($ATOM_line, 78, 2) ),
        );

    croak "It looks like you're trying to parse a non-ATOM line: $ATOM_line"
        if $record{ATOM} ne 'ATOM';

    delete $record{ATOM};
    
    foreach my $value (keys %record) {
        next if $record{$value} eq '' ;
        $self->$value( $record{$value} );
    }
}

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
