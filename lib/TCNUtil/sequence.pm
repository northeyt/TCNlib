package sequence;
use Moose;
use Moose::Util::TypeConstraints;
use Carp;
use TCNUtil::seqStr;
    
subtype 'SeqStr' => (
    as 'Str',
    where {$_ =~ /^[A-Z-]+$/i}
);

has 'rawInputString' => (
    is => 'rw',
    isa => 'Str',
);

has 'id' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
    lazy => 1,
    builder => '_build_id'
);

has 'description' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
    lazy => 1,
    default => 'No Description Supplied!',
);

has 'string' => (
    is => 'rw',
    isa => 'SeqStr',
    required => 1,
    lazy => 1,
    builder => '_build_seqString'
);

has 'seqType' => (
    is => 'rw',
    isa => enum([ qw(protein    protein-frag dna-linear dna-circular
                     rna-linear rna-fragment unknown) ]),
    default => 'protein',
);

has '_PIRCodeForSequenceType' => (
    isa => 'HashRef',
    is  => 'rw',
    lazy => 1,
    builder => '_buildPIRCodeForSequenceType',
);

sub getFASTAStr {
    my $self = shift;
    return ">" . $self->id() . "\n" . $self->string() . "\n";
}

sub getPIRStr {
    my $self = shift;
    return ">" . $self->getPIRCode()  . ";" . $self->id() . "\n"
        . $self->description() . "\n" . $self->string()   . "*\n";
}

sub getPIRCode {
    my $self = shift;
    return $self->_PIRCodeForSequenceType($self->seqType());
}

# Allows lazy object creation: SeqStringAdapter->($myString)
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    if (@_ == 1) {
        return $class->$orig(rawInputString => $_[0]);
    }
    else {
        return $class->$orig(@_);
    }
};

sub _build_id {
    my $self = shift;

    my $id = seqStr::parseIDFromFASTA($self->rawInputString);

    defined $id ? return $id
        : croak "no id parsed from input string, " . $self->rawInputString();
}

sub _build_seqString {
    my $self = shift;

    my $seq = seqStr::parseSeqFromFASTA($self->rawInputString);
    
    defined $seq ? return $seq
        : croak "no seq parsed from input string, " . $self->inputString();
}

sub _buildPIRCodeForSequenceType {
    return {'protein'    => 'P1', 'protein-frag' => 'F1',
            'dna-linear' => 'DL', 'dna-circular' => 'DC',
            'rna-linear' => 'RL', 'rna-circular' => 'RC',
            'unknown'    => 'XX'};
}

1;
