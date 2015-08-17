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

has 'string' => (
    is => 'rw',
    isa => 'SeqStr',
    required => 1,
    lazy => 1,
    builder => '_build_seqString'
);

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
        
sub getFASTAStr {
    my $self = shift;
    return ">" . $self->id() . "\n" . $self->string() . "\n";
}

1;
