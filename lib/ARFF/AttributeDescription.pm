package ARFF::AttributeDescription;
use Moose;

### Attributes #################################################################
################################################################################

has 'name'  => (
    is => 'rw',
    isa => 'Str',
);

has 'type'  => (
    is => 'rw',
    isa => 'Str',
);

has '_prevAttrDesc' => (
    is => 'rw',
    isa => 'ARFF::AttributeDescription',
    predicate => 'has_prevAttrDesc',
    clearer => 'clear_prevAttrDesc',
);

has '_nextAttrDesc' => (
    is => 'rw',
    isa => 'ARFF::AttributeDescription',
    predicate => 'has_nextAttrDesc',
    clearer => 'clear_nextAttrDesc',
);

has 'is_id' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

has 'is_class' => (
    is => 'rw',
    isa => 'Bool',
    default => 0
);

### Builders ###################################################################
################################################################################

around 'BUILDARGS' => sub {
    my $orig  = shift;
    my $class = shift;

    if (@_ == 1) {
        return $class->$orig(_parseAttributeLine(@_));
    }
    else {
        $class->$orig(@_);
    }
};

sub _parseAttributeLine {
    my $arg = shift;

    my ($name, $type) = split(" ", $arg);
    return (name => $name, type => $type);
}

### Methods ####################################################################
################################################################################

sub string {
    my $self = shift;
    return join(" ", '@attribute', $self->name(), $self->type);
}

__PACKAGE__->meta->make_immutable;

1;
