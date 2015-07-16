package ARFF::Attribute;
use Moose;
use Moose::Util::TypeConstraints;

has 'value' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
);

coerce 'ARFF::AttributeDescription',
    from 'Str',
    via {ARFF::AttributeDescription->new($_)};

has 'description' => (
    is => 'rw',
    isa => 'ARFF::AttributeDescription',
    required => 1,
    coerce => 1,
    handles => [qw(name type is_id is_class)],
);

sub string {
    my $self = shift;
    return join(" ", '@attribute', $self->name(), $self->type);
}

__PACKAGE__->meta->make_immutable;

1;
