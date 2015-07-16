package ARFF::Instance;
use Carp;
use strict;
use warnings;

use Moose;

### Attributes #################################################################
################################################################################

has '_attributeFromNameHref' => (
    is => 'rw',
    isa => 'HashRef[attribute]',
    default => sub { {} },
);

has 'idAttribute' => (
    is => 'rw',
    isa => 'attribute',
);

has 'classAttribute' => (
    is => 'rw',
    isa => 'attribute',
    predicate => 'has_classAttribute',
);

### Methods ####################################################################
################################################################################

sub getAttributeNames {
    my $self = shift;
    return keys %{$self->_attributeFromNameHref};
}

sub getAttributes {
    my $self = shift;
    return values %{$self->_attributeFromNameHref};
}

sub id {
    my $self = shift;
    my ($id)
        = map {$_->value} grep {$_->is_id()} $self->getAttributes;
    return $id;
}

sub class {
    my $self = shift;
    my ($class)
        = map {$_->value} grep {$_->is_class()} $self->getAttributes;
    return $class;
}

sub addAttributes {
    my $self  = shift;
    my @attrs = @_;

    map {$self->addAttribute($_)} @attrs;
}

sub addAttribute {
    my $self = shift;
    my $attr = shift;
    $self->_attributeFromNameHref->{$attr->name} = $attr;
}

sub removeAttribute {
    my $self = shift;
    my $attr = shift;
    $self->removeAttributeWithName($attr->name);
} 

sub removeAttributeWithName {
    my $self     = shift;
    my $attrName = shift;
    delete $self->_attributeFromNameHref->{$attrName};
}

sub setValueForAttributeWithName {
    my $self     = shift;
    my $value    = shift;
    my $attrName = shift;
    $self->_attributeFromNameHref->{$attrName}->value($value);
}

sub getValueForAttributeWithName {
    my $self     = shift;
    my $attrName = shift;
    return $self->_attributeFromNameHref->{$attrName}->value();
}

sub getHashOfAttributeNameToValue {
    my $self = shift;

    return map {$_ => $self->getValueForAttributeWithName($_)}
        $self->getAttributeNames();
}

sub string {
    my $self = shift;
    my @attrDescs = @_;

    return join(",", map {$self->getValueForAttributeWithName($_->name)} @attrDescs);
}

sub hasAttributeWithName {
    my $self     = shift;
    my $attrName = shift;

    return exists $self->_attributeFromNameHref->{$attrName};
}

sub copyValuesFrom {
    my $self      = shift;
    my $inst2Copy = shift;

    my @sharedAttrNames = $self->findAttributeNamesSharedWith($inst2Copy);
    my %valueForName
        = map {$_ => $inst2Copy->getValueForAttributeWithName($_)}
            @sharedAttrNames;

    while (my ($attrName, $value) = each %valueForName) {
        $self->setValueForAttributeWithName($value, $attrName);
    }
}

sub findAttributeNamesSharedWith {
    my $self      = shift;
    my $otherInst = shift;

    return grep {$self->hasAttributeWithName($_)}
        $otherInst->getAttributeNames();
}

__PACKAGE__->meta->make_immutable;
