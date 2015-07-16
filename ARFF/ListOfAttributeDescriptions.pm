package ARFF::ListOfAttributeDescriptions;
use Moose;
use Carp;
use List::Util qw(all);

### Attributes #################################################################
################################################################################

has '_root' => (
    is       => 'rw',
    isa      => 'ARFF::AttributeDescription',
    required => 1,
);

has '_head' => (
    is       => 'rw',
    isa      => 'ARFF::AttributeDescription',
    required => 1,
);


### Builders ###################################################################
################################################################################

around 'BUILDARGS' => sub {
    my $orig  = shift;
    my $class = shift;
    
    if (_allArgsAreAttributeDescriptions(@_)) {
        _linkAttributeDescriptions(@_);
        return $class->$orig(_root => _getRootForList(@_),
                             _head => _getHeadForList(@_));
    }
    else {
        return $class->$orig(@_);
    }
};

sub _getRootForList {
    my @attrDescs = @_;
    my $root = ARFF::AttributeDescription->new();
    _linkAttributeDescriptionPair($root, $attrDescs[0]);
    return $root;
}

sub _getHeadForList {
    my @attrDescs = @_;
    my $head = ARFF::AttributeDescription->new();
    _linkAttributeDescriptionPair($attrDescs[-1], $head);
    return $head;
}

sub _allArgsAreAttributeDescriptions {
    my @args = @_;

    return all {ref $_ eq 'ARFF::AttributeDescription'} @args;
}

### Functions ##################################################################
################################################################################

sub _linkAttributeDescriptions {
    my @attrDescs = @_;

    for (my $i = 1 ; $i < @attrDescs ; ++$i) {
        _linkAttributeDescriptionPair($attrDescs[$i - 1], $attrDescs[$i]);
    }
}

sub _linkAttributeDescriptionPair {
    my @pair = @_;
    
    $pair[0]->_nextAttrDesc($pair[1]);
    $pair[1]->_prevAttrDesc($pair[0]);
}

### Methods ####################################################################
################################################################################

sub addDescription {
    my $self       = shift;
    my $attrDesc   = shift;
    my $addAtIndex = shift;
    
    # We either get the attrDesc at the index specified or get the head
    # of the list (i.e. the new attrDesc will be added at the end, before
    # the head).
    my $currentAttrDescAtIndex
        = defined $addAtIndex ? $self->_attributeAtIndex($addAtIndex)
        : $self->_head();
    
    # from p-C-n to p-A-C-n
    _linkAttributeDescriptionPair($currentAttrDescAtIndex->_prevAttrDesc,
                                  $attrDesc)
        if $currentAttrDescAtIndex->has_prevAttrDesc();
    
    _linkAttributeDescriptionPair($attrDesc, $currentAttrDescAtIndex);
}

sub removeDescriptionWithName {
    my $self     = shift;
    my $attrName = shift;

    $self->removeDescription($self->getDescriptionWithName($attrName));
}

sub removeDescription {
    my $self = shift;
    my $attrDesc = shift;

    # from p-A-n to p-n
    _linkAttributeDescriptionPair($attrDesc->_prevAttrDesc,
                                  $attrDesc->_nextAttrDesc)
        if $attrDesc->has_prevAttrDesc && $attrDesc->has_nextAttrDesc;
}

sub count {
    my $self = shift;
    return scalar $self->getDescriptions();
}

sub getDescriptions {
    my $self = shift;

    my @descs = ();
    my $currentAttrDesc = $self->_root();

    while (my $nextAttrDesc = $self->_getNextDesc($currentAttrDesc)) {
        push(@descs, $nextAttrDesc);
        $currentAttrDesc = $nextAttrDesc;
    }
    return @descs;
}

sub cmpDescriptions {
    my $self  = shift;
    my $descA = shift;
    my $descB = shift;

    return 0 if $descA eq $descB;

    my $currentDesc = $descA;
    while (my $prevDesc = $self->_getPrevDesc($currentDesc)) {
        if ($prevDesc eq $descB) {
            # descB is before descA, so return 1
            return 1;
        }
        $currentDesc = $prevDesc;
    }

    $currentDesc = $descA;
    while (my $nextDesc = $self->_getNextDesc($currentDesc)) {
        if ($nextDesc eq $descB) {
            # descB is after descA, so return 0
            return 0;
        }
        $currentDesc = $nextDesc;
    }    
    croak "cmpDescriptions: could not find $descB before or after $descA!";
}

sub _attributeAtIndex {
    my $self  = shift;
    my $index = shift;

    # Root is not included in index
    my $currentAttrDesc = $self->_root();
    
    my $j = 0;
    while (my $nextAttrDesc = $self->_getNextDesc($currentAttrDesc)) {
        return $nextAttrDesc if $j == $index;
        $currentAttrDesc = $nextAttrDesc;
        ++$j;
    }
    croak "Did not find an attribute at index $index";
}

sub _getPrevDesc {
    my $self     = shift;
    my $attrDesc = shift;

    return $attrDesc->_prevAttrDesc() eq $self->_root ? 0
        :  $attrDesc->_prevAttrDesc();
}

sub _getNextDesc {
    my $self     = shift;
    my $attrDesc = shift;
    
    return $attrDesc->_nextAttrDesc() eq $self->_head ? 0
        :  $attrDesc->_nextAttrDesc();
}

sub getDescriptionWithName {
    my $self = shift;
    my $attrName = shift;

    my $currentAttrDesc = $self->_root;

    while (my $nextAttrDesc = $self->_getNextDesc($currentAttrDesc)) {
        return $nextAttrDesc if $nextAttrDesc->name eq $attrName;
        $currentAttrDesc = $nextAttrDesc;
    }    
}
    
sub getIndexFromName {
    my $self     = shift;
    my $attrName = shift;

    my $currentAttrDesc = $self->_root;

    my $i = 0;
    while (my $nextAttrDesc = $self->_getNextDesc($currentAttrDesc)) {
        return $i if $nextAttrDesc->name eq $attrName;
        $currentAttrDesc = $nextAttrDesc;
        ++$i;
    }
}

__PACKAGE__->meta->make_immutable;

1;
