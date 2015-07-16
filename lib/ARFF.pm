package ARFF;
use Moose;
use Carp;
use WEKA;
use write2tmp;
use types;
use ARFF::Types;

use ARFF::ListOfAttributeDescriptions;
use ARFF::Attribute;
use ARFF::AttributeDescription;
use ARFF::Instance;
use ARFF::FileParser;

use Moose::Util::TypeConstraints;
use Moose::Meta::Attribute::Native::Trait::Array;

use List::Util qw(all none);

### Attributes #################################################################
################################################################################

has 'instances' => (
    traits  => ['Array'],
    isa     => 'ArrayRef[ARFF::Instance]',
    is      => 'rw',
    handles => {
        'allInstances' => 'elements',
        'getInstance'  => 'get',
        'findInstance' => 'first', 
    }
);

has 'attributeDescriptions' => (
    isa     => 'ARFF::ListOfAttributeDescriptions',
    is      => 'rw',
    coerce  => 1,
    handles => {
        getAttributeDescriptions           => 'getDescriptions',
        addAttributeDescription            => 'addDescription',
        removeAttributeDescription         => 'removeDescription',
        removeAttributeDescriptionWithName => 'removeDescriptionWithName',
        attributeDescriptionWithName       => 'getDescriptionWithName',
        attributeDescriptionIndexFromName  => 'getIndexFromName',
        countAttributeDescriptions         => 'count',
    },
);

has 'relation' => (
    is => 'rw',
    isa => 'Str'
);

### Methods ####################################################################
################################################################################

# TODO: mergeArff based on attributes rather than attribute lines
sub mergeArff {
    my $self = shift;
    my $arff2merge = shift;

    # Attribute lines of arff to be merged must match self attribute lines
    croak "arff attributes do not match!"
        if ! $self->_arffAttributesMatch($self, $arff2merge);

    push(@{$self->instances}, @{$arff2merge->instances});
}

sub _arffAttributesMatch {
    my ($arff1, $arff2) = @_;

    
    return 0 if @{$arff1->_attributeLines} ne @{$arff2->_attributeLines};
    
    for (my $i = 0 ; $i < @{$arff1->_attributeLines} ; ++$i) {
        my $arff1Attr = $arff1->_attributeLines->[$i];
        my $arff2Attr = $arff2->_attributeLines->[$i];

        return 0 if $arff1Attr ne $arff2Attr;
    }
    return 1;
}

sub addAttributeValuesToInstances {
    my $self              = shift;
    my $attrDescription   = shift;
    my $attrValues        = shift;
    
    if (ref $attrValues eq 'ARRAY') {
        # We expect values to be a list of values ordered to match the order
        # of instances
        for (my $i = 0 ; $i < @{$attrValues} ; ++$i) {
            my $value = $attrValues->[$i];
            my $attr  = ARFF::Attribute->new(description => $attrDescription,
                                             value       => $value);
            $self->getInstance($i)->addAttribute($attr);
        }
    }
    elsif (ref $attrValues eq 'HASH') {
        # We expect values to be in a hash where each key corresponds to a
        # an instance id
        while (my ($id, $value) = each %{$attrValues}) {
            my $attr  = ARFF::Attribute->new(description => $attrDescription,
                                       value       => $value);
            $self->findInstance(sub {$_->id eq $id} )->addAttribute($attr);
        }
    }
    else {
        croak "addAttributeValuesToInstances must be passed a ref to an array "
            . "or a hash of values";
    }
}

sub addAttributesWithDescriptionToAllInstances {
    my $self     = shift;
    my $attrDesc = shift;

    map {$_->addAttribute(_newAttributeWithDescription($attrDesc))}
        $self->allInstances();
}

sub _newAttributeWithDescription {
    my $attrDesc = shift;
    return ARFF::Attribute->new(description => $attrDesc, value => '?');
}

# TO removeAttributeWithName, we remove the attribute with that name
# from each instance and then remove the attributeDescription with
# that name from attributeDescriptions
sub removeAttributeWithName {
    my $self     = shift;
    my $attrName = shift;

    map {$_->removeAttributeWithName($attrName)} $self->allInstances();

    $self->removeAttributeDescriptionWithName($attrName);
}

sub arff2File {
    my $self = shift;

    return write2tmp->new(data => [$self->arff2String],
                          suffix => '.arff')->file_name();
}

sub arff2String {
    my $self = shift;

    my $string = "";

    $string .= $self->getRelationString() . "\n";
    
    my @attrDescs = $self->getAttributeDescriptions();
    
    # Add attribute lines
    $string .= $_ . "\n" foreach $self->getAttributeDescriptionStrings();
    
    # Add space and data header
    $string .= "\n\@data\n";

    # Add instances
    $string .= $_ . "\n" foreach $self->getInstanceStrings();
    
    return $string;
}

sub getRelationString {
    my $self = shift;

    my $relation = defined $self->relation() ? $self->relation
        : 'dataset';
    return '@relation ' . $relation;
}

sub getAttributeDescriptionStrings {
    my $self = shift;

    return map {$_->string()} $self->getAttributeDescriptions;
}

sub getInstanceStrings {
    my $self = shift;
    
    my @attrDescs = $self->getAttributeDescriptions();
    return map {$_->string(@attrDescs)} @{$self->instances()};
}

sub transAttributeWithNameToBinaryFromNominal {
    my $self      = shift;
    my $attrName  = shift;

    # Column numbers are indexed starting from 1, so increment accordingly
    my $attrColNum
        = $self->attributeDescriptionIndexFromName($attrName) + 1;

    # Run weka to get arff with nominal attribute split into binary attrs
    my $binArffStr = WEKA->new->nominalToBinary($self->arff2File(),
                                                $attrColNum);
    
    # parse arff and find names of new binary attributes
    my $binArff
        = ARFF::FileParser->new(allLines => [split("\n", $binArffStr)])->parse();

    my @binaryAttrDescs
        = grep {$_->name() =~ /$attrName/}
            $binArff->getAttributeDescriptions;
    croak "Could not find binarized attribute descriptions for nominal attribute $attrName!"
        if ! @binaryAttrDescs;
    
    # For each new binary attribute, get an array of values.
    # Add each value its corresponding pre-existing instance.
    for (my $i = 0 ; $i < @binaryAttrDescs ; ++$i) {
        my $binAttrDesc = $binaryAttrDescs[$i];
        
        # Add description to self and then array of attribute values to
        # corresponding instances
        my $addPosition
            = $self->attributeDescriptionIndexFromName($attrName) + 1 + $i;
        $self->addAttributeDescription($binAttrDesc, $addPosition);
        
        my @attrValues
            = map {$_->getValueForAttributeWithName($binAttrDesc->name)}
                @{$binArff->instances};
        
        $self->addAttributeValuesToInstances($binAttrDesc, \@attrValues);
    }
    # Remove old nominal attribute
    $self->removeAttributeWithName($attrName);
}

sub standardize {
    my $self          = shift;
    my $referenceArff = shift;

    # User could pass either another arff object our an arff file
    my $refFile
        = ref ($referenceArff) eq ref ($self) ? $referenceArff->arff2File()
            : $referenceArff;
    
    my ($stdizedRefArffFile, $stdizedSelfArffFile)
        = WEKA->new->standardizeArffFiles($refFile, $self->arff2File());
    
    my $stdizedArff
        = ARFF::FileParser->new(file => $stdizedSelfArffFile)->parse();

    for (my $i = 0 ; $i < @{$stdizedArff->instances} ; ++$i) {
        my $stdizedInstance = $stdizedArff->instances->[$i];
        my $selfInstance    = $self->instances->[$i];

        $selfInstance->copyValuesFrom($stdizedInstance);
    }
}

sub doAllInstancesHaveAttributeWithName {
    my $self     = shift;
    my $attrName = shift;
    
    return all  {$_->hasAttributeWithName($attrName)} $self->allInstances();
}

sub doAllInstancesLackAttributeWithName {
    my $self     = shift;
    my $attrName = shift;
    return none {$_->hasAttributeWithName($attrName)} $self->allInstances();
}

sub getValuesForAttributeWithName {
    my $self     = shift;
    my $attrName = shift;
    return map {$_->getValueForAttributeWithName($attrName)} $self->allInstances();
}

__PACKAGE__->meta->make_immutable;

################################ END OF PACKAGE ################################
################################################################################

1;
