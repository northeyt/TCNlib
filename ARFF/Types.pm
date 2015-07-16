package ARFF::Types;
use Moose::Util::TypeConstraints;
use ARFF::AttributeDescription;
use ARFF::ListOfAttributeDescriptions;

coerce 'ARFF::ListOfAttributeDescriptions',
    from 'ArrayRef[ARFF::AttributeDescription]',
    via {ARFF::ListOfAttributeDescriptions->new(@{$_})};

1;
