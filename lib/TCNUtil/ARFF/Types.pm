package ARFF::Types;
use Moose::Util::TypeConstraints;
use TCNUtil::ARFF::AttributeDescription;
use TCNUtil::ARFF::ListOfAttributeDescriptions;

coerce 'ARFF::ListOfAttributeDescriptions',
    from 'ArrayRef[ARFF::AttributeDescription]',
    via {ARFF::ListOfAttributeDescriptions->new(@{$_})};

1;
