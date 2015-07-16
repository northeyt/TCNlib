package ARFF::FileParser;
use Moose;
use Carp;

use types;
use ARFF;

has 'file' => (
    is  => 'rw',
    isa => 'FileReadable',
);

has 'headerLines' => (
    isa     => 'ArrayRef',
    is      => 'rw',
    lazy    => 1,
    builder => '_parseHeaderLines'
);

has 'attributeLines' => (
    isa     => 'ArrayRef',
    is      => 'rw',
    lazy    => 1,
    builder => '_parseAttributeLines'
);

has 'relationLine' => (
    isa => 'Str',
    is => 'rw',
    lazy => 1,
    builder => '_parseRelationLine',
);

has 'instanceLines' => (
    isa => 'ArrayRef',
    is => 'rw',
    lazy => 1,
    builder => '_parseInstanceLines'
);

has 'allLines' => (
    isa => 'ArrayRef',
    is => 'rw',
    lazy => 1,
    builder => '_getLinesFromFile'
);


has 'idAttributeIndex' => (
    is => 'rw',
    predicate => 'has_idAttributeIndex',
);

has '_realidAttributeIndex' => (
    is => 'rw',
    isa => 'Int',
    lazy => 1,
    builder => '_buildRealIDAttributeIndex'
);

has 'classAttributeIndex' => (
    is => 'rw',
    predicate => 'has_classAttributeIndex'
);

has '_realClassAttributeIndex' => (
    is => 'rw',
    isa => 'Int',
    lazy => 1,
    builder => '_buildRealClassAttributeIndex'
);

sub _parseAttributeLines {
    my $self = shift;

    my @attributeLines = ();
    
    foreach my $line (@{$self->headerLines}) {
        if ($line =~ /^\@attribute/) {
            $line =~ s/\@attribute //;
            push(@attributeLines, $line);
        }
    }
    return \@attributeLines;
}

sub _parseRelationLine {
    my $self = shift;

    foreach my $line (@{$self->headerLines}) {
        if ($line =~ /^\@relation/) {
            $line =~ s/\@relation //;
            return $line;
        }
    }    
}

sub _parseHeaderLines {
    my $self = shift;

    my @headerLines = ();
    foreach my $line (@{$self->allLines}) {
        if ($line =~ /^\@data/) {
            last;
        }
        next if $line =~ /^\n$/;
        push(@headerLines, $line);
    }
    return \@headerLines;
}

sub _parseInstanceLines {
    my $self = shift;

    my @instanceLines = ();

    my $reachedHeader = 0;
    
    foreach my $line (@{$self->allLines}) {
        if ($line =~ /^\@data/) {
            $reachedHeader = 1;
            next;
        }
        elsif ($reachedHeader && $line) {
            push(@instanceLines, $line);
        }
    }
    return \@instanceLines;
}

sub _getLinesFromFile {
    my $self = shift;

    my $inputFile = $self->file;
    open(my $IN, "<", $inputFile) or confess "Cannot open file $inputFile, $!";

    my @contentLines = ();
    
    foreach my $line (<$IN>) {
        next if $line =~ /^\n$/;
        chomp $line;
        push(@contentLines, $line);
    }
    
    return \@contentLines;
}

sub attributeDescriptionsFromLines {
    my $self = shift;

    my @attributeDescriptions = ();
    
    for (my $i = 0 ; $i < @{$self->attributeLines} ; ++$i) {
        my $attrDesc = ARFF::AttributeDescription->new($self->attributeLines->[$i]);
        $attrDesc->is_id(1)
            if $self->has_idAttributeIndex
                && $self->_realidAttributeIndex == $i;

        $attrDesc->is_class(1)
            if $self->has_classAttributeIndex
                && $self->_realClassAttributeIndex == $i;
        
        push(@attributeDescriptions, $attrDesc);
    }
    return @attributeDescriptions;
}

sub instancesFromLines {
    my $self = shift;
    my @attributeDescriptions = @_;
    
    my @instances = ();
    foreach my $line (@{$self->instanceLines}) {
        my @values = split(/,/, $line);

        croak "Number of values does not equal number of attributes! "
            . "Num attributes = " . $self->countAttributeDescriptions()
                . "\nNum values = " . @values
                    if @values != @attributeDescriptions;
        
        my @instanceAttributes = ();
        
        for (my $i = 0 ; $i < @values ; ++$i) {
            my $attr = ARFF::Attribute->new(description => $attributeDescriptions[$i],
                                            value       => $values[$i]);
            push(@instanceAttributes, $attr);
        }
        
        my $inst = ARFF::Instance->new();
        $inst->addAttributes(@instanceAttributes);
        push(@instances, $inst);
    }        
    return @instances;
}

sub _buildRealClassAttributeIndex {
    my $self = shift;

    my $classAttributeIndex = $self->classAttributeIndex;
    return $self->_getRealIndex($classAttributeIndex);
}

sub _buildRealIDAttributeIndex {
    my $self = shift;

    my $idAttributeIndex = $self->idAttributeIndex;
    return $self->_getRealIndex($idAttributeIndex);
}

# This translates 'first' and 'last' indices to their numerical values
sub _getRealIndex {
    my $self = shift;
    my $index = shift;
    
    my $realIndex = $index eq 'first' ? 0
        : $index eq 'last' ? scalar @{$self->attributeLines} - 1
            : $index;

    return $realIndex;
}

sub parse {
    my $self = shift;

    my $relation = $self->relationLine();
    my @attributeDescriptions = $self->attributeDescriptionsFromLines();
    my @instances = $self->instancesFromLines(@attributeDescriptions);
    
    my $arff = ARFF->new(attributeDescriptions => \@attributeDescriptions,
                         instances => \@instances,
                         relation  => $relation);
    return $arff;
}

1;
