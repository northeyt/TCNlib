package WEKA;

=head1 NAME


=cut

=head1 SYNOPSIS

=cut


=head1 DESCRIPTION

=cut

use Moose;
use Moose::Util::TypeConstraints;
use TCNPerlVars;
use TCNUtil::types;
use TCNUtil::GLOBAL qw(rm_trail);
use Carp;
use TCNUtil::confusion_table;
use IO::CaptureOutput qw(qxx);

=head1 Attributes

=over 12

=cut

### Attributes #################################################################
################################################################################

=item C<wekaLib>

Path to weka library.

=cut

has 'wekaLib' => (
    is => 'rw',
    isa => 'FileReadable',
    default => $TCNPerlVars::wekaLib,
    required => 1,
    lazy => 1,
);

=item C<java>

Path to java.

=cut

has 'java' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::javaForWeka,
    required => 1,
    lazy => 1,
);

=item C<memory>

Memory limit for WEKA. Written with m or g notation, e.g.

100m
20g

=cut

has 'memory' => (
    is => 'rw',
    isa => 'Str',
    default => '20g',
    required => 1,
    lazy => 1,
);


=item C<remoteMachine>

Name of remote machine for WEKA to be run on.
If not set, WEKA will be locally.

=cut

has 'remoteMachine' => (
    is => 'rw',
    isa => 'Str', # Could write a type for this?
    predicate => 'has_remoteMachine',
);

foreach my $label (qw(posLabel negLabel undefLabel)) {
    has $label => (
        is => 'rw',
        isa => 'Str',
        predicate => "has_$label",
    );
}

has 'translateUndefLabelTo' => (
    is        => 'rw',
    isa       => enum([qw[ posLabel negLabel ]]),
    predicate => 'has_translateUndefLabelTo',
    clearer   => 'clear_translateUndefLabelTo',
);

has '_labelMap' => (
    is => 'rw',
    isa => 'HashRef',
    lazy => 1,
    builder => '_buildLabelMap',
);

### METHODS ####################################################################
################################################################################

# This string will be used to prefix all WEKA system commands
sub _cmdBasePrefix {
    my $self = shift;
    return join(" ", ($self->java, "-Xmx". $self->memory, "-cp", $self->wekaLib));
}

sub nominalToBinary {
    my $self           = shift;
    my $inputArffFile  = shift;
    my $nomAttrColNum  = shift;

    my $wekaClass      = "weka.filters.unsupervised.attribute.NominalToBinary";
    my $cmd            = join(" ", ($self->_cmdBasePrefix(), $wekaClass,
                                    "-R $nomAttrColNum", "-i $inputArffFile"));

    # Run WEKA
    my($outArffStr, $stderr, $success) = qxx($cmd);
    if ($stderr) {
       croak "WEKA STDERR: $stderr\nCommand Run: $cmd";
    }
    
    return $outArffStr;
}

sub standardizeArffFiles {
    my $self = shift;

    # Primary arff's statistics are used to standardize itself AND secondary.
    my $primaryArffFile   = shift;
    my $secondaryArffFile = shift;
    
    my $cmd;
    $cmd .= $self->_cmdBasePrefix();
    $cmd .= " weka.filters.unsupervised.attribute.Standardize";

    my $secondaryStdizedFile;
    if (defined $secondaryArffFile) {
        # Create file for standardised secondary arff
        $secondaryStdizedFile
            = write2tmp->new(suffix => ".arff")->file_name();
        
        $cmd .= " -b -r $secondaryArffFile -s $secondaryStdizedFile";
    }
    
    # Create file for standardised primary arff
    my $primaryStdizedFile
        = write2tmp->new(suffix => ".arff")->file_name();
    
    $cmd .= " -i $primaryArffFile -o $primaryStdizedFile";
    
    # Run WEKA
    my($stdout, $stderr, $success) = qxx($cmd);
    if ($stderr) {
        croak "WEKA STDERR: $stderr\nCommand Run: $cmd";
    }
    
    return defined $secondaryStdizedFile ? ($primaryStdizedFile, $secondaryStdizedFile)
        : $primaryStdizedFile;
}

sub parseTableFromOutput {
    my $self    = shift;
    my $output  = shift;
    my $outForm = shift;

    $outForm = 'CSV' if ! $outForm;
    
    # If output in array ref form, use this - otherwise,
    # split output into lines
    my $outputAref = ref $output eq 'ARRAY' ? $output :
        [split /(?<=\n)/, $output]; # Split output into lines (keeping newline)
    
    # Get confusion table from output
    my $table = $self->_tableFromLines($outputAref, $outForm);
}

sub _tableFromLines {
    my $self     = shift;
    my $lineAref = shift;
    my $outForm  = shift;
    
    my @instances
        = $outForm eq 'CSV' ? $self->_getInstancesFromCSVLines($lineAref)
        : $self->_getInstancesFromDefaultOutputLines($lineAref);
        
    my $table     = confusion_table->new(item_class => 'instance');
    map {$table->add_datum($_)} @instances;
    return $table;
}

sub  _getInstancesFromDefaultOutputLines {
    my $self         = shift;
    my $lineAref     = shift;
    my @fieldTitles  = qw(inst actual predicted score id);

    my @instances = ();
    foreach my $line (@{$lineAref}) {
        chomp $line;
        next if ! $line || $line =~ /^\n$/;
        
        my $instance = $self->_instanceFromLine($line, \@fieldTitles, 'DEF');
        push(@instances, $instance) if $instance;
    }
    return @instances;
}

sub _getInstancesFromCSVLines {
    my $self         = shift;
    my $lineAref     = shift;
    
    my @instances     = ();
    my @fieldTitles   = ();
    my $reachedHeader = 0;
    
    foreach my $line (@{$lineAref}) {
        chomp $line;
        next if ! $line;
        
        if (! $reachedHeader || $line =~ /^\n$/) {
            if ($line =~ /^inst#/) {
                $reachedHeader = 1;
                @fieldTitles = split(",", $line);
            }
            next;
        }
                
        my $instance = $self->_instanceFromLine($line, \@fieldTitles, 'CSV');
        push(@instances, $instance) if $instance; 
    }
    return @instances;
}

sub _instanceFromLine {
    my $self           = shift;
    my $line           = shift;
    my $fieldTitleAref = shift;
    my $lineForm       = shift;
    
    my ($instNum, $valueLabel, $predLabel, $predValue, @remainingFields)
        = $lineForm eq 'CSV' ? $self->_parseCSVLine($line)
        : $self->_parseDefaulOutputLine($line);

    # Remove numeric label identifier from value and prediction
    # label, e.g. 1:I => I
    ($valueLabel, $predLabel)
        = map {[split(":", $_)]->[1]} ($valueLabel, $predLabel);
    
    if ($valueLabel eq $self->undefLabel) {
        if ($self->has_translateUndefLabelTo) {
            my $labelType = $self->translateUndefLabelTo;
            $valueLabel = $self->$labelType;
        }
        else {
            # If this instance's value is unlabelled then we cannot add this
            # to a table, so do not return an instance
            return 0;
        }
    }
    
    my $labelTest
        = eval{$self->_checkLabels([$valueLabel, $predLabel])};
    croak $@ if ! $labelTest;

    my %valueForField = ();
    @valueForField{@{$fieldTitleAref}} = ($instNum, $valueLabel, $predLabel,
                                          $predValue, @remainingFields);
    
    my $obj = bless \%valueForField, 'instance';
    my $instance = datum->new(object => $obj,
                              value => $self->_labelMap->{$valueLabel},
                              prediction => $self->_labelMap->{$predLabel});
    
    return $instance;
}

sub _parseDefaulOutputLine {
    my $self = shift;
    my $line = rm_trail(shift);

    # example line:      1        2:S        2:S   +   0.9 (2qqk:A:415)
    # Remove + that indicates error, but is not present otherwise
    $line    =~ s/\+//;
    # Remove brackets around patch ID if it is present
    $line    =~ s/[()]//;

    my ($instNum, $valueLabel, $predLabel, $predValue, @remainingFields)
        = split(/\s+/, $line);
    
    return ($instNum, $valueLabel, $predLabel, $predValue, @remainingFields);
}

sub _parseCSVLine {
    my $self = shift;
    my $line = shift;
    
    # example line: 2857,1:I,2:S,+,0.887
    my($instNum, $valueLabel, $predLabel, $err, $predValue, @remainingFields)
        = split(",", $line);
    return ($instNum, $valueLabel, $predLabel, $predValue, @remainingFields);
}
   
sub _checkLabels {
    my $self      = shift;
    my $labelAref = shift;
    
    # Check that passed labels are defined and match the labels specified in
    # the map
    foreach my $label (@{$labelAref}) {
        if (! defined $label) {
            croak "Label is not defined!\n";
        }
        elsif (! exists $self->_labelMap->{$label}) {
            croak "label '$label' does not match user-input labels!\n";
        }
    }
    return 1;
}

sub _buildLabelMap {
    my $self     = shift;
    croak "posLabel must be set to create a label map!" if ! $self->has_posLabel;
    croak "negLabel must be set to create a label map!" if ! $self->has_negLabel;
    return {$self->posLabel => 1, $self->negLabel => 0};
}

__PACKAGE__->meta->make_immutable;

################################ END OF PACKAGE ################################
################################################################################

package WEKA::randomForest;

=head1 NAME


=cut

=head1 SYNOPSIS

=cut


=head1 DESCRIPTION

=cut

use Moose;
use Carp;

use TCNPerlVars;
use TCNUtil::types;
use Moose::Util::TypeConstraints;

use IO::CaptureOutput qw(qxx);
use File::Spec;
use File::Basename;
use TCNUtil::ARFF;

extends 'WEKA';

=head1 Attributes

=over 12

=cut

=item C<filterClass>

Name of name of filter class to be used.

=cut

has 'filterClass' => (
    is => 'rw',
    isa => 'Str',
    default => "weka.classifiers.meta.FilteredClassifier",
    lazy => 1,
);

=item C<isFiltered>

BOOL. If TRUE, model will created/run as a filtered classifier.

=cut

has 'isFilter' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

=item C<class>

Name of RF class to be used

=cut

has 'class' => (
    is => 'rw',
    isa => 'Str',
    default => "weka.classifiers.trees.RandomForest",
    lazy => 1,
    required => 1,
);

=item C<model>

Path of model to be run. If train is run, the output model's path is set to this
attribute.

=cut
    
has 'model' => (
    is => 'rw',
    predicate => 'has_model'
);

=item C<trainArff>

Path to training set .arff file.

=cut

=item C<testArff>

Path to test set .arff file.

=cut

coerce 'ARFF',
    from 'FileReadable',
    via { ARFF::FileParser->new(file => $_)->parse() };

foreach my $set (qw(train test)) {

    has $set . 'Arff' => (
        is => 'rw',
        isa => 'ARFF',
        predicate => "has_" . $set . "Arff",
        coerce => 1,
    );
    
}

=item C<classAttrCol>

Column number of class attribute in train/test data. Keyword strings can also be
used, i.e. 'first' or  'last'. Default = 'last'

=cut

has 'classAttrCol' => (
    isa     => 'Str',
    is      => 'rw',
    lazy    => 1,
    builder => '_buildClassAttrCol'    
);

sub _buildClassAttrCol {
    my $self = shift;

    my $arffClassAttrIndex = $self->_getClassAttrIndexFromArffs();

    # column = index + 1
    return defined $arffClassAttrIndex ? $arffClassAttrIndex + 1
        : 'last';
}

sub _getClassAttrIndexFromArffs {
    my $self = shift;
    
    my $testArffClassAttrIndex  = eval {$self->testArff->getClassAttrIndex};
    my $trainArffClassAttrIndex = eval {$self->trainArff->getClassAttrIndex};
    
    if (defined $testArffClassAttrIndex && defined $trainArffClassAttrIndex) {
        croak "Test and training arff set class attributes indices do not match!"
            if $testArffClassAttrIndex != $trainArffClassAttrIndex; 
    }
    
    return defined $testArffClassAttrIndex  ? $testArffClassAttrIndex
        :  defined $trainArffClassAttrIndex ? $trainArffClassAttrIndex
        :  undef;
}

=item C<numTrees>

Number of trees.

=cut

has 'numTrees' => (
    isa => 'Int',
    is => 'rw',
    required => 1,
    lazy => 1,
    default => 100,
);

=item C<numFeatures>

Number of features per tree.

=cut

has 'numFeatures' => (
    isa => 'Int',
    is => 'rw',
    required => 1,
    lazy => 1,
    default => 3,
);

=item C<removeAttribute>

Index of any attributes to remove during training or testing. Can be a single
index or formatted in a way which WEKA will recognise

e.g. "1"
     "1,9-11"
     "1-10"

=cut

has 'removeAttribute' => (
    is => 'rw',
    isa => 'Str',
    predicate => 'has_removeAttribute',
);

=item C<removeAttrClass>

Path to class used for removing attributes. Required if removeAttribute is set

=cut

has 'removeAttrClass' => (
    is => 'rw',
    isa => 'Str',
    predicate => 'has_removeAttrClass',
    default => "weka.filters.unsupervised.attribute.Remove",
);

### METHODS ####################################################################
################################################################################

sub test {
    my $self = shift;

    my $cmd = $self->buildTestCmd();
    
    my($stdout, $stderr, $success) = qxx($cmd);

    if ($stderr) {
       croak "WEKA STDERR: $stderr\nCommand Run: $cmd";
    }
    return $stdout;
}

sub train {
    my $self = shift;

    my $cmd = $self->buildTrainCmd();
    
    my($stdout, $stderr, $success) = qxx($cmd);

    if ($stderr) {
       croak "WEKA STDERR: $stderr\nCommand Run: $cmd";
    }
    return $stdout;
}

sub buildTestCmd {
    my $self = shift;

    my $metaClass  = $self->filterClass;    
    my $model = $self->model();
    croak "A saved model must be supplied for test. "
        . "Set 'model' to a file path." if ! $model;

    my $testArff = $self->testArff->arff2File();
    croak "A testing .arff must be supplied to test! Set 'testArff' to a "
        . "test .arff file" if ! $testArff;

    my $csvOutputStr = "-classifications \"weka.classifiers.evaluation.output.prediction.CSV -p first \"";

    my $cmdBase = $self->_cmdBasePrefix();
    my $cmd = "$cmdBase $metaClass -T $testArff -l $model $csvOutputStr";

    if ($self->has_remoteMachine) {
        my $machine = $self->remoteMachine;
        # Any embedded quotes within command must be escaped
        $cmd =~ s{"}{\\"}xmsg;
        $cmd = "ssh $machine $cmd";
    }
    return $cmd;
}

sub buildTrainCmd {
    my $self = shift;

    croak "A training .arff must be supplied to train! Set 'trainArff' to a "
        . "training .arff file" if ! $self->has_trainArff;
    my $trainArff = $self->trainArff->arff2File();
    
    croak "A file to save a model to must be supplied for training. "
        . "Set 'model' to a file path." if ! $self->has_model;
    my $model = $self->model();
    
    my $testArff    = $self->testArff->arff2File() if $self->has_testArff;
    my $numFeatures = $self->numFeatures;
    my $numTrees    = $self->numTrees;
    my $rfClass     = $self->class;

    # Create cmd
    my $cmdBase = $self->_cmdBasePrefix();
    my $dataIn  = "-d $model -t $trainArff";
    $dataIn .= defined $testArff ? " -T $testArff" : "";
    
    my $classAttrStr = "-c " . $self->classAttrCol;
    my $filterWOpts = "";

    my $rfOptsStr  = "-I $numTrees -K $numFeatures";

    my $cmd = "";
    
    if ($self->has_removeAttribute){
        # Build up string to include meta class and rAttr class
        my $rAttrClass = $self->removeAttrClass;
        my $rAttrStr   = $self->removeAttribute;
        my $metaClass  = $self->filterClass;
        
        # Quotes required for cmd to work
        my $filterWOpts =  "-F \"$rAttrClass -R $rAttrStr\"";
        my $rfClassStr = "-W $rfClass";
        my $csvOutputStr = "-classifications \"weka.classifiers.evaluation.output.prediction.CSV -p first \"";
        
        $cmd = "$cmdBase $metaClass $dataIn $filterWOpts"
            . " $classAttrStr $csvOutputStr $rfClassStr -- $rfOptsStr"; 
    }
    else {
        croak "train is not yet implemented to run without an attribute filter!";
    }
    
    if ($self->has_remoteMachine) {
        my $machine = $self->remoteMachine;
        # Any embedded quotes within command must be escaped
        $cmd =~ s{"}{\\"}xmsg;
        $cmd = "ssh $machine $cmd";
    }
    return $cmd;
}

### METHOD MODIFIERS ###########################################################
################################################################################

around 'model' => sub {
    
    my $orig = shift;
    my $self = shift;
  
    return $self->$orig()
        unless @_;
    
    my $file = shift;

    # User can specify a pre-existing model file or a path to where they would
    # like to write a model.
    if (-e $file) {
        # If file already exists, get full path
        $file = File::Spec->rel2abs($file);
    }
    else {
        # Check that dir path of currently non-existant file is valid,
        # before getting full path
        my ($name, $path, $suffix) = fileparse($file);
        
        if (-d $path) {
            my $fullPath2Dir = File::Spec->rel2abs($path);
            $file = "$fullPath2Dir/$name";
        }
        else {
            croak "Could not parse directory from filename $file!";
        }
    }
    return $self->$orig($file);
    
};

__PACKAGE__->meta->make_immutable;

################################ END OF PACKAGE ################################
################################################################################

package WEKA::randomForest::partitionCV;

=head1 NAME


=cut

=head1 SYNOPSIS

=cut


=head1 DESCRIPTION

=cut

use Moose;
use Carp;

use TCNPerlVars;
use TCNUtil::types;

extends 'WEKA::randomForest';

### ATTRIBUTES #################################################################
################################################################################

=head1 Attributes

=over 12

=cut

=item C<partitionAref>

Reference to array where element has form :
[trainingArff, testArff, outputModel, outputFile]

=cut

has 'partitionAref' => (
    isa => 'ArrayRef',
    is => 'rw',
    required => 1,
);


### METHODS ####################################################################
################################################################################

=head1 Methods

=over 12

=cut

=item  C<run>

=cut

sub run {
    my $self = shift;

    my @tables = ();
    
    foreach my $partitionAref (@{$self->partitionAref}) {
        
        my ($trainArff, $testArff, $model, $outFile) = @{$partitionAref};
        
        $self->trainArff($trainArff);
        $self->testArff($testArff);
        $self->model($model);

        my $output = $self->train();
        my $table = $self->parseTableFromOutput($output);

        if ($outFile) {
            open(my $OUT, ">", $outFile) or die "Cannot open file $outFile, $!";
            print {$OUT} $output;
        }        
        push(@tables, $table);
    }
    
    return @tables;
}


### !!! NOTE !!! ###############################################################
################################################################################
# The two following subroutines are used to parse WEKA output when the -i option
# is given to WEKA to output summary performance metrics. These subroutines work
# but are currently not implemented, as these classes currently parse CSV WEKA
# output.

sub parsePerformanceMetrics {

    my $self = shift;
    my $outputStr = shift;

    my $sectionMainHeader = $self->has_testArff ? "=== Error on test data ==="
        : "=== Stratified cross-validation ===";

    my $sectionSubHeader = "=== Detailed Accuracy By Class ===";

    my $reachedMainHeader = 0;
    my $reachedSubHeader = 0;

    my @metricTypes = ();
    my @values = ();

    foreach my $line (split(/\n/, $outputStr)) {

        if ($line =~ /^Weighted Avg/ && $reachedMainHeader && $reachedSubHeader) {
            # Accuracy is the weighted average recall, so can be found by
            # parsing the weighted average

            croak "You're parsing the Weighted Avg. line before you've parsed "
                . "the main metrics line, something has gone wrong"
                    if ! @values;
            
            my @avgValues = split(/\s+/, $line);
            my $metric = "Accuracy";
            my $value  = $avgValues[4];
            push(@metricTypes, $metric);
            push(@values, $value);
        }
        
        last if @values;
        
        if ($line =~ /$sectionMainHeader/) {
            $reachedMainHeader = 1;
            next;
        }
        elsif ($reachedMainHeader && $line =~ /$sectionSubHeader/) {
            $reachedSubHeader = 1;
            next;
        }
        
        next if (! ($reachedMainHeader && $reachedSubHeader)) || $line =~ /^\s+$/;

        # Remove trailing white-space
        $line =~ s/^\s+|\s+$//g;

        # Hyphenate metrics names
        $line =~ s/([A-Z])\s([A-Z])/$1-$2/g;
                
        if (! @metricTypes) {
            @metricTypes = split(/\s+/, $line);
            next;
        }
        
        # Avoid performance metrics for 0 class - we want metrics for 1 class
        @values = split(/\s+/, $line) if $line !~ /0$/;
        
    }
    
    my %hash = ();
    @hash{@metricTypes} = @values;
    
    if ($hash{"FP-Rate"}) {
        # Specificity is not supplied by WEKA, but can be calculated from
        # FP-Rate 
        my $metric = "Specificity";
        my $value  = 1 - $hash{"FP-Rate"};
        $hash{$metric} = $value;
    }

    # Remove Class 'metric', which is not actually a metric
    delete $hash{Class};
    
    return %hash;
}

sub averageMetrics {
    my $self = shift;
    my $metricsAref = shift;

    my $n = scalar @{$metricsAref};
    
    my @keys = sort {$a cmp $b} keys %{$metricsAref->[0]};
    
    my %totals = map {$_ => 0} @keys;
    
    # Sum each metric
    foreach my $metricHref (@{$metricsAref}) {
        foreach my $key (keys %{$metricHref}) {
            # Metric values will be '?' when they can't be calculated,
            # e.g. MCC when there are no positives. Therefore replace '?'
            # with 0.
            $totals{$key} += $metricHref->{$key} ne '?' ? $metricHref->{$key}
                : 0;
        }
    }
    
    # Get averages
    my %averages = map {$_ => $totals{$_} / $n} keys %totals;

    return %averages;
}


__PACKAGE__->meta->make_immutable;

1;
