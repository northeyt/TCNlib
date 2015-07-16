package confusion_table;

use Moose;
use Carp;

has 'item_class' => (
    is => 'rw',
    isa => 'Str',
    required => 1,
);

has 'data' => (
    traits => ['Array'],
    is  => 'ro',
    isa => 'ArrayRef[datum]',
    default => sub { [] },
    handles => {
        list_data  => 'elements',
        _add_datum => 'push',
        },
);

# Methods
sub true_pos {
    my $self = shift;
    return scalar grep {$_->value == 1 && $_->prediction == 1} $self->list_data();
}

sub false_neg {
    my $self = shift;
    return scalar grep {$_->value == 1 && $_->prediction == 0} $self->list_data();
}

sub false_pos {
    my $self = shift;
    return scalar grep {$_->value == 0 && $_->prediction == 1} $self->list_data();
}

sub true_neg {
    my $self = shift;
    return scalar grep {$_->value == 0 && $_->prediction == 0} $self->list_data();
}

sub predicted {
    my $self = shift;
    return $self->true_pos + $self->false_pos;
}

sub actual {
    my $self = shift;
    return $self->true_pos + $self->false_neg;
}

sub total {
    my $self = shift;
    return scalar $self->list_data();
}

sub sensitivity {
    my $self = shift;

    return ( $self->true_pos / $self->actual );
}

sub specificity {
    my $self = shift;

    my $denom = $self->total - $self->actual;

    # Avoid dividing by zero
    return 0  if ! $denom;
    
    return ( $self->true_neg / $denom );
    
}

sub PPV {
    my $self = shift;

    my $denom = $self->true_pos + $self->false_pos;

    return 0 if ! $denom; # Avoid dividing by zero
    
    return ( $self->true_pos / $denom );
}

# False Positive Rate
sub FPR {
    my $self = shift;

    my $denom = $self->false_pos + $self->true_neg;

    return ($self->false_pos / $denom);
}

# False discovery rate
sub FDR {
    my $self = shift;

    my $denom = $self->false_pos + $self->true_pos;

    return ($self->false_pos / $denom);
}

sub accuracy {
    my $self = shift;

    my $total = $self->total();
    my $tp = $self->true_pos();
    my $tn = $self->true_neg();

    return (($tp + $tn) / $total);
}

sub MCC {
    my $self = shift;
    my $tp = $self->true_pos;
    my $tn = $self->true_neg;
    my $fp = $self->false_pos;
    my $fn = $self->false_neg;

    my $denom = sqrt(  ($tp + $fp)
                     * ($tp + $fn)
                     * ($tn + $fp)
                     * ($tn + $fn) );

    # if denom is 0, set to 1 (shown to be correct in limit)
    $denom = 1 if ! $denom;
    
    my $numer = ($tp * $tn) - ($fp * $fn);

    my $mcc = $numer / $denom;

    return $mcc;
}

sub print_table {
    my $self  = shift;
    print "\t"      . "Positive"     . "\t"   . "Negative"     . "\n"
        . "True\t"  . $self->true_pos() . "\t" . $self->true_neg() . "\n"
            . "False\t" . $self->false_pos() . "\t" . $self->false_neg() . "\n";
}

sub add_datum {
    my $self  = shift;
    my $datum = shift;
    
    croak "add_datum: passed arg is not a datum"
        if ref $datum ne 'datum';

    my $datum_obj_ref    = ref $datum->object();
    my $table_item_class = $self->item_class();
    
    croak "data object does not match table item class.\n"
        . "datum ref: $datum_obj_ref.\ntable item class: $table_item_class\n"
            if ref $datum->object ne $self->item_class;

    $self->_add_datum($datum);
}

sub metrics_array {
    my $self = shift;

    return qw(sensitivity specificity accuracy FPR PPV MCC FDR);
}


sub print_all {
    my $self = shift;
    
    $self->print_table();
    print "\n";
    
    foreach my $metric ($self->metrics_array()) {
        my $value;
        my $ret = eval {$value = $self->$metric(); 1};
        
        print "$metric: ", $ret ? $value : "???", "\n"; 
    }
    return 1;
}
    
sub hash_all {
    my $self = shift;

    my @methods = $self->metrics_array();

    my %hash = ();
    
    foreach (@methods) {
        $hash{$_} = $self->$_;
    }

    return %hash;
}

### Functions

# AUC ###
# This function calculates the Area under Receiver Operating Characteristic
# Curve (AUC) from an array of confusion tables
sub AUC {
    my @confusionTables = @_;

    # x = 1 - specificity, y = sensitivity
    my $xyPairs
        = [map {[1 - $_->specificity(), $_->sensitivity]} @confusionTables];

    # Sort pairs by x values in ascending order
    $xyPairs = [sort {$a->[0] <=> $b->[0]} @{$xyPairs}];

    # Prepend and append 0,0 and 1,1 values, if these pairs don't exist already
    if (! arePairsIdentical($xyPairs->[0], [0,0])){
        $xyPairs = [[0,0], @{$xyPairs}];
    }

    if (! arePairsIdentical($xyPairs->[-1], [1,1])) {
        $xyPairs = [@{$xyPairs}, [1,1]];
    }

    my $AUC = 0;
    
    # Use trapezoid forumla to calculate area under curve
    for (my $i = 1 ; $i < @{$xyPairs} ; ++$i) {
        my($x0, $y0) = @{$xyPairs->[$i - 1]};
        my($x1, $y1) = @{$xyPairs->[$i]};

        $AUC += (($y0 + $y1) / 2) * ($x1 - $x0);
    }

    return $AUC;
}

sub arePairsIdentical {
    my $pairA = shift;
    my $pairB = shift;

    if ($pairA->[0] == $pairB->[0] && $pairA->[1] == $pairB->[1]){
        return 1;
    }
    else {
        return 0;
    }
}

sub mergeTables {
    my @tables = @_;

    # Get item class from first table
    my $itemClass = $tables[0]->item_class();
    
    my $mergedTable = confusion_table->new(item_class => $itemClass);

    foreach my $table (@tables) {
        foreach my $datum (@{$table->data()}) {
            $mergedTable->add_datum($datum);
        }
    }
    return $mergedTable;
}

__PACKAGE__->meta->make_immutable;

################################ END OF PACKAGE ################################
################################################################################

package datum;
    
use Moose;
use Carp;

has 'object' => (
    is => 'rw',
    isa => 'Object',
    required => 1,
);

has 'prediction' => (
    is => 'rw',
    isa => 'Bool',
    required => 1,
);

has 'value' => (
    is => 'rw',
    isa => 'Bool',
    required => 1,
);


1;
__END__

=head1 NAME

confusion_table - Perl extension for blah blah blah

=head1 SYNOPSIS

   use confusion_table;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for confusion_table, 

Blah blah blah.

=head2 EXPORT

None by default.

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Tom Northey, E<lt>zcbtfo4@acrm18E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
