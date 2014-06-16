package cluster;

use Moose;
use Carp;

has 'array' => (
    isa => 'ArrayRef',
    is  => 'rw',
    required => 1,
);

has 'condition' => (
    isa => 'CodeRef',
    is  => 'rw',
    required => 1,
);

sub array_of_clusters {
    my $self = shift;

    my @array = @{ $self->array() };
    my $condition = $self->condition();

    my %clustered = ();
    my @clusters = ();

    for ( my $i = 0 ; $i < @array ; ++$i ){

        # Skip if i has already been clustered in prev loops
        next if exists $clustered{$i};
        
        my @cluster = ( $array[$i] );
        
        for ( my $j = $i + 1 ; $j < @array ; ++$j ) {
            # Compare cluster to j
            if ( &$condition( \@cluster, $array[$j] ) ) {
                
                # If condition is true, cluster j with i
                push( @cluster, $array[$j] );
                
                # Add index to indices to be skipped 
                $clustered{$j} = 1;
            }
        }
    
        # Add cluster to array of clusters
        push(@clusters, \@cluster);
    
    }
    return @clusters;
}

1;

__END__

=head1 NAME

cluster - Cluster some things!

=head1 SYNOPSIS

   use cluster;

   # My array of items to cluster
   $array_ref = [ qw ( some values to cluster ) ];

   # My subroutine that compares a cluster to an item. Returns TRUE if
   # item is to be added to cluster. FALSE if not.
   sub comparison_sub {
      my( $cluster_array, $item ) = @_;
          ...
   }

   # Send these to cluster ...
   cluster->new( array => $array_ref, condition => \&comparison_sub, );

   # To have your items clustered!
   my @clusters = cluster->array_of_clusters();

=head1 DESCRIPTION

cluster.pm supplies a generic method to cluster an array of items. The
items are clustered according to the descision made by the 

=head2 EXPORT

None by default.

=head1 SEE ALSO

=head1 AUTHOR

   # unique_grouping ensures that an item is only clustered
THOMAS NORTHEY, E<lt>tcn@dhcp-62-177.xlate.eduroam.ucl.ac.ukE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by THOMAS NORTHEY

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
