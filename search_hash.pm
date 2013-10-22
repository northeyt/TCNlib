package search_hash;

use Moose;
use Moose::Util::TypeConstraints;
use types;

use Data::Dumper;

use Carp;

# subtypes

# attributes

has 'hash' => (
    is => 'rw',
    isa => 'HashRef',
    required => 1,
);

has 'depth' => (
    is => 'rw',
    isa => 'Int',
    required => 1,
);

has 'term_hash' => (
    is => 'rw',
    isa => 'HashRef',
    required => 1,
);

has 'wildcard_value' => (
    is => 'rw',
    isa => 'Str',
    default => '',
);


# methods

sub search {
    my $self = shift;

    my $max_depth = $self->depth;
    my $curr_depth = 0;
    my $wildcard = $self->wildcard_value;
    
    my @results = ();

    @results = _parse_level( $self->hash, $self->term_hash,
                             $curr_depth, $max_depth,
                             [ @results ], $wildcard );

    print "@results\n";
    return @results;
    
}

sub _parse_level {
    my( $subj_hash, $term_hash, $curr_depth,
        $max_depth, $results,   $wildcard ) = @_;
    
    ++$curr_depth;
    
    print "Parsing hash level $curr_depth ... \nTerm hash\n";
    print Dumper $term_hash;
    print "Subject hash\n";
    print Dumper $subj_hash;
    
    # Return all values if at max depth of hash
    if ( $max_depth == $curr_depth ) {

        print "Max depth reached: returning all values. \n";

        my @values = ();
        
        if ( $term_hash eq $wildcard ) {
            print "Wildcard at max depth - returning all values\n";

            @values = ( values %{ $subj_hash } );
            
            push( @{ $results }, @values );
        }
        else {

            my @final_depth_keys = ();

            if ( ref $term_hash eq 'ARRAY' ) {
                push(@final_depth_keys, @{ $term_hash } );
                
            }elsif ( ! ref $term_hash ){
                push(@final_depth_keys, $term_hash) ;
            }
            
            foreach my $key ( @final_depth_keys ) {
                croak "Wildcard value found within subject hash"
                    if exists $subj_hash->{$wildcard};

                if ( exists $subj_hash->{$key} ) {
                    push( @{ $results }, $subj_hash->{$key} );
                }
                
            }
        }
    }
    else { # Parse subject hash for hash refs that match term hash

        my @matches = ();

        foreach my $key ( keys %{ $term_hash } ) {

            if ( $key eq $wildcard ) {

                if ( ! exists $subj_hash->{$key} ){
                    
                    # Send all values of subj hash to sub
                    foreach my $value ( keys %{ $subj_hash} ) {

                        next if ref $subj_hash->{$value} ne 'HASH';
                        
                        print "Sending all keys of subj hash to sub\n";
                        
                        _parse_level( $subj_hash->{$value},
                                      $term_hash->{$wildcard},
                                      $curr_depth,
                                      $max_depth,
                                      $results,
                                      $wildcard, );   
                    }
                }
                else {
                    croak "Wildcard value found in subject hash\n"
                        . "Please specify a different wildcard value";
                }
                
            } elsif ( exists $subj_hash->{$key} ) {
                print "Term $key found in subject hash\n";
                
                _parse_level( $subj_hash->{$key},
                              $term_hash->{$key},
                              $curr_depth,
                              $max_depth,
                              $results,
                              $wildcard, );
            }
        }        
    }

    return @{ $results };
}


__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

search_hash - Perl extension for blah blah blah

=head1 SYNOPSIS

   use search_hash;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for search_hash, 

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
