package rotate2pc;

use strict; 
use warnings;

use IO::CaptureOutput qw(capture);
    
use Exporter;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = 1.00;
@ISA     = qw(Exporter);
@EXPORT  = ();
@EXPORT_OK = qw(rotate2pc xyzmean);



use Math::MatrixReal;
use Math::VectorReal qw(:all);
use Math::Trig;
use Statistics::PCA;

use Carp;

sub rotate2pc {
    my(@vector) = @_;
    
    foreach (@vector) {
        croak "rotate2pc: $_ is not a vector"
            if ref $_ ne 'Math::VectorReal';
    }
    
    @vector = meancenter(@vector);

    my @pc = get_eigenvectors(@vector);

    my @axis = ( X, Y );

    my @R = ();
    
    for my $i ( 0 .. 1 ) {
        
        my ($plane_vector) = plane( O, $axis[$i], $pc[$i] );

        my $inner_product = innerproduct($axis[$i], $pc[$i] );
        my $angle = acos($inner_product);

        $R[$i] = RM_about_vector($plane_vector, $angle);

        if ($i == 0) {
            $pc[1] = $pc[1] * $R[0];
        }
    }
    return $R[0]->multiply($R[1]);
}


sub get_eigenvectors {

    my(@vector) = @_;

    foreach (@vector) {
        croak "get_eigenvectors: passed value $_ is not a vector"
            if ref $_ ne 'Math::VectorReal';
    }
    
    # Get pcs
    my $data = [ map { [ $_->x(), $_->y(), $_->z() ] } @vector  ];

    my $pca = Statistics::PCA->new();

    capture sub {$pca->load_data( { format => 'table', data => $data } ) };
    
    capture sub {$pca->pca( { eigen => 'M' } ) };
                                    
    return map { vector( $_->[0], $_->[1], $_->[2] ) }
        @{ $pca->results('eigenvector')};    
}

# Returns x y and z means for mean centering
sub xyzmean {
    my(@vectors) = @_;

    croak "No vectors passed to meancenter" if ! @vectors;
   
    my $N = scalar @vectors;
    
    # Means
    my %total = ( x => 0, y => 0, z => 0 );
    
    foreach my $v (@vectors) {
        
        croak "meancenter: passe value $v is not a vector object"
            if ref $v ne 'Math::VectorReal';
        
        foreach my $coord (keys %total) {
            $total{$coord} += $v->$coord;
        }
    }

    my %mean = ( x => $total{x} / $N,
                 y => $total{y} / $N,
                 z => $total{z} / $N,
             );
 
    return %mean;
    
}

sub meancenter {
    my(@vectors) = @_;

    croak "No vectors passed to meancenter" if ! @vectors;

    foreach (@vectors) {
        croak "meancenter: passe value $_ is not a vector object"
            if ref $_ ne 'Math::VectorReal';
    }

    my %mean = xyzmean(@vectors);
    
    my @centered_v = map { vector( ( $_->x - $mean{x} ),
                                   ( $_->y - $mean{y} ),
                                   ( $_->z - $mean{z} ), )
                       } @vectors;

    return @centered_v;
                               
}

sub innerproduct {
    my(@vector) = @_;

    croak "Two vectors must be passed to innerproduct"
        if @vector ne 2;

    foreach (@vector) {
        croak "innerproduct: passed value $_ is not a vector"
            if ref $_ ne 'Math::VectorReal';
        return 0 if $_->length == 0; # Avoid division by zero
    }

    my $dot_prod = $vector[0] . $vector[1];
    
    my $inner_prod = ( $vector[0]->norm . $vector[1]->norm );

    return $inner_prod;
}

sub RM_about_vector {
    my($vector, $angle) = @_;

    croak "RM_about_vector: passed value $vector is not a vector"
        if ref $vector ne 'Math::VectorReal';

    my $v = $vector->norm;

    my $a = ($v->x + $v->y + $v->z);

    my $f = 1 - cos($angle);

    my $R = Math::MatrixReal->new_from_rows( [
        [  cos($angle) + $v->x**2 * $f,
           $v->x * $v->y * $f  - $v->z * sin($angle),
           $v->y * sin($angle) + $v->x * $v->z * $f, ],
        [  $v->z * sin($angle) + $v->x * $v->y * $f,
           cos($angle) + $v->y**2 * $f,
          -$v->x * sin($angle) + $v->y * $v->z * $f, ],
        [ -$v->y * sin($angle) + $v->x * $v->z * $f,
           $v->x * sin($angle) + $v->y * $v->z * $f,
           cos($angle) + $v->z**2 * $f, ]
    ] );
                                             
    return $R;
          
}


1;
__END__

=head1 NAME

rotate2pc - Perl extension for blah blah blah

=head1 SYNOPSIS

   use rotate2pc;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for rotate2pc, 

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
