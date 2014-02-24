package GLOBAL;

use strict; 
use warnings;
use Exporter;
use File::Temp;
use Carp;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = 1.00;
@ISA     = qw(Exporter);
@EXPORT  = ();
@EXPORT_OK = qw(array_is_identical array_is_subset is_int rm_trail one2three_lc three2one_lc);

### Modules

sub array_is_identical{
    my $arr1 = shift;
    my $arr2 = shift;

    if (scalar @{$arr1} eq scalar @{$arr2}) {
        foreach my $ele (@{$arr1}) {
            return 0 unless grep { /^$ele$/ } @{$arr2};
        }
        return 1;
    }
    else {
        return 0;
    }
}

sub array_is_subset {
    my $arr1 = shift;
    my $arr2 = shift or die "is_subset must be sent two array refs";
    
    my @sorted = sort { scalar @{$a} <=> scalar @{$b} } ($arr1, $arr2);

    foreach my $ele ( @{ $sorted[0] } ) {
        return 0 if ! grep /^$ele$/, @{ $sorted[1] };
    }

    return 1;
}

sub is_int {
    int($_[0]) eq $_[0] ? return 1 : return 0;
}

sub rm_trail {
    my $str = shift;
    
    $str =~ s{\A \s+|\s+ \z}{}gxms;

    return $str;
}

sub one2three_lc {
    my($onelc) = @_;

    $onelc = rm_trail($onelc);
    
    croak "1to3lc must be passed a single-character string"
        if ( ! $onelc || length $onelc != 1 );

    my %one2three
        = ( Q => 'GLN',
            W => 'TRP',
            E => 'GLU',
            R => 'ARG',
            T => 'THR',
            Y => 'TYR',
            I => 'ILE',
            P => 'PRO',
            A => 'ALA',
            S => 'SER',
            D => 'ASP',
            F => 'PHE',
            G => 'GLY',
            H => 'HIS',
            K => 'LYS',
            L => 'LEU',
            C => 'CYS',
            V => 'VAL',
            N => 'ASN',
            M => 'MET',
        );

    exists $one2three{ uc $onelc } ? return $one2three{ uc $onelc }
        : croak "$onelc is not a valid one-letter code";
     
}


sub three2one_lc {
    my($threelc) = @_;

    $threelc = rm_trail($threelc);
    
    croak "3to1lc must be passed a three-character string"
        if ( ! $threelc || length $threelc != 3 );

    my %three2one
        = ( GLN => 'Q',
            TRP => 'W',
            GLU => 'E',
            ARG => 'R',
            THR => 'T',
            TYR => 'Y',
            ILE => 'I',
            PRO => 'P',
            ALA => 'A',
            SER => 'S',
            ASP => 'D',
            PHE => 'F',
            GLY => 'G',
            HIS => 'H',
            LYS => 'K',
            LEU => 'L',
            CYS => 'C',
            VAL => 'V',
            ASN => 'N',
            MET => 'M',
        );

    exists $three2one{ uc $threelc } ? return $three2one{ uc $threelc }
        : croak "$threelc is not a valid three-letter code";

}

1;
__END__

=head1 NAME

GLOBAL - Common subroutines

=head1 SYNOPSIS

   use GLOBAL;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for GLOBAL, 

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
