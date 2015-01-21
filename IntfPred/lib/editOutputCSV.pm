package IntfPred::lib::editOutputCSV;

use strict; 
use warnings;

use Carp;

sub getLinesFromCSVFile {
    my $file = shift;

    open(my $IN, "<", $file) or die "Cannot open file $file, $!";

    my @lines = ();

    my $reachedHeader = 0;
    
    while (my $line = <$IN>) {
        if (! $reachedHeader) {
            $reachedHeader = 1 if $line =~ /^inst#/;
            next;
        }
        next if $line =~ /^\n$/;

        push(@lines, $line);
    }
    return @lines;
}

sub getCSVHeader {
    my $inputCSV = shift;

    open(my $CSV, "<", $inputCSV) or die "Cannot open file $inputCSV, $!";
    
    while (my $line = <$CSV>) {
        if ($line =~ /^inst#/){
            return $line;
        }
    }
    close $CSV;

    die "Did not parse header from CSV file!";
}

sub parseCSVLine {
    my $line = shift;
    chomp $line;
    my ($inst, $value, $prediction, $err, $score, $patchID) = split(",", $line);

    return [$patchID, $value, $prediction, $score];
}


1;
__END__

=head1 NAME

IntfPred::lib::editOutputCSV - Perl extension for blah blah blah

=head1 SYNOPSIS

   use IntfPred::lib::editOutputCSV;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for IntfPred::lib::editOutputCSV, 

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

Copyright (C) 2015 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
