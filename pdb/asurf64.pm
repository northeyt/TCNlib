package pdb::asurf64;

use Moose;
use Carp;
use types;
use TCNPerlVars;
use GLOBAL qw(rm_trail);
use pdb::pdbFunctions;
use IO::CaptureOutput qw(capture_exec);
use File::Spec;
use File::Basename;
use Cwd;

# Attributes

has 'asurf64exec' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::asurf64,
    required => 1
);

has 'radiiFile' => (
    is => 'rw',
    isa => 'FileReadable',
    default => $TCNPerlVars::radii_file,
    required => 1,
);

has 'probeRadius' => (
    is => 'rw',
    isa => 'Num',
    default => 1.40,
    required => 1,
);

has 'input' => (
    is => 'rw',
);

# Methods

sub getOutput {
    my $self = shift;

    my $outFile = $self->runExec();

    open(my $fh, "<", $outFile) or die "Cannot open file $outFile, $!";

    my %atomSerialHash = ();

    while (my $line = <$fh>) {
        my($serial, $ASA, $radius) = $self->parseLine($line);
        $atomSerialHash{$serial} = [$ASA, $radius];
    }
    return \%atomSerialHash;
}

sub parseLine {
    my $self = shift;
    my $line = shift;
    
    my $serial = rm_trail(substr($line,  6, 5));
    my $ASA    = rm_trail(substr($line, 54, 8));
    my $radius = rm_trail(substr($line, 64, 4));

    return ($serial, $ASA, $radius);
}

sub runExec {
    my $self = shift;

    my $exec = $self->asurf64exec;
    my $radiiFile = $self->radiiFile;
    my $probeRadius = $self->probeRadius;
    
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());

    # Get basename of input file
    my($baseName, $dir, $ext)  = fileparse($inputFile, '\..*');

    # Ensure that inputFile has full path to file
    $inputFile = File::Spec->rel2abs($inputFile);
    
    # Get current working directory, then change to temp dir
    my $oldCwd = cwd();
    my $newCwd = "/tmp/";
    chdir($newCwd);

    my $cmd = "$exec -r $radiiFile -p $probeRadius -h $inputFile";
    my $outputFile = File::Spec->rel2abs($baseName) . ".asa";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd); 

    chdir($oldCwd);
   
    if (! $success) {
        croak "asurf64 failed: $stderr\n";
    }
    elsif (! -e $outputFile) {
        croak "asurf64 failed to create an output file";
    }
    return $outputFile;
}
1;
__END__

=head1 NAME

pdb::asurf64 - Perl extension for blah blah blah

=head1 SYNOPSIS

   use pdb::asurf64;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::asurf64, 

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

Copyright (C) 2014 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
