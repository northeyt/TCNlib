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
use write2tmp;
write2tmp->Cache_Limit(100);

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
    default => $TCNPerlVars::asurf_radii_file,
    required => 1,
);

has 'standardDataFile' => (
    is => 'rw',
    isa => 'FileReadable',
    default => $TCNPerlVars::standardData,
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

has 'resid2RelASAHref' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
    lazy => 1,
);

# Methods

sub getOutput {
    my $self = shift;

    my ($outASAFile, $outRSAFile, $outLogFile) = $self->runExec();

    my %atomSerialHash = $self->getAtomSerialHash($outASAFile);

    my %resid2RelASA = $self->getResid2RelASAHash($outRSAFile);
    $self->resid2RelASAHref(\%resid2RelASA);

    # Remove output ASA and RSA files
    unlink($outASAFile);
    unlink($outRSAFile);
    unlink($outLogFile);
    
    return \%atomSerialHash;
}

sub getResid2RelASAHash {
    my $self = shift;
    my $outRSAFile = shift;
    
   open(my $fh, "<", $outRSAFile) or die "Cannot open file $outRSAFile, $!";

   my %resid2RelASA = ();

   my $headerStr
       = "REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL\n";

   my $reachedHeader = 0;
   
   while (my $line = <$fh>) {
       if ($line eq $headerStr) {
           $reachedHeader = 1;
           next;
       }

       next unless $reachedHeader;
       last if substr($line, 0, 3) eq 'END';

       chomp $line;
       
       my $chainID = rm_trail(substr($line, 7, 2));
       my $resSeq  = rm_trail(substr($line, 9, 5));
                                                     
       my $resid   = $chainID . "." . $resSeq;

       $resid2RelASA{$resid} = {allAtoms => rm_trail(substr($line, 22, 6))};
   }   
   return %resid2RelASA;
}


sub getAtomSerialHash {
    my $self = shift; 
    my $outASAFile = shift;

    open(my $fh, "<", $outASAFile) or die "Cannot open file $outASAFile, $!";

    my %atomSerialHash = ();

    while (my $line = <$fh>) {
        my($serial, $ASA, $radius) = $self->parseLine($line);
        $atomSerialHash{$serial} = [$ASA, $radius];
    }
    return %atomSerialHash;
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

    # Create symbolic link to standard data file (if not present already)
    my $link = "standard.data";
    unless (-e $link) {
        symlink($self->standardDataFile(), $link)
            or croak "asurf64.pm : unable to create symbolic link to standard data file";
    }
    
    my $cmd = "$exec -r $radiiFile -p $probeRadius -h $inputFile";
    
    my $outputASAFile = File::Spec->rel2abs($baseName) . ".asa";
    my $outputRSAFile = File::Spec->rel2abs($baseName) . ".rsa";
    my $outputLogFile = File::Spec->rel2abs($baseName) . ".log";
    
    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd); 

    chdir($oldCwd);
   
    if (! $success) {
        if (exists write2tmp->Cache->{$inputFile}) {
            write2tmp->retain_file(file_name => $inputFile);
        }
        croak "asurf64 failed: $stderr\nCmd run: $cmd\n";
    }
    elsif (! -e $outputASAFile) {
        croak "asurf64 failed to create an output file";
    }
    return ($outputASAFile, $outputRSAFile, $outputLogFile);
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
