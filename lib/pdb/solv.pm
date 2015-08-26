package pdb::solv;

use Moose;
use Carp;
use TCNUtil::types;
use TCNPerlVars;
use TCNUtil::GLOBAL qw(rm_trail);
use pdb::pdbFunctions;
use IO::CaptureOutput qw(capture_exec);
use File::Spec;
use File::Basename;
use Cwd;
use TCNUtil::write2tmp;
write2tmp->Cache_Limit(100);

# Attributes

has 'pdbsolvExec' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::pdbsolv,
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

has 'resid2RelASAHref' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
    lazy => 1,
);

# Methods

sub getOutput {
    my $self = shift;

    my ($outASAFile, $outRSAFile) = $self->runExec();

    my %atomSerialHash = $self->getAtomSerialHash($outASAFile);

    my %resid2RelASA = $self->getResid2RelASAHash($outRSAFile);
    $self->resid2RelASAHref(\%resid2RelASA);

    # Remove output ASA and RSA files
    unlink($outASAFile);
    unlink($outRSAFile);
    
    return \%atomSerialHash;
}

sub getResid2RelASAHash {
    my $self = shift;
    my $outRSAFile = shift;
    
   open(my $fh, "<", $outRSAFile) or die "Cannot open file $outRSAFile, $!";

   my %resid2RelASA = ();

   my $headerStr
       = "#       RESIDUE  AA   ACCESS  RELACC  SCACC   SCRELACC\n";

   my $reachedHeader = 0;
   
   while (my $line = <$fh>) {
       if ($line eq $headerStr) {
           $reachedHeader = 1;
           next;
       }

       next unless $reachedHeader;

       chomp $line;
       
       my $chainID = rm_trail(substr($line, 8,  2));
       my $resSeq  = rm_trail(substr($line, 10, 4));
                                                     
       my $resid   = $chainID . "." . $resSeq;

       $resid2RelASA{$resid} = {allAtoms => rm_trail(substr($line, 30, 7))};
   }   
   return %resid2RelASA;
}


sub getAtomSerialHash {
    my $self = shift; 
    my $outASAFile = shift;

    open(my $fh, "<", $outASAFile) or die "Cannot open file $outASAFile, $!";

    my %atomSerialHash = ();

    while (my $line = <$fh>) {
        next unless $line =~ /^(ATOM|HETATM)/;
        my($serial, $ASA) = $self->parseLine($line);
        $atomSerialHash{$serial} = $ASA;
    }
    return %atomSerialHash;
}

sub parseLine {
    my $self   = shift;
    my $line   = shift;
    my $serial = rm_trail(substr($line,  6, 5));
    my $ASA    = rm_trail(substr($line, 60, 6));
    return ($serial, $ASA);
}

sub runExec {
    my $self = shift;

    my $exec        = $self->pdbsolvExec;
    my $radiiFile   = $self->radiiFile;
    my $probeRadius = $self->probeRadius;
    
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());

    # Get basename of input file
    my($baseName, $dir, $ext)  = fileparse($inputFile, '\..*');

    # Ensure that inputFile has full path to file
    $inputFile = File::Spec->rel2abs($inputFile);
    
    my $outputASAFile = File::Spec->rel2abs($baseName) . ".asa";
    my $outputRSAFile = File::Spec->rel2abs($baseName) . ".rsa";
    
    my $cmd = "$exec -f $radiiFile -p $probeRadius -r $outputRSAFile "
        . "$inputFile > $outputASAFile";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd); 
   
    if (! $success) {
        if (exists write2tmp->Cache->{$inputFile}) {
            write2tmp->retain_file(file_name => $inputFile);
        }
        croak "pdbsolve failed: $stderr\nCmd run: $cmd\n";
    }
    elsif (! -e $outputASAFile) {
        croak "pdbsolv failed to create an output file. STDERR: $stderr\n";
    }
    return ($outputASAFile, $outputRSAFile);
}

1;
__END__

=head1 NAME

pdb::solv - Perl extension for blah blah blah

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
