package MSA;

use Moose;
use types;
use Carp;
use TCNPerlVars;
use TryCatch;
use write2tmp;
use IO::CaptureOutput qw(capture_exec);
use File::Basename;

### Attributes #################################################################

has 'clustalwExec' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::clustalw
);

has 'input' => (
    is => 'rw',
);

has 'processedInput' => (
    is => 'ro',
    isa => 'ArrayRef[Str]',
    builder => '_processInput',
    lazy => 1,
);

has 'processedInputIDs' => (
    is => 'ro',
    isa => 'ArrayRef[Str]',
    builder => '_getInputIDs',
    lazy => 1,
);


### Methods ####################################################################

# This method aligns the sequences given in the input. The method returns an
# array of sequences, where each element corresponds to the alignment sequence
# for an input sequence. The sequences will be in the same order that was given
# in the input
sub align {
    my $self = shift;

    # Get input
    my $inputAref = $self->processedInput();

    # If array only contains one sequence, return that sequence
    # (no need to align!)
    if (scalar @{$inputAref} == 1) {
        return (seqStrFromFASTAStr($inputAref->[0]));
    }
    
    my $inputFile = writeInput2File($inputAref);

    my @output = $self->_getOutput($inputFile);

    return @output;
}

sub _processInput {
    my $self = shift;

    my @processedInput = ();
    
    if (ref $self->input() eq 'ARRAY') { 
        foreach my $ele (@{$self->input()}) {
            # If ele is a string, assume it is a FASTA seq string
            if (! ref $ele) {
                push(@processedInput, $ele);
            }
            else {
                # Else try and get FASTAStr from element
                try {
                    push(@processedInput, $ele->getFASTAStr());
                }
                catch {
                    croak "Could not get FASTA seq from input element $ele, $@";
                }
            }
        }
    }
    elsif (-e $self->input()) {
        # Assume that input is a file containing FASTA seqs
        push(@processedInput, parseFASTAFile($self->input()));
    }
    elsif (! ref $self->input()) {
        # Assume that input is a string containing multiple FASTA seqs
        push(@processedInput, splitFASTAsStr($self->input()));
    }

    return \@processedInput;
}

sub _getInputIDs {
    my $self = shift;

    my @IDs = ();
    
    foreach my $FASTA (@{$self->processedInput()}) {
        my ($id) = $FASTA =~ />(\S+)/g;
        push(@IDs, $id);
    }

    return \@IDs;
}


sub parseFASTAFile {
    my $fileName = shift;
    
    open(my $FH, "<", $fileName) or croak "Cannot open file $fileName";
    my $FASTAsString = join("", <$FH>);

    return splitFASTAsStr($FASTAsString);
}


sub splitFASTAsStr {
    my $FASTAsString = shift;
    
    my @FASTAs = split(/(?=>)/, $FASTAsString);

    return @FASTAs;
}

sub writeInput2File {
    my $inputAref = shift;

    my $w2t = write2tmp->new(data => $inputAref, suffix => ".fasta");

    return $w2t->file_name();
}

sub _getOutput {
    my $self = shift;
    my $inputFile = shift;

    my $rawOutputFile = $self->_runClustalw($inputFile);

    return $self->_parseClustalwOutput($rawOutputFile);
}


sub _runClustalw {
    my $self = shift;
    my $inputFile = shift;    
    my $exec = $self->clustalwExec();

    my $inputFileBaseName = basename($inputFile);
    my $outputFile = '/tmp/' . $inputFileBaseName . '.aln';
    
    my $cmd = "$exec -OUTFILE=$outputFile -INFILE=$inputFile";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd);
    
    if (! $success) {
        my $err = "clustalw run failed.\nCommand run: $cmd\nSTDERR: $stderr";
        croak $err;
    }

    return $outputFile;
}

sub _parseClustalwOutput {
    my $self = shift;
    my $rawOutputFile = shift;

    open(my $FH, "<", $rawOutputFile)
        or croak "Cannot open file $rawOutputFile";

    my %id2AlnStr = ();
    
    while (my $line = <$FH>) {
        # Skip header and blank lines
        next if $. < 4 || $line =~ /^\s+$/;

        my($id, $alnStr) = $line =~ /(\S+)\s+(\S+)/g;

        if (! exists $id2AlnStr{$id}) {
            $id2AlnStr{$id} = $alnStr;
        }
        else {
            $id2AlnStr{$id} .= $alnStr;
        }
    }

    my @orderedStrs = ();

    foreach my $id (@{$self->processedInputIDs()}) {
        push(@orderedStrs, $id2AlnStr{$id});
    }

    return @orderedStrs;
}

sub seqStrFromFASTAStr {
    my $FASTAStr = shift;

    my ($seq) = $FASTAStr =~ /> \S* \s+ (.*) /gxms;

    # Remove any whitespace from seq
    $seq =~ s/\s//g;

    return $seq;
}

__PACKAGE__->meta->make_immutable;

1;
__END__

=head1 NAME

alignment - Perl extension for blah blah blah

=head1 SYNOPSIS

   use alignment;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for alignment, 

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
