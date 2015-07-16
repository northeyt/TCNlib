package BioBlastHitAdapter;

use Moose;
use Moose::Util::TypeConstraints;

subtype 'BioBlastHit',
    as 'Object',
    where {$_->can("hsp") && $_->hsp->can("hit_string")},
    message {"Object does not look like a Bio Blast hit object!'"};


has 'toAdapt' => (
    isa => 'BioBlastHit',
    is  => 'rw',
    required => 1,
);

# Allow calling like so: BioBlastHitAdapter::new($bioBlastObj) 
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    if (@_ == 1) {
        return $class->$orig(toAdapt => $_[0]);
    }
    else {
        return $class->$orig(@_);
    }
};

sub getFASTAStr {
    my $self = shift;

    my $id = $self->toAdapt->name();
    my $seqStr = $self->hsp->hit_string();

    # Remove any "-" characters from string that represent an alignment
    $seqStr =~ s/-//g;
    
    return "$id\n$seqStr";
}

package MSAligner;
use Moose::Role;
use Moose::Util::TypeConstraints;
use Carp;
use write2tmp;
use roles::consScoreCalculating;

requires 'align';
requires 'getAlignedSeqStringAref';

subtype 'canGetFASTAStr',
    as 'Object',
    where {$_->can("getFASTAStr")},
    message {"Must be able to perform getFASTAStr method on passed seqs!"};

coerce 'canGetFASTAStr',
    from 'BioBlastHit',
    via {BioBlastHitAdapter->new($_)};

has 'seqs' => (
    is      => 'rw',
    isa     => 'ArrayRef[canGetFASTAStr]',
    lazy    => 1,
    builder => '_buildSeqsFromOtherInput',
);

has 'consScoreCalculator' => (
    is => 'rw',
    does => 'roles::consScoreCalculating', # Can do consScoreCalculator
    predicate => 'hasConsScoreCalculator',
);

sub _buildSeqsFromOtherInput {
    croak "_buildSeqsFromOtherInput has not been implemented for this class!";
}

sub getInputIDs {
    my $self = shift;

    return map {seqStr::parseIDFromFASTA($_->getFASTAStr)} @{$self->seqs};
}

sub calculateConsScores {
    my $self = shift;
    
    croak "No consScoreCalculator supplied!"
        if ! $self->hasConsScoreCalculator;

    my @seqs
        = map {sequence->new($_)} $self->alignedSeqFASTStrs();
    
    $self->consScoreCalculator->seqs(\@seqs);
    
    return $self->consScoreCalculator->calcConservationScores();
}

sub alignedSeqFASTStrs {
    my $self = shift;
    my @seqIDs = $self->getInputIDs;
    my @seqs   = @{$self->getAlignedSeqStringAref};

    my @FASTAs = ();
    
    for (my $i = 0 ; $i < @seqIDs ; ++$i) {
        my $seqID = $seqIDs[$i];
        my $seq   = $seqs[$i];
        
        push(@FASTAs, ">" . $seqID . "\n" . $seq);
    }
    return @FASTAs;
}

package MSAexec;
use Moose::Role;
use Carp;
use sequence;

with ('roles::fileExecutor', 'MSAligner');

sub align {
    my $self = shift;
    return $self->runExec();
}

has 'inputSeqsStr' => (
    isa => 'Str',
    is  => 'rw',
    lazy => 1,
    builder => '_buildStrFromSeqs',
);

has 'inputSeqsFile' => (
    isa => 'FileReadable',
    is  => 'rw',
    lazy => 1,
    builder => '_writeSeqsStr2File',
);

before 'align' => sub {
    my $self = shift;
    
    if (@{$self->seqs} < 2) {
        croak "Only one sequence supplied: you must supply more than one sequence to align!";
    }
};

sub _writeSeqsStr2File {
    my $self = shift;

    # Write string to temporary file
    return write2tmp->new(data => [$self->inputSeqsStr])->file_name();
}

sub _buildStrFromSeqs {
    my $self = shift;
    return join("\n", map {$_->getFASTAStr()} @{$self->seqs});
}

# If seqs has not been supplied directly (but rather in string or file form),
# attempt to reverse build seqs from alternative input
sub _buildSeqsFromOtherInput {
    my $self = shift;

    # Try to parse input fasta sequences from inputSeqsFile
    # The inputSeqsFile attribute builder mean that we only have to try and
    # parse from file.
    my @FASTASeqs = parseFASTASeqsFromFile($self->inputSeqsFile());
    return [map {sequence->new($_)} @FASTASeqs];
}

sub parseFASTASeqsFromFile {
    my $inputFile = shift;
    open(my $IN, "<", $inputFile)
        or die "Cannot open $inputFile, $!";

    my $fData;
    {
        local $/;
        $fData = <$IN>;
    }
    return split(/(?=>)/, $fData);
}

package MSA::Clustal;
use Moose;
use Carp;
use File::Basename;
use Bio::SeqIO;

with 'MSAexec';

# Build clustal format cmd
sub cmdStringFromInputs {
    my $self = shift;

    my $inputFile = $self->inputSeqsFile();
    my $outputFile = $self->getTmpOutputFile();
    my $exec = $self->execFilePath();

    my $flags = $self->getFlags;
    my $opts  = $self->getOpts;

    return "$exec -OUTFILE=$outputFile -INFILE=$inputFile";
}

sub getTmpOutputFile {
    my $self = shift;  
    my $inputFile = $self->inputSeqsFile;
    my $inputFileBaseName = basename($inputFile);

    # TODO: create output file through write2tmp
    my $outputFile = '/tmp/' . $inputFileBaseName . '.aln';

    # Assign output file to tmp Cache
    write2tmp->Cache->{$outputFile} = 1;
    return $outputFile;
}

sub getAlignedSeqStringAref {
    my $self = shift;

    if (! eval {$self->align}) {
        if ($@ =~ /Only one sequence supplied/) {
            # Only one sequence has been supplied, so simply parse sequence
            # from input file
            return [map {$_->string()} @{$self->seqs()}];
        }
        else {
            croak $@;
        }
    }
    else {
        # Parse aligned sequence strings from output file
        return [$self->getAlignedSeqStringsFromOutFile()];
    }
}

sub getAlignedSeqStringsFromOutFile {
    my $self = shift;

    my $file2parse = $self->getTmpOutputFile();
    open(my $FH, "<", $file2parse)
        or croak "Cannot open file $file2parse";

    my %id2AlnStr = ();
    
    while (my $line = <$FH>) {
        # Skip header and any non-alignment lines
        next if $. < 4 || $line =~ /^\s/;

        my($id, $alnStr) = $line =~ /(\S+)\s+(\S+)/g;

        if (! exists $id2AlnStr{$id}) {
            $id2AlnStr{$id} = $alnStr;
        }
        else {
            $id2AlnStr{$id} .= $alnStr;
        }
    }
    # Return aligned sequence strings, ordered by order of input
    return map {$id2AlnStr{$_}} $self->getInputIDs();    
}

sub getSequenceStringsFromInputFile {
    my $self = shift;

    my $seqIO = Bio::SeqIO->new(-file => $self->inputSeqsFile(),
                                '-format' => 'Fasta');
    my @seqStrings = ();
    while (my $seq = $seqIO->next_seq) {
        push(@seqStrings, $seq->seq);
    }
    return @seqStrings;
}

package MSA::Clustalw;
use Moose;
use TCNPerlVars;

extends 'MSA::Clustal';

sub _buildExecPath {
    return $TCNPerlVars::clustalw;
}

package MSA::ClustalO;

use Moose;
use TCNPerlVars;

extends 'MSA::Clustal';

sub _buildExecPath {
    return $TCNPerlVars::clustalO;
}

package MSA::Muscle;
use Moose;
use Carp;

with 'MSAexec';

sub cmdStringFromInputs {
    my $self = shift;

    # Build muscle format cmd
    my $exec = $self->execFilePath();
    my $inputFile = $self->inputSeqsFile();

    # Ensure that quiet mode is on
    push(@{$self->flags}, "-quiet") if ! grep {/^-quiet$/} @{$self->flags};

    my $flags = $self->getFlags;
    my $opts  = $self->getOpts;
    
    return "$exec $flags $opts -in $inputFile";
}

sub _buildExecPath {
    return $TCNPerlVars::muscle;
}

sub getAlignedSeqStringAref {
    my $self = shift;

    croak "alignment was not successful"
        if ! $self->success();

    my %id2AlnStr = ();
    
    # Parse aligned sequence strings from output
    my %id2seq = $self->stdout =~ />(.*?)\n # ID line, starting with >
                                   (.*?)    # Sequence, including newlines
                                   (?=>|\z) # Up to next fasta entry /gxms;
    
    # Parse id and remove whitespace from each fasta sequence
    while (my ($id, $seq) = each %id2seq){
        # Remove whitespace from sequence
        $seq =~ s/\s//gms;
        $id2AlnStr{$id} = $seq;
    }

    # Return aligned sequence strings, ordered by order of input
    return [map {$id2AlnStr{$_}} $self->getInputIDs()];
}
