package cdhit;
use Moose;
use Carp;
use IO::CaptureOutput qw(capture_exec);
use File::Basename;
use types;
use TCNPerlVars;
use TryCatch;
use write2tmp;
use XML::DOM;

### Attributes

has 'execPath' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::cdhit,
);

has 'clstr2xmlExecPath' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::clstr2xml,
);

has 'input' => (
    is => 'rw',
);

has 'seqIDThreshold' => (
    isa => 'Num',
    is => 'rw',
    default => 0.9,
);

### Methods

sub getClusters {
    my $self = shift;

    # Get cdhit output cluster file
    my $outClusterFile = $self->_runExec();

    # Run output through clstr2xml
    my $c2XMLExec = $self->clstr2xmlExecPath();
    my $clusterXMLStr = `$c2XMLExec $outClusterFile`;

    my @clusters = $self->_processClusterXML($clusterXMLStr);

    return @clusters;
}

sub _processClusterXML {
    my $self = shift;
    my $clusterXMLStr = shift;

    my $parser = new XML::DOM::Parser;
    my $doc = $parser->parse($clusterXMLStr);
    
    my $clstrEleName = "representativeMember";
    my @clusterEles = $doc->getElementsByTagName($clstrEleName);

    my @clusters = ();

    foreach my $clusterEle (@clusterEles) {
        my @cluster = ();
        foreach my $member ($clusterEle->getElementsByTagName('member')) {
            my $seqName = $member->getFirstChild->getNodeValue();
            my $seqID = $member->getAttribute('Identity');
            push(@cluster, [$seqName, $seqID]);
        }
        push(@clusters, \@cluster);
    }
    return @clusters;
}

sub _getFASTAFileFromInput {
    my $self = shift;

    my @FASTAStrs = ();

    # Is input a file?
    if (-e $self->input()) {
        return $self->input();
    }

    # In input an array ref?
    if (ref $self->input() eq 'ARRAY'){
        my @FASTAStrs = ();
        foreach my $ele (@{$self->input()}) {
            try {
                my $FASTAStr = $ele->getFASTAStr();
                push(@FASTAStrs, $FASTAStr);
            }
            catch {
                # Is input element a string?
                if (! ref $ele) {
                    # Assume string is a FASTA Str
                    push(@FASTAStrs, $ele . "\n");
                }
                else {
                    croak "Could not get FASTA String from input: $ele";
                }
            };
        }

        # Write FASTA strings to a file
        my $w2t = write2tmp->new(data => \@FASTAStrs, suffix => ".fasta");
        return $w2t->file_name();
    }
    else {
        croak "Unable to process input " . $self->input();
    }
}

sub _runExec {
    my $self = shift;

    my $execPath = $self->execPath;
    my $seqThresh = $self->seqIDThreshold();
    my $FASTAFile = $self->_getFASTAFileFromInput();
    my $FASTAFileBName = basename($self->_getFASTAFileFromInput());
    
    my $outFile = '/tmp/' . $FASTAFileBName . ".out";
    my $outClusterFile = $outFile . ".clstr";
    
    my $cmd = "$execPath -c $seqThresh -i $FASTAFile -o $outFile";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd);

    if (! $success) {
        my $err = "cdhit run failed.\nCommand run: $cmd\nSTDERR: $stderr";
        croak $err;
    }
    return $outClusterFile;
}

__PACKAGE__->meta->make_immutable;

1;
