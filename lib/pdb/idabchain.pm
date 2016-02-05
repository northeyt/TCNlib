package pdb::idabchain;

use Moose;
use Carp;
use IO::CaptureOutput qw(capture_exec);
use TCNPerlVars;
use TCNUtil::types;
use TryCatch;
use TCNUtil::write2tmp;
use pdb::pdbFunctions;
use TCNUtil::local::error;

### Attributes ################################################################

# File path to idabchain executable
has 'execPath' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::idabchain,
    required => 1,
);

# Kabat alignment data dir
has 'kabatAlignDir' => (
    is => 'rw',
    isa => 'Directory',
    default => $ENV{'KABATALIGN'},
    required => 1,
);

has 'kabatConcDir' => (
    is => 'rw',
    isa => 'Directory',
    default => $ENV{'DATADIR'},
    required => 1,
);

# User-defined input
has 'input' => (
    is => 'rw',
);

### Methods ###################################################################

sub chainIs {
    my $self = shift;

    my $chainID = "";
    
    my $input = $self->input();
    try {
        $chainID = $input->chain_id();
    }
    catch ($err){
        croak "input '$input' has no chain_id! Is input definitely a chain?";
    };

    my $outputHref = $self->getOutput();

    if (exists $outputHref->{$chainID}) {
        my $result = $outputHref->{$chainID};
        if ($result eq 'Heavy' && $self->testForscFv($input)){
            return 'scFv';
        }
        else {
            return $result;
        }
    }
    else {
        my $errStr = "Chain $chainID was not identified in idabchain output!";
        my $newErr = local::error->new(message => $errStr,
                                       data => {chainID => $chainID},
                                       type => 'NoOutputForChainID');
    }
}

sub testForscFv {
    my $self = shift;
    my $chain = shift;

    if ($chain->chain_length() > 120 + 90) {

        my @truncated = $chain->seq_range_atoms(120, -1);

        my $w2t = write2tmp->new(data => [map {"$_"} @truncated]);

        $self->input($w2t->file_name());
        my $outputHref = $self->getOutput();
        if ($outputHref->{$chain->chain_id()} eq 'Light') {
            return 1;
        }
        else {
            return 0;
        }
    }
}

# Runs idabchain and returns hashref to hash of form
#  chainID => ChainStatus
# here ChainStatus can be: Antigen, Light, or Heavy.
sub getOutput {
    my $self = shift;

    my $output = "";
    
    try {
        $output = $self->_runExec();
    }
    catch ($err where {$_->type() eq 'runExecError'}){
        process_runExecError($err);
    };
    
    my $outputHref = $self->_parseOutput($output);

    return $outputHref;
}

# Processes error returned by _runExec by checking input file
sub process_runExecError {

    my $err = shift;
    
    # Investigate input file
    my $inputFile = $err->data->{inputFile};

    # Is input file empty?
    if (-z $inputFile) {
        # Construct error
        my $errStr =  "idabchain: Input file $inputFile is empty!";
        my $newErr = local::error->new(message => $errStr,
                                       parent => $err,
                                       data => {inputFile => $inputFile},
                                       type => 'EmptyInputFile');
        croak $newErr;
    }
    else {
        # Check to see if input file is all HETATMs
        open(my $FH, "<", $inputFile)
            or croak "Cannot open input file $inputFile";
        my $allHETATM = 0;
        while (my $line = <$FH>) {
            $allHETATM = 1 if $line =~ /^HETATM/;
            $allHETATM = 0 if $line =~ /^ATOM/;
        }

        my $errStr = "";

        if ($allHETATM) {
            $errStr
                = "Input file $inputFile is invalid - only HETATM lines present"
                    . " (inputfile must contain ATOM lines)";

            # Construct error
            my $newErr = local::error->new(message => $errStr,
                                           type => 'AllHETATMInputFile',
                                           data => {inputFile => $inputFile},
                                           parent => $err);
            croak $newErr;
        }
        else {
            # _runExec failed for an unknown reason, so return original error
            croak $err;
        }
    }
}


# Run the exectutable on user-defined input
sub _runExec {
    my $self = shift;

    # Define env variables
    $ENV{'KABATALIGN'} = $self->kabatAlignDir();
    $ENV{'DATADIR'} = $self->kabatConcDir();
    
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());
    
    my $exec = $self->execPath();
    my $cmd = "$exec $inputFile";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd);

    if (! $success) {
        # Retain temporary file for error checking
        if (exists write2tmp->Cache->{$inputFile}){
            write2tmp->retain_file(file_name => $inputFile);
        }
        
        my $errStr = "idabchain run failed.\nCmd run: $cmd\nSTDERR: $stderr";
        my $err = local::error->new(message => $errStr, type => 'runExecError',
                                    data => {inputFile => $inputFile});
        
        croak $err;
    }
    return $stdout;
}

# Given idabchain output, returns a hashref to hash of form
#  chainID => ChainStatus
# Where ChainStatus can be: Antigen, Light, or Heavy.
sub _parseOutput {
    my $self = shift;
    my $inputStr = shift;

    # Example inputStr:
    #  Chain 1, A: Antigen
    #  Chain 2, B: Antigen
    #  Chain 3, L: Light
    #  Chain 4, H: Heavy
    #  Chain 5, M: Light
    #  Chain 6, K: Heavy
                            
    my %hash = ();
    
    my @lines = split("\n", $inputStr);
    foreach my $line (@lines) {
        next if $line =~ /not found in matrix/; # Don't need these lines
        chomp $line;
        my($chainID, $chainType) = $line =~ /^Chain \d+, (\w+): (.*)$/;
        $chainType = 'scFv' if $chainType eq 'Heavy/Light';
        $hash{$chainID} = $chainType;
    }
    return \%hash;
}


1;
__END__

=head1 NAME

pdb::idabchain - Perl extension for running idabchain and parsing output

=head1 SYNOPSIS

   use pdb::idabchain;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::idabchain, 

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

=head1 BUGS

None reported... yet.

=cut
