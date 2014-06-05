#!/acrm/usr/local/bin/perl
package pdb::idabchain;

use Moose;
use Carp;
use IO::CaptureOutput qw(capture_exec);
use TCNPerlVars;
use types;
use TryCatch;

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
        return $outputHref->{$chainID};
    }
    else {
        croak "Chain $chainID was not identified in idabchain output!";
    }
}

# Runs idabchain and returns hashref to hash of form
#  chainID => ChainStatus
# here ChainStatus can be: Antigen, Light, or Heavy.
sub getOutput {
    my $self = shift;

    my $output = $self->_runExec();
    my $outputHref = $self->_parseOutput($output);

    return $outputHref;
}


# Run the exectutable on user-defined input
sub _runExec {
    my $self = shift;

    # Define env variables
    $ENV{'KABATALIGN'} = $self->kabatAlignDir();
    $ENV{'DATADIR'} = $self->kabatConcDir();
    
    my $inputFile = _getInputFile($self->input());
    my $exec = $self->execPath();
    my $cmd = "$exec $inputFile";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd);

    if (! $success ) {
        my $err = "idabchain run failed.\nCommand run: $cmd\nSTDERR: $stderr";
        croak $err;
    }
    
    return $stdout;
}

# This sub attempts to obtain an input file name from an input variable
sub _getInputFile {
    my $input = shift;

    my $inputFile = "";
    
    # Is input a pdb or chain object?
    try {
        $inputFile = $input->pdb_file();
    }
    catch ($err) {
        # Is input a file path?
        try{
            if (-e $input) {
                $inputFile = $input;
            }
            else {
                croak "Input file $input does not exist!";
            }
        }
        catch ($err) {
            croak $err;
        }
    };
    return $inputFile;
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
        chomp $line;
        my($chainID, $chainType) = $line =~ /^Chain \d+, (\w+): (.*)$/;
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
