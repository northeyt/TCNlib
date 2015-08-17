package pdb::princip;

use strict; 
use warnings;

use Moose;
use Carp;
use pdb::pdbFunctions;
use IO::CaptureOutput qw(qxx);
use TCNUtil::types;

### ATTRIBUTES #################################################################
################################################################################

has 'execFile' => (
    is => 'rw',
    isa => 'FileExecutable',
    required => 1,
    lazy => 1,
    default => $TCNPerlVars::princip64,
);

has 'input' => (
    is => 'rw',
);

has 'output' => (
    is  => 'rw',
    isa => 'Str',
    predicate => 'has_output',
);

### METHODS ####################################################################
################################################################################


sub run {
    my $self = shift;

    my $pdbFile = pdb::pdbFunctions::getPDBFile($self->input);
    
    # input must be piped to princip
    my $cmdStr = "echo $pdbFile | " . $self->execFile();

    my ($stdout, $stderr, $success) = qxx($cmdStr);

    if ($success) {
        $self->output($stdout);
    }
    else {
        croak "princip was not successful: STDOUT: $stdout\nSTDERR: $stderr\n";
    }
}

sub getPlanarity {
    my $self = shift;

    $self->run if ! $self->has_output;
    
    if ($self->output =~ /RMS difference from best-fit plane:\s+(\S+)/){
        return $1;
    }
    else {
        croak "Could not parse planarity from princip output. OUTPUT: "
            . $self->output();
    }
}

1;


__END__

=head1 NAME

pdb::princip - run the princip program for pdb::pdb objects

=head1 SYNOPSIS

   use pdb::princip;
   

=head1 DESCRIPTION

The princip program does ...

=cut
