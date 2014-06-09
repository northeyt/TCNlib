#!/acrm/usr/local/bin/perl
package pdb::kabatnum;

use Moose;
use Carp;
use IO::CaptureOutput qw(capture_exec);
use write2tmp;
use TCNPerlVars;
use types;
use pdb::pdb;
use pdb::pdbFunctions;

### Attributes ################################################################

# File path to kabatnum executable
has 'execPath' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::kabatnum,
    required => 1,
);

# Numbering scheme to use
has 'scheme' => (
    is => 'rw',
    isa => 'Str',
    default => 'kabat',
);

# Maps chosen scheme to executable opt flag
has 'scheme2option' => (
    is => 'ro',
    isa => 'HashRef',
    default => sub { {kabat => '-k',
                      chothia => '-c',
                      ichothia => '-a'} },
);

has 'CDRDefinitions' => (
    is => 'ro',
    isa => 'HashRef',
    default => sub { {kabat => {light => {map {$_ => 1} [24 .. 34, 50 .. 56,
                                                         89 .. 97]},
                                heavy => {map {$_ => 1} [31 .. 35, 50 .. 65,
                                                         96 .. 102]} } }
                 }
);

# User-defined input
has 'input' => (
    is => 'rw',
);

### Methods ####################################################################

# This method assigns numbering and CDR status to atoms of chain input
sub sequenceChain {
    my $self = shift;

    # Ensure that input is a chain object
    croak "sequenceChain must be passed an object with an atom array"
        if ref $self->input() ne 'chain';

    my $chain = $self->input();
    
    # Write atoms to tmp file
    my @atoms = map {"$_"} @{$self->input->atom_array()};
    my $w2t = write2tmp->new(data => \@atoms);

    # Temporarily set input as this tmp file
    my $origInput = $self->input();
    $self->input($w2t->file_name());
  
    my $output = $self->getOutput($w2t);

    # Reset input
    $self->input($origInput);
    
    my @atomLines = split("\n", $output);

    # Create numbered atoms and label corresponding input atoms
    $self->labelAtoms($chain, \@atomLines);
}

# This method labels CDR atoms of a chain, given atom lines
# output by kabat num (and a chain!)
sub labelAtoms {
    my($self, $chain, $atomLinesAref) = @_;

    # Get atom attribute name that you are going to set
    my $attr = $self->scheme() . "Seq";
    
    my $i = 0;
    my @inputAtoms = @{$chain->atom_array()};
    
    foreach my $line (@{$atomLinesAref}) {
        
        next unless  $line =~ /^HETATM|ATOM/;
        my $atom = atom->new(ATOM_line => $line);
        
        # Loop through input atoms until atom with corresponding
        # serial is found
        until ($atom->serial() == $inputAtoms[$i]->serial()){
            ++$i;
        }
        
        $inputAtoms[$i]->$attr($atom->resSeq());
        
        # Label atom as CDR if resSeq is found in CDR definition
        if ($self->isAtomCDR($inputAtoms[$i], $chain->is_ab_variable())){
            $inputAtoms[$i]->is_CDR(1);
        }
    }
}

# This method checks if an atom is a CDR atom. Must be passed an atom
# to check and a chain type (light or heavy) 
sub isAtomCDR {
    my ($self, $atom, $chainType) = @_;
    
    # Get numbering scheme and atom attribute to check
    my $scheme = $self->scheme();
    my $attr = $scheme . 'Seq';

    # Get CDR resSeqs hash
    my $resSeqsHref = $self->CDRDefinitions->{$scheme}->{lc $chainType};

    # Get numeric characters of atom's resSeq to check against
    my $atomRSNumber = $atom->$attr() =~ /(\d+)/g;

    if (exists $resSeqsHref->{$atomRSNumber}){
        return 1;
    }
    else {
        return 0;
    }
}
    
# Wrapper for _runExec
sub getOutput {
    my $self = shift;
    
    my $output = $self->_runExec();

    return $output;
}

# Runs kabatnum and returns output
sub _runExec {
    my $self = shift;
        
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());
    my $exec = $self->execPath();
    my $opt = $self->scheme2option()->{$self->scheme()};
    
    my $cmd = "$exec $opt $inputFile";

    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd);

    if (! $success ) {
        my $err = "kabatnum run failed.\nCommand run: $cmd\nSTDERR: $stderr";
        croak $err;
    }
    return $stdout;
}


1;
__END__

=head1 NAME

pdb::kabatnum - Perl extension to run kabatnum

=head1 SYNOPSIS

   use pdb::kabatnum;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::kabatnum, 

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
