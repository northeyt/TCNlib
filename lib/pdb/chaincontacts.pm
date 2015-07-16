package pdb::chaincontacts;

use Moose;
use Carp;
use IO::CaptureOutput qw(capture_exec);
use write2tmp;
use TCNPerlVars;
use types;
use pdb::pdbFunctions;

### Attributes ################################################################

# File path to chaincontacts executable
has 'execPath' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::chaincontacts,
    required => 1,
);

# Distance threshold that defines if resiudes are contacting
has 'threshold' => (
    isa => 'Num',
    is => 'rw',
    default => 4,
);

# User-defined input
has 'input' => (
    is => 'rw',
);

### Methods ####################################################################

# Runs exec and returns results object
sub getOutput {
    my $self = shift;
    
    my $output = $self->_runExec();

    return $self->_parseOutput($output);
}

# This method parses output of chaincontacts to create a result object
sub _parseOutput {
    my $self = shift;
    my $output = shift;

    my @lines = split("\n", $output);

    my %contactHash = ();

    # Parse contact information from relevant lines
    foreach my $line (@lines) {
            
        if ($line =~ /^Chain/) {

            # Process line
            my ($chIDA, $resA, $chIDB, $resB, $contacts) = _parseLine($line);

            # Add contact to hash
            _addContactToHash(\%contactHash, $chIDA, $resA, $chIDB, $resB,
                             $contacts);

            # Reverse order of chains
            _addContactToHash(\%contactHash, $chIDB, $resB, $chIDA, $resA,
                             $contacts);
            
        }
    }

    my $result = pdb::chaincontacts::result->new(resultHash => \%contactHash);

    return $result;
}

sub _addContactToHash {
    my($contactHref, $chIDA, $resA, $chIDB, $resB, $contacts) = @_;
            
    # Create hash if not created already
    $contactHref->{$chIDA} = {}
        if ! exists $contactHref->{$chIDA};
    
    # Create array if not created already
    $contactHref->{$chIDA}->{$chIDB} = []
        if ! exists $contactHref->{$chIDA}->{$chIDB};
    
    my $contactArr = [$resA, $resB, $contacts];
    
    push(@{$contactHref->{$chIDA}->{$chIDB}}, $contactArr); 
}


sub _parseLine {
    my $line = shift;

    my ($chIDA, $resA, $chIDB, $resB, $contacts)
        = $line =~ m{Chain: \s* (\S+) \s*
                     Res: \s* (\S+) \s* - \s*
                     Chain: \s* (\S+) \s*
                     Res: \s* (\S+) \s*
                     Contacts: \s* (\d+)
                }gxms;
    
    return($chIDA, $resA, $chIDB, $resB, $contacts);
}


# Runs chaincontacts and returns output
sub _runExec {
    my $self = shift;
    
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());
    my $exec = $self->execPath();

    my $dist = $self->threshold();
    
    my $cmd = "$exec -r $dist $inputFile";
    
    my($stdout, $stderr, $success, $exit_code) = capture_exec($cmd);

    if (! $success ) {
        my $err = "chaincontacts run failed.\nCommand run: $cmd\nSTDERR: $stderr";
        croak $err;
    } 
    return $stdout;
}

__PACKAGE__->meta->make_immutable;

package pdb::chaincontacts::result;

use Moose;
use TryCatch;

### Attributes #################################################################

has 'resultHash' => (
    is => 'rw',
    isa => 'HashRef',
);

### Methods ####################################################################

# This method accepts two groups of chain ids and returns an arrayref of the
# resids from the SECOND group of chains that are in contact with residues
# the FIRST chain group. Chain IDs or chains can be input.
# e.g. $result->chain2chainContacts([$chainA, $chainB], [$chainC, $chainD]);
# Returns an REF to array of group 2 resids that in contact with group 1 chains
sub chain2chainContacts {
    my $self = shift;
    my($groupA, $groupB) = @_;

    my @groupA = @{$groupA};
    my @groupB = @{$groupB};
    
    # If any input is a chain or atom then get chain id
    foreach my $group (\@groupA, \@groupB) {
        for (my $i = 0 ; $i < @{$group} ; $i++) {
            try {
                # Chain
                $group->[$i] = $group->[$i]->chain_id();
            }
            catch  {
                try {
                    # Atom
                    $group->[$i] = $group->[$i]->chainID();
                }
                catch {
                    # No need to do anything, will assume var is a chain id
                    # string
                };
            };
        }
    }
    
    my %resids = ();

    # Loop through chains that have been supplied to test contacts for
    foreach my $chB (@groupB) {
        foreach my $chA (@groupA){
            # Get relevant contacts
            my $contactAref = $self->resultHash()->{$chA}->{$chB};
        
            # Hash contacts, as resSeqs from B can be in contact with multiple
            # resSeqs from A
            map { $resids{"$chB." . $_->[1]} = 1 } @{$contactAref};
        }
    }
        return [keys %resids];
}


__PACKAGE__->meta->make_immutable;

1;
__END__

=head1 NAME

pdb::chaincontacts - Perl extension to provide interface to chaincontacts

=head1 SYNOPSIS

   use pdb::chaincontacts;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::chaincontacts, 

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
