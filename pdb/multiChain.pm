package pdb::multiChain;

use strict; 
use warnings;

use pdb::asurf64;

# This method will assign ASAb values to the atoms of the input chains,
# where ASAb values are based on the surface accessibility values of the atoms
# when all input chains are considered.
# This should be used when you want consider the ASA of atoms in the context of
# a subset of the chains of of a pdb.
# Input: $chain->multiChainReadASA(\@otherChains)
sub readASAb {
    
    my $chainAref = shift;
    # Get array of all atoms from all chains
    my $atomAref = pdb::pdbFunctions::generateAtomAref(@{$chainAref});

    my $asurf = pdb::asurf64->new(input => $atomAref);
    my $atomSerial2ASARadHref = $asurf->getOutput();

    foreach my $atom (@{$atomAref}) {
        my ($ASA, $radius) = @{$atomSerial2ASARadHref->{$atom->serial()}};
        $atom->ASAb($ASA);
        $atom->radius($radius);
    }
}

# Given a ref to an array of chains, this functions return a ref to a hash
# where key = atom serial, value = atom for atoms from all chains
sub multiChainAtomSerialHref {
    my $chainAref = shift;

    my %atomSerialHash = ();

    my $atomAref = pdb::pdbFunctions::generateAtomAref(@{$chainAref});
    
    foreach my $atom (@{$atomAref}) {
        $atomSerialHash{$atom->serial()} = $atom;
    }
    return \%atomSerialHash;
}


1;
__END__

=head1 NAME

pdb::multiChain - Perl extension for blah blah blah

=head1 SYNOPSIS

   use pdb::multiChain;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::multiChain, 

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
