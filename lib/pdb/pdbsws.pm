package pdb::pdbsws::querier;
use Moose::Role;

requires 'getACsFromPDBCodeAndChainID';
requires 'getIDsFromPDBCodeAndChainID';
requires 'mapResSeq2SwissProtNum';

package pdb::pdbsws::Factory;
use Moose;

has 'remote' => (
    is       => 'rw',
    isa      => 'Bool',
    required => 1,
    default  => 0,
);

sub getpdbsws {
    my $self = shift;
    my @args = @_;

    if ($self->remote) {
        require pdb::pdbsws::Remote;
        return pdb::pdbsws::Remote->new(@args);
    }
    else {
        require pdb::pdbsws::Local;
        return pdb::pdbsws::Local->new(@args);
    }
}


1;
__END__

=head1 NAME

pdb::pdbsws - Perl extension for access to pdbsws, including methods for
              common searches

=head1 SYNOPSIS

   use pdb::pdbsws;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::pdbsws, 

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

Copyright (C) 2013 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
