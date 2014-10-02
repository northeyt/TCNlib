package parse_how_file;

require Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw( parse_how_file );

use strict; 
use warnings;
use Carp;

sub parse_how_file {
    my ($file) = @_;

    croak "You must pass parse_how_file a .how file!" if ! $file;
    
    open(my $fh, '<', $file) or die "Cannot open file $file, $!";

    my $record;

    {
        local $/;
        $record = <$fh>;
    }

    croak "Nothing read from .how file $file" if ! $record;
    
    close $fh;

    my @entry_hash_arr = ();
    
    # Records are seperated by a blank line
    my @entries = split("\n\n", $record );
    
    foreach my $entry (@entries) {
        my($antigen_seqlength, $pdb_code, $antigen_chid,
           $antigen_seq, $contacts, @antibody_chainids)
            = _parse_entry($entry);

        push (@entry_hash_arr,
              { pdb_code => $pdb_code,
                antibody_chain_ids => \@antibody_chainids,
                antigen_chain_id   => $antigen_chid,
                antigen_seq_length => $antigen_seqlength, } );
    }
    return @entry_hash_arr;
}

# Returns array of ( antigen_seqlength, pdb_code, antigen_chainid,
# antigen_sequence_string, contacts_string, antibody_chainids )
# Contacts string corresponds to antigen sequence string,
# where '.' = no contact and [A-Z] = contacting antibody chain id.
sub _parse_entry {
    my ($entry) = @_;

    my ($ant_seqlength, $pdb_code, $ant_chainid, $antAndContactSeq)
        = $entry =~ m{ \s+ (\d+)     # Antigen seq length
                       \s  (\w+)     # PDB code
                       \.  (\S)      # Antigen chain id
                       \n  ([A-Z.\n]+) # Antigen + Contact Sequence w/ newlines
                 }xms;
    
    croak "Something went wrong trying to parse .how entry"
        if ! defined ($ant_seqlength && $pdb_code && $ant_chainid
                          && $antAndContactSeq);

    # Remove newlines from antigen and contact sequence
    foreach (\$antAndContactSeq) {
        ${$_} =~ s{\s}{}gxms;
    }

    # Antigen sequence length determines the substring length of the
    # antigen and contact sequence that corresponds to the antigen sequence
    my $antigenSeq = substr($antAndContactSeq, 0, $ant_seqlength);
    my $contactSeq = substr($antAndContactSeq, $ant_seqlength);
    
    # Get antibody chains
    my @contacting_chains = $contactSeq =~ m{[A-Z]}gxms;

    croak "No antibody chains parsed from .how entry"
        if ! @contacting_chains;
    
    my %ab_chains = map {$_ => 1} @contacting_chains;

    return($ant_seqlength, $pdb_code, $ant_chainid, $antigenSeq, $contactSeq,
            sort keys (%ab_chains)); 
}


1;

__END__

=head1 NAME

parse_how_file - Perl extension for parsing .how files. This record format

is used by DiscoTope team to describe datasets. parse_how_file parses
information on the antigen and antibody from these records.

=head1 SYNOPSIS

   use parse_how_file;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for parse_how_file, 

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
