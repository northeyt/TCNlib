#!/usr/bin/perl -w
# create_patch_files.pl --- Create patch files for a dataset of antigens
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 09 Jan 2014
# Version: 0.01

use warnings;
use strict;

use Carp;
use Data::Dumper;

use pdb::pdb;
use pdb::makepatch;
use pdb::xmas2pdb;
use write2tmp;

my $USAGE = <<EOF;
create_patch_files.pl:
    [patch radius] [pdb dir] [xmas dir] [patch.centres file]
EOF

my $patch_radius = shift or croak $USAGE;
$patch_radius == int ($patch_radius)
    or croak "Patch radius must be a positive integer";

my $pdb_dir = shift or croak $USAGE;
(-d $pdb_dir) or croak "$pdb_dir is not a directory";

my $xmas_dir = shift or croak $USAGE;
(-d $xmas_dir) or croak "$xmas_dir is not a directory";

my $patch_cent_file = shift or croak $USAGE;
(-e $patch_cent_file) or croak $USAGE; 

# Create patch directories if they do not exist already
# NOTE: script originally worked on multiple radii, now not the case
my %radii = ( $patch_radius => '', );

foreach my $radius (keys %radii) {
    my $dir_name = 'patches_' . $radius;

    if ( ! -d $dir_name ) {
        mkdir($dir_name) or die "Cannot create dir '$dir_name', $!";
    }
    $radii{$radius} = $dir_name;
}

open( my $IN, '<', $patch_cent_file )
    or die "Cannot open $patch_cent_file, $!";

my @entries = ();

{
    local $/ = ')';

    @entries = <$IN>;
}

# Avoid inclusion of trailing whitespace as an entry
if ( $entries[-1] =~ / \A \s+ \z /xms ) {
    pop @entries;
}

close $IN;

croak "No entries read from file $patch_cent_file!"
    if ! @entries;

foreach my $entry (@entries ) {

    $entry =~ s/\s+//g;

    my $pdb_code    = substr($entry, 0, 4);
    my $ag_chain_id = substr($entry, 5, 1);

    print "pdb: $pdb_code, antigen chain: $ag_chain_id ...";

    my @patch_centre_res = $entry =~ / $ag_chain_id : (\d+) /gxmsi;

    my $file_prefix = $pdb_code . "_$ag_chain_id";

    my $pdb_file  = "$pdb_dir/$file_prefix.pdb";
    my $xmas_file = "$xmas_dir/$file_prefix.xmas";
    
    my $chain = chain->new( pdb_code => $pdb_code,
                            chain_id => $ag_chain_id,
                            pdb_file  => $pdb_file,
                            xmas_file => $xmas_file,
                            hydrogen_cleanup => 1,
                            solvent_cleanup => 1,
                        );

    $chain->read_ASA();

    print "ASA has been read\n";
    
    my $replace = { tempFactor => 'ASAb', occupancy => 'radius' };
    
    # ATOM lines formatted for makepatch
    my @atom_lines
        = map { $_->stringify($replace) } @{ $chain->atom_array };
    
    my $mod_pdb = write2tmp->new(data => [@atom_lines])->file_name(); 

    # Open output files
    foreach my $radius ( keys %radii ) {
        my $dir = $radii{$radius};
        my $file_name = lc ($pdb_code) . uc ($ag_chain_id) . '.patches';
        open( my $OUT, '>>', "$dir/$file_name" );
        
        foreach my $resSeq (@patch_centre_res) {
            my $resid = "$ag_chain_id.$resSeq";
            
            my $cent_atom = $chain->highestASA($resid, 'ASAb');

            # Check residue of central atom also has a C-alpha atom
            # (required for makepatch)
            if (! exists $chain->resid_index->{$cent_atom->resid()}->{'C'}) {
                print " Residue " . $cent_atom->resid()
                    . " has no C-atom, skipping central atom "
                        . $cent_atom->serial() . "\n";

                next;
            }
            
            my $makepatch = makepatch->new( pdb_file => $mod_pdb,
                                            pdb_code => $pdb_code,
                                            central_atom => $cent_atom,
                                            patch_type => 'normal',
                                            radius => $radius,
                                        );
            
            my $patch = eval { patch->new($makepatch) };
            if ($@) {
                print Dumper $@;
                exit;
            }
            print {$OUT} $patch->summary();
        }
        close $OUT;
    }
    print "done\n";
}

__END__

=head1 NAME

create_patch_files.pl - Describe the usage of script briefly

=head1 SYNOPSIS

create_patch_files.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for create_patch_files.pl, 

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
