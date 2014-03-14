package pdb::get_files;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use TCNPerlVars;

use Carp;

# Subtypes

# Attributes

has 'pdb_code' => (
    is => 'rw',
    isa => 'ValidPDB',
    required => 1
);

for my $name ( 'pdb', 'xmas' ) {
    my $att_name = $name . '_file';
    my $builder = "_build_$name" . "_fname";
    
    has $att_name => ( is => 'ro',
                       isa => 'FileReadable',
                       builder => $builder,
                       lazy => 1,
                   );
    
}

# Methods

sub _build_pdb_fname {
    my $self = shift;

    my $pdbprep = $TCNPerlVars::pdbprep;
    my $pdbext  = $TCNPerlVars::pdbext;
    my $pdbdir  = $TCNPerlVars::pdbdir;
    
    my $pdb_code = $self->pdb_code();
    my $fname = $pdbprep . lc $pdb_code . $pdbext ;
    croak "No file for $pdb_code was found in $pdbdir. File name: $fname"
        if ! -e $fname;
    
    return $fname;
}
sub _build_xmas_fname {
    my $self = shift;
        
    my $xmasprep = $TCNPerlVars::xmasprep;
    my $xmasext  = $TCNPerlVars::xmasext;
    my $xmasdir  = $TCNPerlVars::xmasdir;
    
    my $pdb_code = $self->pdb_code();

    my $fname = $xmasprep . lc $pdb_code . $xmasext;
    croak "No file for $pdb_code was found in $xmasdir. File name: $fname"
        if ! -e $fname;
    
    return $fname; 
}



__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

pdb::get_files - Simple class to obtain locally stored pdb and xmas files for
a given pdb code.

=head1 SYNOPSIS

   use pdb::get_files;

   $get_my_file = pdb::get_files->new( pdb_code => $my_pdb_code );

   $pdb_file_name   = $get_my_file->pdb_code();
   $xmas_file_name  = $get_my_file->xmas_code();

=head1 DESCRIPTION

Stub documentation for pdb::get_files, 

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
