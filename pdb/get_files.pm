package pdb::get_files;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use TCNPerlVars;

use Carp;
use TryCatch;
use local::error;
use pdb::file_cache;

has 'pdb_code' => (
    is => 'rw',
    isa => 'ValidPDB',
);

foreach my $name ('pdb', 'xmas') { 
    my $att_name = $name . '_file';
    my $builder = "_build_$name" . "_fname";
    
    has $att_name => ( is => 'ro',
                       isa => 'FileReadable',
                       builder => $builder,
                       lazy => 1,
                   );
}

# Attributes for finding pdb file in local depo
has 'pdbprepname' => (is => 'rw', isa => 'Str',
                      default => $TCNPerlVars::pdbprepname);
has 'pdbext'      => (is => 'rw', isa => 'Str',
                      default => $TCNPerlVars::pdbext);
has 'pdbdir'      => (is => 'rw', isa => 'Str',
                      default => $TCNPerlVars::pdbdir);

# Attributes for finding local pdb files
has 'xmasprep' => (is => 'rw', isa => 'Str', default => $TCNPerlVars::xmasprep);
has 'xmasext'  => (is => 'rw', isa => 'Str', default => $TCNPerlVars::xmasext);
has 'xmasdir'  => (is => 'rw', isa => 'Str', default => $TCNPerlVars::xmasdir);

has 'local_cache' => (
    is      => 'rw',
    isa     => 'pdb::file_cache',
    lazy    => 1,
    default => sub {pdb::file_cache->new()},
);

# Methods

sub _build_pdb_fname {
    my $self = shift;

    my $fname;
    
    try {
        $fname = $self->_fileFromLocalDepo();
    }
    catch (local::error $e where {$_->type eq 'FileNotFoundInLocalDepo'}) {
        $fname = $self->_fileFromLocalCache();
    };

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

sub _fileFromLocalDepo {
    my $self  = shift;
    my $fname = $self->pdbdir() . '/' . $self->pdbprepname()
        . lc $self->pdb_code() . $self->pdbext();

    if (! -e $fname) {
        my $message = "No file for " . $self->pdb_code() . " was found in "
            . $self->pdbdir() . ". File name: $fname";
        my $err = local::error->new(type => 'FileNotFoundInLocalDepo',
                                    data => $self, message => $message);
        croak $err;
    }
    return $fname;
}

sub _fileFromLocalCache {
    my $self = shift;

    $self->local_cache->pdb_code($self->pdb_code());
    return $self->local_cache->get_file();
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
