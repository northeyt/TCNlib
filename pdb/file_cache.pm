package pdb::file_cache;

use Moose;
use TCNPerlVars;
use pdb::remote_file_getter;

has 'cache_dir'   => (is => 'rw', isa => 'Str',
                      default => $TCNPerlVars::pdb_file_cache_dir);

has 'file_ext'    => (is => 'rw', isa => 'Str',
                      default => $TCNPerlVars::pdbext);

has 'file_prefix' => (is => 'rw', isa  => 'Str',
                      default => $TCNPerlVars::pdbprepname);

has 'pdb_code'    => (is => 'rw', isa => 'Str');

has 'remote_file_getter' => (
    is      => 'rw',
    isa     => 'pdb::remote_file_getter',
    lazy    => 1,
    default => sub {pdb::remote_file_getter->new()},
);

sub get_file {
    my $self = shift;

    my $file_name = $self->_build_file_name();

    if (! -e $file_name) {
        $self->_add_pdb_file_from_remote_to_cache();
    }
    return $file_name;
}

sub _build_file_name {
    my $self = shift;
    return $self->cache_dir . '/' . $self->file_prefix . lc ($self->pdb_code) . $self->file_ext;
}

sub _add_pdb_file_from_remote_to_cache {
    my $self = shift;
    
    $self->remote_file_getter->pdb_code($self->pdb_code);
    my $pdb_data = $self->remote_file_getter->get_pdb_file_data();
    $self->_write_pdb_data_to_cache($pdb_data);
}

sub _write_pdb_data_to_cache {
    my $self     = shift;
    my $pdb_data = shift;

    my $cache_file = $self->_build_file_name();
    open(my $OUT, ">", $cache_file) or die "Cannot open file $cache_file, $!";
    print {$OUT} $pdb_data;
}

__PACKAGE__->meta->make_immutable;

1;
