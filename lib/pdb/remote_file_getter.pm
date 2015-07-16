package pdb::remote_file_getter;
use Moose;
use Carp;
use LWP::UserAgent;

has 'pdb_code' => (is => 'rw', isa => 'Str');

sub get_pdb_file_data {
    my $self = shift;
    my $url  = $self->_build_URL();
    my $ua   = LWP::UserAgent->new();
    
    my $response = $ua->get($url);

    croak "Error trying to get pdb file remotely: " . $response->status_line()
        if ! $response->is_success();

    return $response->decoded_content();
}

sub _build_URL {
    my $self = shift;
    croak "You have not set a pdb code!" if ! defined $self->pdb_code();
    return "http://www.rcsb.org/pdb/files/" . uc ($self->pdb_code) . ".pdb";
}

__PACKAGE__->meta->make_immutable;

1;
