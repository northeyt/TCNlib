package FOSTA::FEPFinder;
use Moose::Role;
use TCNUtil::sequence;

requires 'getReliableFEPSequencesFromSwissProtID';

package FOSTA::Factory;
use Moose;

has 'remote' => (
    is       => 'rw',
    isa      => 'Bool',
    required => 1,
    default  => 0,
);

sub getFOSTA {
    my $self = shift;
    my @args = @_;

    if ($self->remote) {
        require TCNUtil::FOSTA::Remote;
        return FOSTA::Remote->new(@args);
    }
    else {
        require TCNUtil::FOSTA::Local;
        return FOSTA::Local->new(@args);
    }
}

1;
