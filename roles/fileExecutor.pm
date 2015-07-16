package roles::fileExecutor;
use Moose::Role;
use types;
use IO::CaptureOutput qw( qxx );

requires 'cmdStringFromInputs';

has 'execFilePath' => (
    is => 'rw',
    isa => 'FileExecutable',
    required => 1,
    lazy => 1,
    builder => '_buildExecPath',
);

foreach (qw(stdout stderr)){  
    has $_ => (
        is => 'rw',
        isa => 'Str',
        predicate => 'has_' . $_,
    );
}

has 'success' => (
    is => 'rw',
    isa => 'Bool',
    lazy => 1,
    builder => 'runExec',
);

has 'flags' => (
    is => 'rw',
    isa => 'ArrayRef',
    default => sub { [] },
    lazy => 1,
);

has 'opts' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
    lazy => 1,
);

sub getFlags {
    my $self = shift;

    return join(" ", @{$self->flags});
}

sub getOpts {
    my $self = shift;

    my @optStrings = map {join(" ", $_, $self->opts->{$_})} keys %{$self->opts};

    return join(" ", @optStrings);
}

sub runExec {
    my $self = shift;

    my ($stdout, $stderr, $success) = qxx($self->cmdStringFromInputs());

    $self->stdout($stdout) if $stdout;
    $self->stderr($stderr) if $stderr;

    return $success;
}

1;
