package pdb::getresol;

use Moose;
use Carp;

use TCNPerlVars;
use types;
use local::error;
use IO::CaptureOutput qw( capture_exec );

# Attributes

has 'getresol_file' => (
    is => 'rw',
    default => $TCNPerlVars::getresol,
    required => 1,
);

has 'pdb_file' => (
    isa => 'FileReadable',
    is  => 'rw',
    predicate => 'has_pdb_file',
);

has 'experimental_method' => (
    isa => 'Str',
    is => 'rw',
);

has 'resolution' => (
    isa => 'Num',
    is  => 'rw',
);

has 'r_value' => (
    isa => 'Num',
    is  => 'rw',
);

# Methods

sub run {
    my $self = shift;

    croak "getresol: must set 'pdb_file' attribute before using run method"
        if ! $self->has_pdb_file();

    my $prog = $self->getresol_file();
    my $file = $self->pdb_file();

    my $cmd = "$prog $file";

    my ($stdout, $stderr, $success) = capture_exec($cmd);

    if ($success) {
        my($method, $res, $r_factor) = _parse_output($stdout);

        if ( defined $method && defined $res && defined $r_factor ) {
            $self->experimental_method( $method );
            $self->resolution( $res );
            $self->r_value( $r_factor / 100 );

            return 1;
        }
        else {
            my $message = "Unable to parse getresol output: '$stdout'\n";

            my $error = local::error->new( message => $message,
                                          type => 'parse_getresol_failed',
                                          data => {
                                              getresol_output => $stdout,
                                              cmd => $cmd,
                                          } );
            return $error;
        }
    }
    else {
        my $message
            = "Something went wrong trying to geteresol on file '$file'\n";

        my $error = local::error->new( message => $message,
                                       type => 'getresol_failed',
                                       data => { stderr => $stderr,
                                                 cmd => $cmd } );

        return $error;
    }

}

# Methods

sub _parse_output {
    my $output = shift;
    
    my($method, $resol, $r_fact)
        = $output =~ m{(\w+) , \s (\d+\.*?\d*?)A/(\d+\.*?\d*?)%}xms ;

    return($method, $resol, $r_fact);
}

1;
__END__

=head1 NAME

pdb::getresol - Perl wrapper to parse output of getresol program

=head1 SYNOPSIS

   use pdb::getresol;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::getresol, 

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
