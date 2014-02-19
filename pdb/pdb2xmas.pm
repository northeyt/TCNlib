package pdb::pdb2xmas;

use Moose;
use File::chdir;
use IO::CaptureOutput qw(capture_exec);
use Carp;

use TCNPerlVars;
use types;
use write2tmp;

has 'pdb_file' => (
    is => 'rw',
    isa => 'FileReadable',
    required => 1,
);

has 'pdb2xmas_file' => (
    is => 'rw',
    isa => 'FileExecutable',
    default => $TCNPerlVars::pdb2xmas,
);

has 'process_bin' => (
    is => 'rw',
    isa => 'Directory',
    default => $TCNPerlVars::pdb2xmas_bin,
);

has 'log' => (
    is => 'rw',
    isa => 'Str',
);

foreach my $process ( qw( solv ss hb  ) ) {
    has $process . '_file' => (
        is => 'rw',
        isa => 'FileExecutable',
        default => $TCNPerlVars::pdb2xmas_bin . $process,
    );
}

has 'last_output' => (
    is => 'rw',
    isa => 'ArrayRef',
);


# Methods

sub output {
    my $self = shift;
    
    my $abs_pdb_path = File::Spec->rel2abs( $self->pdb_file() );
    local $CWD = $self->process_bin();

    my $pdb2xmas = $self->pdb2xmas_file();
    my $solv = $self->solv_file();
    my $ss   = $self->ss_file();
    my $hb   = $self->hb_file();
    
    my $cmd = "$pdb2xmas $abs_pdb_path | $solv | $ss | $hb";

    my ($stdout, $stderr) = capture_exec( $cmd );
    
    croak "pdb2xmas produced no output. Command run:\n$cmd\n"
        if ! $stdout;
    
    $self->log($stderr);
 
    my @output = split ( /(?<=\n)/, $stdout );

    $self->last_output( [ @output ] );
    
    return @output;
}


1;
__END__

=head1 NAME

pdb::pdb2xmas - class to run AMs pdb2xmas program

=head1 SYNOPSIS

   use pdb::pdb2xmas;

   # Get array of pdb2xmas output lines
   $pdb2xmas->new(pdb_file => $pdb_file)->output();

   # Get array ref to output of last pdb2xmas process
   $output_arr = $pdb2xmas->last_output();

   # Get string of any stderr from last process
   $stderr_str = $pdb2xmas->log();

=head1 DESCRIPTION

Stub documentation for pdb::pdb2xmas, 

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

Copyright (C) 2014 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
