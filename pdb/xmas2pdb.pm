package xmas2pdb;

use Moose;
use Moose::Util::TypeConstraints;
use types;
use TCNPerlVars;
use write2tmp;

use Carp;

# subtypes

subtype 'ValidForm',
    as 'Str',
    where { $_ =~ m{ \A \s* ( monomer| multimer ) }xmsi },
    message { "$_ is not a valid Chain Id" };


# attributes

has 'radii_file' => (
    is => 'rw',
    isa => 'FileReadable',
    required => 1,
    default => $TCNPerlVars::radii_file,
);

has 'xmas2pdb_file' => (
    is => 'rw',
    isa => 'FileExecutable',
    required => 1,
    default => $TCNPerlVars::xmas2pdb,
);

has xmas_file => (
    is => 'rw',
    isa => 'FileReadable',
    predicate => 'has_xmas_file',
);

has all_lines => (
    is => 'rw',
    isa => 'Int',
    default => 0,
);

has form => (
    is => 'rw',
    isa => 'ValidForm',
    required => 1,
);

has 'chain_ids' => (
    is => 'rw',
    isa => 'ArrayRef[ValidChar]',
    predicate => 'has_chain_ids',
);

has output => (
    is => 'ro',
    isa => 'ArrayRef[Str]',
    lazy => 1,
    builder => '_run_xmas2pdb',
);

has output_file => (
    is => 'ro',
    isa => 'FileReadable',
    lazy => 1,
    builder => '_write_output',
);


# methods

sub _run_xmas2pdb {
    my $self = shift;

    croak "xmas file must be specified before output is read"
        if ! $self->has_xmas_file();

    my $xmas2pdb = $self->xmas2pdb_file();

    if ( -l $xmas2pdb ) {
        $xmas2pdb = readlink $xmas2pdb;
    }
    
    my $xmas_file = $self->xmas_file();
    my $radii_file = $self->radii_file();
    
    my $atoms_only = $self->all_lines ? '' : '-a' ;

    my $form = $self->form() eq 'monomer' ? '-m' : '-s';

    my $cmd = "$xmas2pdb $atoms_only $form -r $radii_file $xmas_file";

    my @output = `$cmd`;

    if ( $self->has_chain_ids() ) {
        my @filtered_output = ();
        foreach my $line (@output) {
            
            if ( $line =~ /^(ATOM|HETATM)/ ){
                my $chain = substr($line, 21, 1);
                
                next unless grep { /$chain/  }
                    @{ $self->chain_ids() };
            }
            push(@filtered_output, $line);
        }
        @output = @filtered_output;
    }
    
    croak "xmas2pdb produced no output. Attempted command:\n$cmd"
        if ! @output;

    return [ @output ];
    
}

sub _write_output {
    my $self = shift;

    my $file;
    
    if (@_) {
        $file = $_[0];
        open(my $fh, '>', $file)
            or die "Cannot write xmas2pdb output to $file, $!";

        print {$fh} @{ $self->output };
    }
    else {
        my $tmp = write2tmp->new( suffix => '.pdb', data => $self->output );
        $file = $tmp->file_name;
    }
    return $file;
}

__PACKAGE__->meta->make_immutable;


1;
__END__

=head1 NAME

xmas2pdb - Perl extension for blah blah blah

=head1 SYNOPSIS

   use xmas2pdb;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for xmas2pdb, 

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
