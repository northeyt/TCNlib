package uniprot::uniprot;

use Moose;
use Moose::Util::TypeConstraints;
use TCNUtil::types;

use TCNUtil::local::error;

use LWP::Simple;
use Carp;

# Subtypes

# Attributes

has 'accession_code' => (
    is => 'rw',
    isa => 'ValidAC',
    predicate => 'has_accession_code',
);

has 'file_name' => (
    is => 'rw',
    isa => 'FileReadable',
    predicate => 'has_file_name',
);

has 'data_array' => (
    is => 'ro',
    isa => 'local::error|ArrayRef[Str]',
    lazy => 1,
    builder => '_get_entry_data',
);

# Array lines hashed by line id
has 'data_hash' => (
    is => 'ro',
    isa => 'local::error|HashRef',
    lazy => 1,
    builder => '_build_data_hash',
);

# Local mode can be used to try and get entry data locally - currently not
# implemented
has 'local_mode' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

has 'organism_NCBI_TaxID' => (
    is => 'ro',
    isa => 'local::error|ArrayRef[Int]',
    builder => '_build_organism_NCBI_TaxID',
    lazy => 1,
);

with 'TCNUtil::local::can_error';

# Methods

sub _build_organism_NCBI_TaxID {
    my $self = shift;

    if ( ref $self->data_hash eq 'local::error' ) {
        return $self->data_hash;
    }
    
    my @OX_line = @{ $self->data_hash->{OX} };

    croak "Can't get organism NCBI TaxID: no OX lines parsed from entry data"
        if ! @OX_line;

    my @TaxID = ();

    foreach (@OX_line ) {
        my @line_TaxID = m{ NCBI_TaxID= (\d+) }gxms;
        push(@TaxID, @line_TaxID);
    }

    croak "No TaxIDs parsed from OX lines" if ! @TaxID;

    return [ @TaxID ];
}

sub _get_entry_data {
    my $self = shift;

    # A file for the entry or an ac code must be specified
    if ($self->has_file_name) {
        open(my $fh, '<', $self->file_name)
            or die "Cannot not open file '" . $self->file_name . "' $!\n";

        my @data = <$fh>;

        croak "File '" . $self->file_name . "' contains no data!" if ! @data;
        
        return [ @data ];
    }
    elsif ($self->local_mode) {
        croak "Local mode has not been implemented yet!";
    }
    else {
        croak "If no file is specified for a uniprot entry you must specify "
            . "an Accession Code"
                if ! $self->has_accession_code;

        my $ac = $self->accession_code;
     
        my $url = "http://www.uniprot.org/uniprot/$ac" . '.txt' ;
        
        my $response = get($url);

        if ( ! $response ) {
            my $message = "Request failed for url $url";
            
            my $error
                = local::error->new(
                    message => $message,
                    type => 'url_request_failed',
                    data => { url => $url,
                              accession_code => $ac, } );

            $self->set_error( $error->id() => $error );

            return $error;
        }
        
        my @line = map { $_ . "\n" } split("\n", $response);

        croak "No data found in reponse from request for url $url"
            if ! @line;

        return [ @line ];
    }        
}
    
sub _build_data_hash {
    my $self = shift;

    if ( ref $self->data_array eq 'local::error' ) {
        return $self->data_array;
    }
    
    my %hash = ();

    my $SQ_flag = 0;
    
    foreach my $line ( @{ $self->data_array } ){
        my $line_id = substr($line, 0, 2);
    
        if ($line_id eq '  ' ) {
            if ($SQ_flag) {
                push( @{ $hash{SQ} }, $line );
            }
            else {
                croak "Badly formatted line '$line' found in line array";
            }
        }
        elsif ($line_id eq '//') {
            last;
        }
        else{
            if ( ! exists $hash{$line_id} ) {
                $hash{$line_id} = [ $line ];
            }
            else {
                push( @{ $hash{$line_id} }, $line );
            }

            if ($line_id eq 'SQ') {
                $SQ_flag = 1;
            }else {
                $SQ_flag = 0;
            }
            
        }
    }
    return { %hash };    
}

__PACKAGE__->meta->make_immutable;

1;
__END__

=head1 NAME

uniprot::uniprot - Module to access uniprot entries via LWP, with parsers
                   added as needed

=head1 SYNOPSIS

   use uniprot::uniprot;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for uniprot::uniprot, 

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
