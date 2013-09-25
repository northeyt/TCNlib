package types;

use Moose;
use Moose::Util::TypeConstraints;
use Carp;
use TCNPerlVars;
use File::Temp;
use write2tmp;

### File and Dir checks

subtype 'Directory',
    as 'Str',
    where { -d $_ },
    message { "$_ is not a directory" };


subtype 'FileReadable',
        as 'Str',
        where { -r $_ },
    message { "$_ is not a readable file!" };

# Prints array of strings to tmp file and then returns tmp file name
coerce 'FileReadable',
    from 'ArrayRef[Str]',
    via { my %arg = ( data => $_,
                      dir  => $TCNPerlVars::tmpdir, );

          #my $parent = ( caller(10) )[3];
          #print "$parent\n";
          my $tmp = write2tmp->new(%arg);
          
          return $tmp->file_name; 
};


subtype 'FileWritable',
        as 'Str',
        where { -w $_ },
    message { "$_ is not a writable file!" };

subtype 'FileExecutable',
    as 'Str',
    where { -X $_ },
    message { "$_ is not exectuable" };


### Subs for PDB code checks

subtype 'ValidChar',
    as 'Str',
    where { $_ =~ m{ \A [A-Z0-9] \z }xmsi },
    message { "$_ is not a valid pdb chain id " };


subtype 'ValidPDB',
    as 'Str',
    where { $_ =~ m{ \A \d [A-Z0-9]{3}  \z }xmsi },
    message { "$_ is not a valid PDB code" };

subtype 'ValidPDBID',
    as 'Str',
    where { $_ =~ m{ \A \d [A-Z0-9]{3} [A-Z] \z }xmsi },
    message{ "$_ is not a valid pdbid!" };

1;
__END__

=head1 NAME

types - Perl extension for blah blah blah

=head1 SYNOPSIS

   use types;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for types, 

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
