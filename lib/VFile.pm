#
# Copyright (c) 2003-2019 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


package VFile;

    use strict;
    use warnings;
    use P3Utils;

=head1 CheckV Virus Description File

This is the base class for a virus description file.  It takes as input a list of virus IDs and populates a hash mapping
them to [taxID, name] pairs.

The subclass constructor takes as input the name of the CheckV database directory.  The subclass must also
override the L</processLine> and L</processHeader> methods.

The object needs the following fields.

=over 4

=item fileName

Name of the file containing the virus definition data.

=back

=head2 Processing Methods

=head3 Process

    $vFile->Process(\@virusList, $retHash);

Read the virus description file and fill the identifying information for the listed viruses into the return hash.

=over 4

=item virusList

Reference to a list of virus IDs.

=item retHash

Reference to a hash that maps virus IDs to [taxonID, name] pairs.

=back

=cut

sub Process {
    my ($self, $virusList, $retHash) = @_;
    # Get a hash of the virus IDs.
    my %vidHash = map { $_ => 1 } @$virusList;
    # Open the file and parse the header.
    open(my $ih, '<', $self->{fileName}) || die "Could not open virus description file $self->{fileName}: $!";
    my $cols = $self->processHeader($ih);
    # Loop through the file, searching for viruses.
    while (! eof $ih) {
        my ($virusID, @fields) = P3Utils::get_cols($ih, $cols);
        if ($vidHash{$virusID}) {
            $retHash->{$virusID} = $self->processLine(\@fields);
        }
    }
    # Verify we found everything.
    my @missing = grep { ! exists $retHash->{$_} } @$virusList;
    if (scalar @missing) {
        print STDERR "WARNING:  Missing virus definitions for " . join(", ", @missing) . ".\n";
        for my $missing (@missing) {
            $retHash->{$missing} = [10239, "Unknown virus $missing"];
        }
    }
}

=head2 Virtual Methods

=head3 processHeader

    my $cols = $vFile->processHeader($ih);

Parse the headers and return the IDs of the useful columns.

=over 4

=item ih

Open file handle for the input file, positioned on the header record.

=item RETURN

Returns a reference to a list of IDs for the useful columns.  The first must identify the column containing the virus ID.

=back

=cut

sub processHeader {
    die "VFile::processHeader method not overridden.";
}


=head3 processLine

    my [$taxID, $name] = $vFile->processLine($fields);

Process the input fields from a virus description line and return the taxonomic ID and the name.

=over 4

=item fields

Reference to a list of the useful fields in the current input line.

=item RETURN

Returns a reference to a 2-tuple consisting of (0) the virus taxonomic ID and (1) the virus name.

=back

=cut

sub processLine {
    die "VFile::processLine method not overridden.";
}


1;


