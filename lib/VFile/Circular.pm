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


package VFile::Circular;

    use strict;
    use warnings;
    use P3Utils;
    use base qw(VFile);

=head1 Virus Description Processor for DTR Viruses

This is a subclass of L<VFile> that processes the C<checkv_circular.tsv> file.  This file contains unknown viruses.  The name is in the
C<original_id> field and the taxonomic ID is always C<10239>.

=head2 Special Methods

=head3 new

    my $vFile = VFile::Circular->new($dbDir);

=over 4

=item dbDir

The name of the directory containing the CheckV database.

=back

=cut

sub new {
    my ($class, $dbDir) = @_;
    my $retVal = {
        fileName => "$dbDir/genome_db/checkv_circular.tsv"
    };
    bless $retVal, $class;
    return $retVal;
}

=head2 Virtual Overrides

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
    my ($self, $ih) = @_;
    my (undef, $retVal) = P3Utils::find_headers($ih, virus_description => 'checkv_id', 'original_id');
    return $retVal;
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
    my ($self, $fields) = @_;
    return [10239, $fields->[0]];
}

1;


