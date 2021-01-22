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


package VBin;

    use strict;
    use warnings;

=head1 Virus Bin Descriptor

This object contains data about a virus bin.  In addition to a list of contig IDs, it contains the total completeness, the weighted error,
the weighted coverage, the virus ID, and the total length.

The fields of this object are as follows.

=over 4

=item contigs

Reference to a list of the IDs for the contigs in this bin.

=item completeness

The total percent completeness representing how much of the virus was found.

=item total_error

The weighted sum of the error amount.

=item total_coverage

The weighted sum of the coverage.

=item virus_id

The ID of the virus found in the bin.

=item len

The total length of the virus in base pairs.

=item num

The bin ID number.

=back

=head2 Special Methods

=head3 new

    my $vBin = VBin->new($id);

Create a new, empty virus bin.

=over 4

=item id

The ID of the virus in this bin.

=back

=cut

sub new {
    my ($class, $id) = @_;
    my $retVal = {
        total_error => 0,
        completeness => 0,
        total_coverage => 0,
        virus_id => $id,
        len => 0,
        contigs => [],
        num => 0
    };
    bless $retVal, $class;
    return $retVal;
}

=head2 Update Methods

=head2 AddContig

    $vBin->AddContig($contigID, $len, $completeness, $error, $covg);

Add a contig to this bin.

=over 4

=item contigID

ID of the contig to add

=item len

length of the contig being added.

=item completeness

completeness percent of the contig being added

=item error

error percent of the contig being added

=item covg

coverage of this contig

=back

=cut

sub AddContig {
    my ($self, $contigID, $len, $completeness, $error, $covg) = @_;
    push @{$self->{contigs}}, $contigID;
    $self->{len} += $len;
    $self->{completeness} += $completeness;
    $self->{total_error} += $error * $len;
    $self->{total_coverage} += $covg * $len;
}

=head3 num

    my $num = $vBin->num;

or

    $vBin->num($newNum);

Return or set the ID number of this bin.

=cut

sub num {
    my ($self, $newNum) = @_;
    if (defined $newNum) {
        $self->{num} = $newNum;
    }
    return $self->{num};
}

=head2 Query Methods

=head3 contigs

    my $contigList = $vBin->contigs;

Return a reference to the list of contigs in this bin.

=cut

sub contigs {
    my ($self) = @_;
    return $self->{contigs};
}

=head3 len

    my $len = $vBin->len;

Return the number of base pairs in this bin.

=cut

sub len {
    my ($self) = @_;
    return $self->{len};
}

=head3 covg

    my $covg = $vBin->covg;

Return the mean coverage for this bin.

=cut

sub covg {
    my ($self) = @_;
    return ($self->{total_coverage} / $self->{len});
}

=head3 err

    my $err = $vBin->err;

Return the mean confidence error for the bin.

=cut

sub err {
    my ($self) = @_;
    return ($self->{total_error} / $self->{len});
}

=head3 percent

    my $percent = $vBin->percent;

Return the percent of the virus covered by this bin.

=cut

sub percent {
    my ($self) = @_;
    return $self->{completeness};
}

1;


