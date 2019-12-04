#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
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

package EvalCom::Rep;

    use strict;
    use warnings;
    use RepGenome;
    use ScriptUtils;
    use base qw(EvalCom);

=head1 Evaluate Completeness Based on Representative Genomes

This is a subclass of L<EvalCom> that classifies genomes according to their closest representative genome rather than their taxonomic
grouping.  The theory is that this will reduce errors due to incorrect taxonomic classification and genome distribution bias.

As a fallback, there is a universal role table for all genomes.

The input files are found in the working directory (usually C<CheckR> in the SEEDtk global directory).  C<comp.tbl> contains
the universal roles for each representative genome group, including the backup.  For each group, the header line will contain
the group ID, the minimum similarity score for group membership, the group name, and the seed protein sequence.  Each data
line will contain a role ID.  The trailer line will be two slashes (C<//>).

The goal of this object is to, given a GEO, find the appropriate universal-role hash to use for computing completeness and
contamination.  A universal-role hash maps role IDs to weights (see L<EvalCom>).

The object contains the following fields.

=over 4

=item repObjects

A hash mapping representative-genome group IDs to L<RepGenome> objects.  The genome name will be the group name.

=item repRoles

A hash mapping representative-genome group IDs to universal-role hashes.

=item repScores

A hash mapping representative-genome group IDs to minimum similarity scores.

=item repList

A sorted list of representative groups, with the highest-scoring ones first.

=back

=head2 Special Methods

=head3 new

    my $checker = EvalCom::Rep->new($inputDir, %options);

Create a new consistency/completeness checker from the specified input directory.

=over 4

=item inputDir

The input directory containing the C<comp.tbl> file.

=item options

A hash containing zero or more of the following keys.

=over 8

=item K

The kmer length for the similarity checks.  The default is C<8>.

=item logH

Open handle for an output stream to contain log messages. The default is to not write log messages.

=item stats

A L<Stats> object for tracking statistics. If none is specified, one will be created.

=back

=back

=cut

sub new {
    my ($class, $inputDir, %options) = @_;
    # Create and bless the object.
    my $retVal = EvalCom::new($class, %options);
    # Get the statistics object.
    my $stats = $retVal->stats;
    # Get the kmer size.
    my $K = $options{K} // 8;
    $retVal->{K} = $K;
    # Read the representative-genome hashes.
    my (%repObjects, %repRoles, %repScores);
    open(my $ih, '<', "$inputDir/comp.tbl") || die "Could not open comp.tbl: $!";
    my $groupCount = 0;
    while (! eof $ih) {
        # Read the header line.
        my ($id, $score, $name, $prot) = ScriptUtils::get_line($ih);
        $repObjects{$id} = RepGenome->new($id, name => $name, prot => $prot);
        $repScores{$id} = $score;
        $stats->Add("group${score}In" => 1);
        # Read the universal role hash.
        my $line = <$ih>;
        my %weights;
        while (substr($line, 0, 2) ne '//') {
            $line =~ s/\s+//g;
            $weights{$line} = 1.0;
            $line = <$ih>;
        }
        $repRoles{$id} = \%weights;
        $groupCount++;
    }
    $retVal->Log("$groupCount representative-genome groups read.\n");
    close $ih; undef $ih;
    $retVal->{repObjects} = \%repObjects;
    $retVal->{repScores} = \%repScores;
    $retVal->{repRoles} = \%repRoles;
    # Now we sort the groups by score limit.
    my @repList = sort { ($repScores{$b} <=> $repScores{$a}) && $a cmp $b } keys %repRoles;
    $retVal->{repList} = \@repList;
    # Return the object.
    return $retVal;
}

=head2 Virtual Methods

=head3 Choose

    my ($group, $roleH) = $checker->Choose($geo);

Select the evaluation group for this genome.  The group name and the universal role hash are returned.

This method MUST be overridden by the subclass.

=over 4

=item geo

The L<GEO> for the genome to evaluate.

=item RETURN

Returns a two-element list consisting of (0) the name of the evaluation group and (1) a hash mapping each universal role ID to its weight.

=back

=cut

sub Choose {
    my ($self, $geo) = @_;
    # These will be the return values.
    my $repGroup = '';
    my $roleHash;
    # Get the statistics object.
    my $stats = $self->{stats};
    # Now we do the representative-genome search.  Everything matches the fall-back group, which is last in the list.
    # Get the seed protein.
    my $prot = $geo->seed;
    # Get the hashes.
    my $repScores = $self->{repScores};
    my $repObjects = $self->{repObjects};
    my $repRoles = $self->{repRoles};
    my $repList = $self->{repList};
    # Find the best-scoring group.
    my $bestScore = 0;
    for my $group (@$repList) {
        my $score = $repObjects->{$group}->check_genome($prot);
        if ($score > $bestScore && $score >= $repScores->{$group}) {
            $repGroup = $repObjects->{$group}->name;
            $roleHash = $repRoles->{$group};
            $bestScore = $score;
        }
    }
    # Fall back to the last group.
    if ($bestScore == 0) {
        $repGroup = $repObjects->{0}->name;
        $roleHash = $repRoles->{0};
    }
    return ($repGroup, $roleHash);
}

1;