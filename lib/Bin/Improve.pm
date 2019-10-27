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


package Bin::Improve;

    use strict;
    use warnings;
    use GEO;
    use FastA;
    use EvalCon;
    use Stats;
    use P3DataAPI;
    use RoleParse;

=head1 Use Quality Information to Produce an Improved Bin

This object operates on a L<Bin> object for which a FASTA file and a C<GenomeTypeObject> exists.  The latter must have come from
RASTtk, so that it contains quality information.  We will search for contigs with no good roles and remove them, creating a
new FASTA file.

This object contains the following fields.

=over 4

=item workDir

The name of the working directory containing the bins.

=item roleHashes

A 2-tuple containing (0) a map from each role ID to its name and (1) a map from each role checksum to a role ID.

=item p3

A L<P3DataAPI> object for accessing PATRIC.

=item stats

A L<stats> object for tracking statistics.

=item minComplete

The minimum completeness for a genome to be eligible for improvement.

=item variantMap

Reference to a hash keyed on subsystem name.  Each subsystem maps to a list of lists, each sub-list containing a list of roles forming a subsystem variant possibility.

=back

=head2 Special Methods

=head3 new

    my $improver = Bin::Improve->new($workDir, %options);

Create a new bin improvement object.

=over 4

=item workDir

The name of the working directory containing the binning files, or C<undef> if no binning files are present.
In this case, reference genomes are pulled from the PATRIC store.

=item options

A hash containing zero or more of the following options.

=over 8

=item p3

A L<P3DataAPI> object for accessing PATRIC.  If none is specified, one will be created.

=item stats

A L<Stats> object for tracking statistics.  If none is specified, one will be created.

=item minComplete

The minimum completeness for a genome to be eligible for improvement.  The default is C<90>.

=item roleFile

The C<roles.in.subsystems> file containing the mapping between role IDs, checksums, and role names.  The default is
C<roles.in.subsystems> in the P3 data directory.

=item variantMap

The C<variantMap.tbl> file containing the subsystem definitions. The default is C<variantMap.tbl> in the P3 data directory.

=back

=back

=cut

sub new {
    my ($class, $workDir, %options) = @_;
    # Get the helper objects.
    my $p3 = $options{p3} // P3DataAPI->new();
    my $stats = $options{stats} // Stats->new();
    my $roleFile = $options{roleFile} // "$FIG_Config::p3data/roles.in.subsystems";
    my $variantMap = $options{variantMap} // "$FIG_Config::p3data/variantMap.tbl";
    # Get the tuning options.
    my $min = $options{minComplete} // 90;
    # Get the role hashes.
    my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::p3data/roles.in.subsystems", $stats);
    # Get the variant map.
    my %vMap;
    open(my $vh, '<', $variantMap) || die "Could not open $variantMap: $!";
    while (! eof $vh) {
        my $line = <$vh>;
        chomp $line;
        my ($sub, $var, $roles) = split /\t/, $line;
        my @roles = split /\s+/, $roles;
        push @{$vMap{$sub}}, \@roles;
        $stats->Add(variantIn => 1);
    }
    # Create the object.
    my $retVal = { workDir => $workDir, roleHashes => [$nMap, $cMap], stats => $stats, p3 => $p3, minComplete => $min,
        variantMap => \%vMap };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Query Methods

=head3 eligible

    my $okFlag = $improver->eligible($gto);

Return TRUE if the specified genome is eligible for improvement, else FALSE.

=over 4

=item gto

A L<GenomeTypeObject> representing the genome to be checked.  It must include quality information.  (This will be the case if it
is output from RASTtk.)

=item RETURN

Return TRUE if the genome is eligible, else FALSE.

=back

=cut

sub eligible {
    my ($self, $gto) = @_;
    my $retVal;
    my $stats = $self->{stats};
    my $quality = $gto->{quality};
    if (! $quality) {
        $stats->Add(improveNoData => 1);
    } elsif ($quality->{completeness} >= $self->{minComplete} && ! GEO::contamX($quality->{contamination})) {
        $retVal = 1;
    } else {
        $stats->Add(improveNotEligible => 1);
    }
    return $retVal;
}


=head2 Public Manipulation Methods

=head3 Improve

    my $improvedFlag = $improver->Improve($refs, $gto);

Attempt to improve a bin.  The L<GenomeTypeObject> will be updated in place and its quality data deleted.

=over 4

=item refs

Reference to a list of IDs for the reference genomes.

=item gto

The L<GenomeTypeObject> containing the genome and its quality information.

=item RETURN

Returns TRUE if the bin was modified, else FALSE.

=back

=cut

sub Improve {
    my ($self, $refs, $gto) = @_;
    # Get the stats object.
    my $stats = $self->{stats};
    # Get the work directory.
    my $workDir = $self->{workDir};
    # This will be set to a TRUE if we improve the bin.
    my $retVal;
    # Create the GEO options.
    my %gOptions = (roleHashes => $self->{roleHashes}, detail => 2, p3 => $self->{p3}, stats => $stats);
    # Create the GEO for the sample bin.
    my $geo = GEO->CreateFromGto($gto, %gOptions);
    # Get the reference genome IDs.
    for my $refGenome (@$refs) {
        my $refGeo;
        if ($workDir && -s "$workDir/$refGenome.json") {
            $refGeo = GEO->CreateFromGto("$workDir/$refGenome.json", %gOptions);
        } else {
            my $gHash = GEO->CreateFromPatric($refGenome, %gOptions);
            $refGeo = $gHash->{$refGenome};
        }
        $geo->AddRefGenome($refGeo);
    }
    # Do a full quality analysis.
    $geo->AnalyzeQualityData();
    # Find the bad contigs.
    my $badHash = $geo->FindBadContigs();
    my $badFound = scalar keys %$badHash;
    if (! $badFound) {
        $stats->Add(improveNoBadContigs => 1);
    } else {
        # Remove the bad contigs.
        $self->TrimGto($gto, $badHash);
        $retVal = 1;
    }
    return $retVal;
}

=head3 TrimGto

    $improver->TrimGto($gto, $badH);

Remove bad contigs from a L<GenomeTypeObject> and adjust the features and subsystems accordingly.  The quality data will need to be regenerated.

=over 4

=item gto

The L<GenomeTypeObject> to be modified.

=item badH

Reference to a hash that maps the ID of each contig to be removed to a TRUE value.

=back

=cut

sub TrimGto {
    my ($self, $gto, $badH) = @_;
    my $stats = $self->{stats};
    # Physically remove the contigs.
    my $oldContigs = $gto->{contigs};
    my @newContigs;
    for my $contig (@$oldContigs) {
        if (! $badH->{$contig->{id}}) {
            push @newContigs, $contig;
            $stats->Add(contigKept => 1);
        }
    }
    $gto->{contigs} = \@newContigs;
    # This hash will track the removed features.
    my %lostFids;
    # Physically remove the features on the contigs.
    my $oldFids = $gto->{features};
    my @newFids;
    for my $fid (@$oldFids) {
        my $locList = $fid->{location};
        my $badContigFound;
        for my $loc (@$locList) {
            $badContigFound ||= $badH->{$loc->[0]};
        }
        if ($badContigFound) {
            $lostFids{$fid->{id}} = 1;
            $stats->Add(fidRemoved => 1);
        } else {
            push @newFids, $fid;
            $stats->Add(fidKept => 1);
        }
    }
    $gto->{features} = \@newFids;
    # Now physically remove the deleted features from the subsystems.  If the subsystem is no longer sufficiently
    # large, delete the subsystem too.
    my $variantMap = $self->{variantMap};
    my $oldSubs = $gto->{subsystems};
    my @newSubs;
    for my $sub (@$oldSubs) {
        my $bindings = $sub->{role_bindings};
        # We need not only the new role bindings, but the checksums of the roles kept.
        my (@newRoles, %checkSums, $changes);
        for my $binding (@$bindings) {
            my $oldFids = $binding->{features};
            my @newFids;
            for my $fid (@$oldFids) {
                if ($lostFids{$fid}) {
                    $stats->Add(subFidDeleted => 1);
                    $changes = 1;
                } else {
                    push @newFids, $fid;
                }
            }
            if (@newFids) {
                $binding->{features} = @newFids;
                push @newRoles, $binding;
                $checkSums{RoleParse::Checksum($binding->{role_id})} = 1;
            }
        }
        if (! $changes) {
            push @newSubs, $sub;
            $stats->Add(subsystemUnchanged => 1);
        } else {
            # Try to find a matching variant.
            my $subID = $sub->{name};
            $subID =~ tr/_/ /;
            my $variants = $variantMap->{$subID};
            if (! $variants) {
                print STDERR "Could not find $subID.\n"; ## TODO delete this
                $stats->Add(subsystemNotFound => 1);
            } else {
                my $found;
                for my $variant (@$variants) {
                    my $missing = 1;
                    for my $role (@$variant) {
                        $missing &&= $checkSums{$role};
                    }
                    $found ||= ! $missing;
                }
                if ($found) {
                    $stats->Add(subsystemUpdated => 1);
                    push @newSubs, $sub;
                } else {
                    $stats->Add(subsystemDeleted => 1);
                }
            }
        }
    }
    $gto->{subsystems} = \@newSubs;
    # Delete the old quality data.
    delete $gto->{quality};
}


1;


