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


package EvalHelper;

    use strict;
    use warnings;
    use EvalCon;
    use EvalCom::Tax;
    use EvalCom::Rep;
    use P3DataAPI;
    use P3Utils;
    use BinningReports;
    use GEO;
    use File::Temp;
    use File::Copy::Recursive;
    use GenomeTypeObject;
    use Data::Dumper;

=head1 Genome Evaluation Helper

This package is designed to help in evaluating genomes. It is not as efficient as the command-line scripts because it
only evaluates a single genome at a time; however, it is useful in web environments where a single genome evaluation is
all that is needed. The main L</Process> method takes as input a GTO file name or a PATRIC genome ID and an optional
reference genome ID and produces an output web page.  There is also L</ProcessGto> that evaluates a GTO already in
memory.

=head2 Special Methods

=head3 ProcessGto

    my $geo = EvalHelper::ProcessGto($gto, %options);

Evaluate a L<GenomeTypeObject> in memory. The C<quality> member of the GTO will be updated with the results of the operation.

=over 4

=item gto

A L<GenomeTypeObject> representing a genome to be evaluated.

=item options

A hash containing zero or more of the following options.

=over 8

=item ref

The PATRIC ID of a reference genome to use for comparison. If specified, the C<deep> option is implied.

=item deep

Compares the genome to a reference genome in order to provide more details on problematic roles. If this
option is specified and C<ref> is B<not> specified, a reference genome will be computed.

=item checkDir

The name of the directory containing the reference genome table and the completeness data files. The default
is C<CheckG> in the SEEDtk global data directory.

=item predictors

The name of the directory containing the role definition files and the function predictors for the consistency
checking. The default is C<FunctionPredictors> in the SEEDtk global data directory.

=item roleFile

The name of the C<roles.in.subsystems> file.  The default is to look in the C<FunctionPredictors> directory
and then the SEEDtk evaluation directory.

=item rolesToUse

The name of the C<roles.to.use> file. The default is to look in the C<FunctionPredictors> directory
and then the SEEDtk evaluation directory.

=item p3

A L<P3DataAPI> object for accessing the PATRIC database. If omitted, one will be created internally.

=item template

The name of the template file. The default is C<p3_code/lib/BinningReports/webdetails.tt> in the SEEDtk module directory.

=item outFile

The name of an optional output file. If specified, will contain the evaluation tool output.

=item outHtml

The name of an optional web page file. If specified, will contain the HTML output.

=item improve

The name of an optional FASTA file.  If specified, it will contain the contigs for a proposed improved version of the
genome.

=item external

If TRUE, the incoming genome is presumed to be external, and no contig links will be generated on the web page.

=item binned

If TRUE, the incoming genome is presumed to have user-specified contig IDs, which are stored as descriptions in the PATRIC genome_sequence records.
This overrides C<external>.

=item parallel

Number of parallel processes to run when evaluating the function predictors. The default is C<16>.

=item workDir

Name of a working directory for creating the matrix.  If none is specified, a temporary directory will be used.

=item RETURN

Returns the GEO used for evaluation.

=back

=back

=cut

sub ProcessGto {
    my ($gto, %options) = @_;
    # Connect to PATRIC.
    my $p3 = $options{p3} // P3DataAPI->new();
    # Create the work directory.
    my $workDir = $options{workDir};
    my $tmpObject;
    if (! $workDir) {
        $tmpObject = File::Temp->newdir();
        $workDir = $tmpObject->dirname;
    } elsif (! -d $workDir) {
        File::Copy::Recursive::pathmk($workDir) || die "Could not create work directory: $!";
    }
    # Get the role files.
    my $roleFile = $options{roleFile};
    my $rolesToUse = $options{rolesToUse};
    # Create the consistency helper.
    my $evalCon = EvalCon->new(predictors => $options{predictors}, roleFile => $roleFile, rolesToUse => $rolesToUse);
    # Get access to the statistics object.
    my $stats = $evalCon->stats;
    # Create the completeness helper.
    my $checkDir = $options{checkDir} // "$FIG_Config::p3data/CheckG";
    my ($nMap, $cMap) = $evalCon->roleHashes;
    my %evalOptions = (stats => $stats);
    my $evalG;
    if (-s "$checkDir/REP") {
        open(my $xh, '<', "$checkDir/REP") || die "Could not open REP file: $!";
        my $k = <$xh>;
        chomp $k;
        $evalG = EvalCom::Rep->new($checkDir, %evalOptions, K => $k);
    } else {
        $evalG = EvalCom::Tax->new($checkDir, %evalOptions, roleHashes=> [$nMap, $cMap]);
    }
    # Compute the detail level.
    my $detailLevel = (($options{deep} || $options{ref}) ? 2 : 1);
    # Set up the options for creating the GEOs.
    my %geoOptions = (roleHashes => [$nMap, $cMap], p3 => $p3, stats => $stats, detail => $detailLevel);
    # Start the predictor matrix for the consistency checker.
    $evalCon->OpenMatrix($workDir);
    my $retVal = GEO->CreateFromGto($gto, %geoOptions, external => $options{external}, binned => $options{binned});
    my $genomeID = $retVal->id;
    # Do we have a reference genome ID?
    my $refID = $options{ref};
    my %refMap;
    if (! $refID && $options{deep}) {
        # Here we must compute it, so we need to load the reference map.
        open(my $rh, "<$checkDir/ref.genomes.tbl") || die "Could not open reference genome table: $!";
        while (! eof $rh) {
            my $line = <$rh>;
            if ($line =~ /^(\d+)\t(\d+\.\d+)/) {
                $refMap{$1} = $2;
            }
        }
        # Get the lineage ID list.
        my $taxResults = [[$retVal->lineage || []]];
        $refID = _FindRef($taxResults, \%refMap, $genomeID);
    }
    if ($refID) {
        my $gHash = GEO->CreateFromPatric([$refID], %geoOptions);
        my $refGeo = $gHash->{$refID};
        if ($refGeo) {
            $retVal->AddRefGenome($refGeo);
        }
    }
    # Open the output file for the quality data.
    my $qFile = "$workDir/$genomeID.out";
    open(my $oh, '>', $qFile) || die "Could not open work file: $!";
    # Output the completeness data.
    my ($complete, $contam, $group, $seedFlag) = $evalG->Check2($retVal, $oh);
    close $oh;
    # Create the eval matrix for the consistency checker.
    $evalCon->OpenMatrix($workDir);
    $evalCon->AddGeoToMatrix($retVal);
    $evalCon->CloseMatrix();
    # Evaluate the consistency.
    my $par = $options{parallel} // 16;
    my $rc = system('eval_matrix', "-p", $par, '-q', $evalCon->predictors, $workDir, $workDir);
    if ($rc) {
        die "EvalCon returned error code $rc.";
    }
    # Store the quality metrics in the GEO.
    $retVal->AddQuality($qFile);
    # If there is an improvement file, create the improved FASTA.
    my $improveFile = $options{improve};
    if ($improveFile && $complete >= 90 && ! GEO::contamX($contam)) {
        # Find the bad contigs.
        my $badHash = $retVal->FindBadContigs();
        my $badFound = scalar keys %$badHash;
        if ($badFound) {
            # Write the good contigs.
            open(my $oh, '>', $improveFile) || die "Could not open FASTA file $improveFile: $!";
            my $contigs = $gto->{contigs};
            for my $contig (@$contigs) {
                my ($contigID) = $contig->{id};
                if (! $badHash->{$contigID}) {
                    print $oh ">$contigID\n$contig->{dna}\n";
                }
            }
        }
    }
    # If there is an output file, copy into it.
    if ($options{outFile}) {
        File::Copy::Recursive::fcopy($qFile, $options{outFile}) || die "Could not copy output file: $!";
    }
    # If there is an output web page file, we create the detail page.
    if ($options{outHtml}) {
        my $detailFile = $options{template} // "$FIG_Config::mod_base/p3_code/lib/BinningReports/webdetails.tt";
        my $page = BinningReports::Detail(undef, undef, $detailFile, $retVal, $nMap, undef, $options{outputGenomeBaseURL});
        open(my $oh, '>', $options{outHtml}) || die "Could not open HTML output file: $!";
        print $oh $page;
    }
    # Copy the quality metrics to the GTO.
    $retVal->UpdateGTO($gto);
    # Return the geo.
    return $retVal;
}


=head3 Process

    my $geo = EvalHelper::Process($genome, %options);

Evaluate a single genome.

=over 4

=item genome

The ID of a PATRIC genome or the file name of a L<GenomeTypeObject> for the genome to evaluate.

=item options

A hash containing zero or more of the following options.

=over 8

=item ref

The PATRIC ID of a reference genome to use for comparison. If specified, the C<deep> option is implied.

=item deep

Compares the genome to a reference genome in order to provide more details on problematic roles. If this
option is specified and C<ref> is B<not> specified, a reference genome will be computed.

=item checkDir

The name of the directory containing the reference genome table and the completeness data files. The default
is C<CheckR> in the SEEDtk evaluation directory.

=item predictors

The name of the directory containing the function predictors for the consistency
checking. The default is C<FunctionPredictors> in the SEEDtk evaluation directory.

=item roleFile

The name of the C<roles.in.subsystems> file.  The default is to look in the C<FunctionPredictors> directory
and then the SEEDtk evaluation directory.

=item rolesToUse

The name of the C<roles.to.use> file. The default is to look in the C<FunctionPredictors> directory
and then the SEEDtk evaluation directory.

=item p3

A L<P3DataAPI> object for accessing the PATRIC database. If omitted, one will be created internally.

=item template

The name of the template file. The default is C<p3_code/lib/BinningReports/webdetails.tt> in the SEEDtk module directory.

=item outFile

The name of an optional output file. If specified, will contain the evaluation tool output.

=item outHtml

The name of an optional web page file. If specified, will contain the HTML output.

=item external

If TRUE, the incoming genome is presumed to be external, and no contig links will be generated on the web page.

=item binned

If TRUE, the incoming genome is presumed to have user-specified contig IDs, which are stored as descriptions in the PATRIC genome_sequence records.
This overrides C<external>.

=back

=item RETURN

Returns a L<GEO> for the genome in question with embedded quality data.

=back

=cut

sub Process {
    my ($genome, %options) = @_;
    # Connect to PATRIC.
    my $p3 = $options{p3} // P3DataAPI->new();
    # Create the work directory.
    my $tmpObject = File::Temp->newdir();
    my $workDir = $tmpObject->dirname;
    # Get the role files.
    my $roleFile = $options{roleFile};
    my $rolesToUse = $options{rolesToUse};
    # Create the consistency helper.
    my $evalCon = EvalCon->new(predictors => $options{predictors}, roleFile => $roleFile, rolesToUse => $rolesToUse);
    # Get access to the statistics object.
    my $stats = $evalCon->stats;
    # Create the completeness helper.
    my $checkDir = $options{checkDir} // "$FIG_Config::p3data/Eval/CheckR";
    my ($nMap, $cMap) = $evalCon->roleHashes;
    my $evalG = EvalCom::Rep->new($checkDir, roleHashes=> [$nMap, $cMap], stats => $stats);
    # Compute the detail level.
    my $detailLevel = (($options{deep} || $options{ref}) ? 2 : 1);
    # Set up the options for creating the GEOs.
    my %geoOptions = (roleHashes => [$nMap, $cMap], p3 => $p3, stats => $stats, detail => $detailLevel);
    # Start the predictor matrix for the consistency checker.
    $evalCon->OpenMatrix($workDir);
    # This will be the two GEOs and actual genome ID (which is different if we have a GTO).
    my ($retVal, $refGeo, $genomeID);
    # Find out what kind of genome we have. If we have a GTO, we load it here. We wait on the PATRIC
    # load in case we need a reference genome ID.
    my @genomes;
    my $patricIn = 0;
    if ($genome =~ /^\d+\.\d+$/) {
        $patricIn = 1;
        $genomeID = $genome;
        push @genomes, $genomeID;
    } else {
        my $gHash = GEO->CreateFromGtoFiles([$genome], %geoOptions, external => $options{external}, binned => $options{binned});
        # A little fancy dancing is required because we don't know the genome ID, and it's the key to the hash we got
        # back. Thankfully, the hash is at most a singleton.
        ($genomeID) = keys %$gHash;
        if ($genomeID) {
            $retVal = $gHash->{$genomeID};
        } else {
            die "Could not load genome from $genome.";
        }
    }
    # Do we have a reference genome ID?
    my $refID = $options{ref};
    my %refMap;
    if (! $refID && $options{deep}) {
        # Here we must compute it, so we need to load the reference map.
        open(my $rh, "<$checkDir/ref.genomes.tbl") || die "Could not open reference genome table: $!";
        while (! eof $rh) {
            my $line = <$rh>;
            if ($line =~ /^(\d+)\t(\d+\.\d+)/) {
                $refMap{$1} = $2;
            }
        }
        # Get the lineage ID list.
        my $taxResults;
        if ($retVal) {
            $taxResults = [[$retVal->lineage || []]];
        } else {
            $taxResults = P3Utils::get_data_keyed($p3, genome => [], ['taxon_lineage_ids'], [$genome]);
        }
        $refID = _FindRef($taxResults, \%refMap, $genomeID);
        if ($refID) {
            push @genomes, $refID;
        }
    }
    # Do we have PATRIC genomes to read?
    if (@genomes) {
        my $gHash = GEO->CreateFromPatric(\@genomes, %geoOptions);
        if (! $retVal) {
            $retVal = $gHash->{$genome};
        }
        if ($refID) {
            $refGeo = $gHash->{$refID};
        }
    }
    # Now we have the two GEOs. Attach the ref to the main.
    if (! $retVal) {
        die "Could not read $genome from PATRIC.";
    } elsif ($refGeo) {
        $retVal->AddRefGenome($refGeo);
    }
    # Open the output file for the quality data.
    my $qFile = "$workDir/$genomeID.out";
    open(my $oh, '>', $qFile) || die "Could not open work file: $!";
    # Output the completeness data.
    $evalG->Check2($retVal, $oh);
    close $oh;
    # Create the eval matrix for the consistency checker.
    $evalCon->OpenMatrix($workDir);
    $evalCon->AddGeoToMatrix($retVal);
    $evalCon->CloseMatrix();
    # Evaluate the consistency.
    my $rc = system('eval_matrix', '-q', $evalCon->predictors, $workDir, $workDir);
    if ($rc) {
        die "EvalCon returned error code $rc.";
    }
    # Store the quality metrics in the GEO.
    $retVal->AddQuality($qFile);
    # If there is an output file, copy into it.
    if ($options{outFile}) {
        File::Copy::Recursive::fcopy($qFile, $options{outFile}) || die "Could not copy output file: $!";
    }
    # If there is an output web page file, we create the detail page.
    if ($options{outHtml}) {
        my $detailFile = $options{template} // "$FIG_Config::mod_base/p3_code/lib/BinningReports/webdetails.tt";
        my $retVal = BinningReports::Detail(undef, undef, $detailFile, $retVal, $nMap);
        open(my $oh, '>', $options{outHtml}) || die "Could not open HTML output file: $!";
        print $oh $retVal;
    }
    # Return the result.
    return $retVal;
}

=head2 Internal Utilities

=head3 _FindRef

    my $refID = EvalHelper::_FindRef($taxResults, $refMap, $genome);

Compute the reference genome ID from the results of a taxonomy lineage ID search.

=over 4

=item taxResults

The results of a query for the taxonomic lineage IDs of the current genome.

=item refMap

A reference to a hash mapping taxonomic IDs to reference genomes.

=item genome

The ID of the genome whose reference is desired.

=item RETURN

Returns the ID of the reference genome, or C<undef> if none could be found or the genome is its own reference.

=back

=cut

sub _FindRef {
    my ($taxResults, $refMap, $genome) = @_;
    # This will be the return value.
    my $retVal;
    # Insure we have a result.
    if ($taxResults && @$taxResults) {
        my $lineage = $taxResults->[0][0];
        # Loop through the lineage until we find something. Note that sometimes the lineage ID list comes back
        # as an empty string instead of a list so we need an extra IF.
        my $refFound;
        if ($lineage) {
            while (! $refFound && (my $tax = pop @$lineage)) {
                $refFound = $refMap->{$tax};
            }
        }
        # If we found something and it's not the genome of interest, run with it.
        if ($refFound && $refFound ne $genome) {
            $retVal = $refFound;
        }
    }
    return $retVal;
}


1;
