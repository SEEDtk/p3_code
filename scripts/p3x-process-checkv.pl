#!/usr/bin/env perl
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


use strict;
use warnings;
use FIG_Config;
use ScriptUtils;
use Stats;
use P3Utils;
use File::Copy::Recursive;
use CVUtils;

=head1 Process CheckV Output for Binning

    p3x-process-checkv.pl [ options ] checkv_db binDir outDir

This command processes the output from a CheckV job.  The output should be in the C<CheckV> sub-directory of the binning directory.  This
command will produce a tab-delimited summary of the viruses found and a C<vbin>I<NN>C<.fa> file for each virus bin.

The CheckV output directory will contain a C<completeness.csv> that describes each of the contigs that were processed by CheckV.  The
following columns are significant.

=over 4

=item contig_id

ID of the relevant contig

=item contig_length

length of the relevant contig

=item aai_completeness

percent of the virus covered by this contig

=item aai_error

error rating, used to filter poor matches

=item aai_top_hit

ID of the virus found

=back

The C<aai_top_hit> value can be found in the C<checkv_reps.tsv> file of the C<genome_db> sub-directory of the CheckV database directory.
In this file, the C<checkv_id> column contains a virus ID, and the C<type> column describes the file containing the virus information,
either C<checkv_circular.tsv> or C<checkv_genbank.tsv>.  It is worth knowing that taxonomic information is only available for the Genbank
viruses.

The following columns are significant in C<checkv_circular.tsv>.  Note that the taxonomic ID for all viruses in this file is the generic
superkingdom ID 10239.

=over 4

=item checkv_id

CheckV ID of the virus

=item original_id

name of the virus

=back

The following columns are significant in C<checkv_genbank.tsv>.

=over 4

=item checkv_id

CheckV ID of the virus

=item ncbi_id

taxonomic ID of the virus

=item ncbi_name

name of the virus

=back

=head2 Parameters

The positional parameters are the name of the CheckV database, the name of the binning directory, and the name of the output directory.
The command-line options are as follows.

=over 4

=item maxError

The maximum permissible error value for a contig to be considered valid.  The default is C<10>.

=item minPct

The minimum percentage of a virus that must be present in order to be binned.  The default is C<10>.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = P3Utils::script_opts('checkv_db binDir outDir',
        ['maxError|M=f', 'maximum permissible confidence error', { default => '10.0' }],
        ['minPct|p=f', 'minimum permissible percent of virus present', { default => '10.0' }]
        );
# Create a statistics object.
my $stats = Stats->new();
my $maxError = $opt->maxerror;
if ($maxError <= 0.0) {
    die "Maximum error must be > 0.";
}
my $minPct = $opt->minpct;
if ($minPct > 100.0) {
    die "Minimum percent must be <= 100.";
}
my ($checkVDb, $binDir, $outDir) = @ARGV;
if (! $checkVDb) {
    die "No CheckV database directory specified.";
} elsif (! -d $checkVDb) {
    die "CheckV database directory $checkVDb not found or invalid.";
}
if (! $binDir) {
    die "No binning directory specified.";
} elsif (! -d $binDir) {
    die "Binning directory $binDir not found.";
} elsif (! -s "$binDir/checkv/completeness.tsv") {
    die "Binning directory $binDir does not have CheckV output.";
}
if (! $outDir) {
    $outDir = $binDir;
    print "Output will be to binning directory $binDir.\n";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory $outDir: $!";
} else {
    print "Output will be to directory $outDir.\n";
}
# Create the processing object.
my $checkV = CVUtils->new($checkVDb, $binDir, $outDir, stats => $stats, maxE => $maxError, logH => \*STDOUT, minP => $minPct);
# Process the bins.
$checkV->CreateBins();
# All done.
print "All done.\n" . $stats->Show();


