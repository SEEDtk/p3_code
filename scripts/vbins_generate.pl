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

=head1 Generate Virus Bins

    vbins_generate.pl [options] checkvDB binDir

This script will generate the virus-based bins for a binning directory.  It will invoke B<checkv> on the C<unbinned.fasta> file and
then use L<CVUtils> to analyze the results and produce virus bins.

=head2 Parameters

The positional parameters are the name of the checkV database directory and the name of the binning directory.

The command-line options are as follows.

=over 4

=item maxError

The maximum permissible error value for a contig to be considered valid.  The default is C<10>.

=item minPct

The minimum percentage of a virus that must be present in order to be binned.  The default is C<10>.

=item threads

The number of compute threads to used, passed down to checkv. The default is C<1>.

=back

=cut

use strict;
use P3Utils;
use Stats;
use P3Utils;
use File::Copy::Recursive;
use CVUtils;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('checkvDB binDir',
			       ['maxError|M=f', 'maximum permissible confidence error', { default => '10.0' }],
			       ['minPct|p=f', 'minimum permissible percent of virus present', { default => '10.0' }],
			       ['threads|t=i', 'number of threads for checkv', { default => 1 }],
    );
my $stats = Stats->new();
my $maxError = $opt->maxerror;
if ($maxError <= 0.0) {
    die "Maximum error must be > 0.";
}
my $minPct = $opt->minpct;
if ($minPct > 100.0) {
    die "Minimum percent must be <= 100.";
}
my ($checkVDb, $binDir) = @ARGV;
if (! $checkVDb) {
    die "No CheckV database directory specified.";
} elsif (! -d $checkVDb) {
    die "CheckV database directory $checkVDb not found or invalid.";
}
if (! $binDir) {
    die "No binning directory specified.";
} elsif (! -d $binDir) {
    die "Binning directory $binDir not found.";
} elsif (! -s "$binDir/unbinned.fasta") {
    die "Binning directory $binDir does not have unbinned contig FASTA file.";
}
# Invoke checkv to identify viral contigs.
my @parms = ('checkv', 'completeness', "-t", $opt->threads, "$binDir/unbinned.fasta", "$binDir/checkv", '-d', $checkVDb);
my $cmd = join(" ", @parms);
print "Running: $cmd\n";
my $rc = system($cmd);
if ($rc) {
    die "Error in CheckV (rc = $rc).\n";
}
if ($rc == 0) {
    # Now we must process the checkv output.
    print "Processing checkV output for $binDir.\n";
    $stats->Add(outputRuns => 1);
    my $checkv = CVUtils->new($checkVDb, $binDir, $binDir, logH => \*STDOUT, stats => $stats, maxE => $maxError, minP => $minPct);
    $checkv->CreateBins();
}
print "All done:\n" . $stats->Show();

