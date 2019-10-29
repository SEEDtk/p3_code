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

=head1 Improve the Quality of a GTO file

    p3x-improve-gto.pl [options] refGenome gtoIn gtoOut

This script removes bad contigs from a L<GenomeTypeObject> to improve its quality.  Any quality data
present in the GTO will be deleted in the new one and must be recomputed.

The standard output will contain a single record with the value C<1> if the GTO was modified, and C<0>
otherwise.

=head2 Parameters

The positional parameters are the ID of the reference genome, the file name of the input L<GenomeTypeObject>,
and the output filename.  If the output filename is omitted, the input file will be overwritten with the
new GTO.  The input GTO must have quality data in it. (This is always the case if it is the output of
the annotation service.)  If the reference genome ID is C<x>, then it is computed internally.

Additional command-line options are as follows.

=over 4

=item workDir

The name of the working directory.  This is used to store temporary files.  If a L<GenomeTypeObject> file
for the reference genome already exists, it should be in this directory with the name I<refGenome>C<.json>.
(This will be the case if the incoming GTO is the output of a binning job in that directory.) If no such
file is found, the reference genome GTO will be downloaded from PATRIC using the web API.  The default is
the current directory.

=item minComplete

The minimum percent completeness for a genome to be eligibile for improvement.  The default is C<90>.

=item roleFile

The C<roles.in.subsystems> file containing the mapping between role IDs, checksums, and role names.  The
default is C<role.in.subsystems> in the P3 data directory.

=item variantMap

The C<variantMap.tbl> file containing the subsystem definitions.  The default is C<variantMap.tbl> in the
P3 data directory.

=item refGenomes

The C<ref.genomes.tbl> file used to look up reference genomes.  The default is C<CheckG/ref.genomes.tbl> in the
PATRIC data directory.  The file is only loaded if no reference genome is provided on the command line.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Bin::Improve;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('refGenome gtoIn gtoOut',
    ['workDir|work|w=s', 'working directory', { default => '.' }],
    ['minComplete|min|m=i', 'minimum acceptable completeness', { default => 90 }],
    ['roleFile|roles|R=s', 'role definition file', { default => "$FIG_Config::p3data/roles.in.subsystems" }],
    ['variantMap|subs|S=s', 'subsystem definition file', { default => "$FIG_Config::p3data/variantMap.tbl" }],
    ['refGenomes|G=s', 'reference-genome lookup file', { default => "$FIG_Config::p3data/CheckG/ref.genomes.tbl" }],
    );
my $stats = Stats->new();
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Get the options.
my $workDir = $opt->workdir;
my $minComplete = $opt->mincomplete;
my $roleFile = $opt->rolefile;
my $variantMap = $opt->variantmap;
my $refGenomes = $opt->refgenomes;
# Get the positional parameters.
my ($refGenome, $gtoIn, $gtoOut) = @ARGV;
if (! $refGenome) {
    die "No reference genome specified.";
} elsif (! $gtoIn) {
    die "Missing input GTO.";
} elsif (! -s $gtoIn) {
    die "Invalid or missing input GTO file $gtoIn.";
} elsif (! $gtoOut) {
    $gtoOut = $gtoIn;
}
if ($refGenome =~ /^\d+\.\d+$/) {
    # Here we know the reference genome ID.
    $refGenomes = undef;
}
# Create the improver.
print STDERR "Working directory is $workDir with role file $roleFile.\n";
my $improver = Bin::Improve->new($workDir, p3 => $p3, stats => $stats, minComplete => $minComplete,
        roleFile => $roleFile, variantMap => $variantMap, refGenomes => $refGenomes);
# This will be the output value.
my $retVal = '0';
# Read the input GTO.
print STDERR "Reading genome from $gtoIn.\n";
my $gto = GenomeTypeObject->create_from_file($gtoIn);
if (! $improver->eligible($gto)) {
    print STDERR "Genome not eligible for improvement.\n";
} else {
    # Get the reference genome ID.
    if ($refGenome eq 'x') {
        $refGenome = $improver->find_ref($gto);
        die "Could not compute reference genome." if ! $refGenome;
        print STDERR "Reference genome is $refGenome.\n";
    }
    print STDERR "Processing genome for improvement.\n";
    my $ok = $improver->Improve([$refGenome], $gto);
    if (! $ok) {
        print STDERR "Genome could not be improved.\n";
    } else {
        print STDERR "Writing improved genome to $gtoOut.\n";
        $gto->destroy_to_file($gtoOut);
        $retVal = '1';
    }
    print "All done.\n" . $stats->Show();
}
print "$retVal\n";
