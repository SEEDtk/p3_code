=head1 Evaluate a GenomeTypeObject

    p3x-eval-gto.pl [options] gtoFile outFile htmlFile

This script evaluates a single L<GenomeTypeObject> and stores the quality data back into it.

=head2 Parameters

The positional parameters are the file name of the input L<GenomeTypeObject>, the name of an output file for status information,
and the name of a web page file for the quality analysis.

Additional command-line options are as follows:

=over 4

=item ref

The PATRIC ID of a reference genome to use for comparison.  If specified, the C<deep> option is implied.

=item deep

If specified, the genome is compared to a reference genome in order to provide more details on problematic roles.
If this option is specified and C<ref> is not specified, a reference genome will be computed.

=item evalDir

The name of the directory containing the evaluation files. The default is C<Eval> in the global data directory.

=item workDir

Name of a working directory for intermediate files.  If none is specified, a temporary directory will be used.

=item test

If specified, the new GTO will be left in the work directory with the name "eval.gto", and the original GTO will
not be updated.

=back

=cut

use strict;
use P3Utils;
use P3DataAPI;
use File::Temp;
use File::Copy::Recursive;

# Get the command-line options.
my $opt = P3Utils::script_opts('gtoFile outFile outHtml',
        ['ref|r=s', 'reference genome ID (implies deep)'],
        ['deep', 'if specified, the genome is compared to a reference genome for more detailed analysis'],
        ['evalDir|eval=s', 'evaluation data directory', { default => "$FIG_Config::p3data/Eval" }],
        ['checkDir=s', 'no longer used'],
        ['predictors=s', 'no longer used'],
        ['template=s', 'no longer used'],
        ['external', 'no longer used'],
        ['binned', 'no longer used'],
        ['improve=s', 'no longer used'],
        ['parallel=i', 'no longer used'],
        ['workDir=s', 'name of a working directory for the java output (erased after use)'],
        ['genomeBaseUrl=s', 'no longer used'],
        ['test', 'test evaluation without updated GTO']
        );
# Get the input parameters.
my ($gtoFile, $outFile, $outHtml) = @ARGV;
if (! $gtoFile) {
    die "No input GTO file specified.";
} elsif (! $outFile) {
    die "No output data file specified.";
} elsif (($opt->ref || $opt->deep) && ! $outHtml) {
    die "No output web page file specified.";
}
# Verify the GTO.
if (! -s $gtoFile) {
    die "Input GTO file is empty or not found.";
}
# Find the work directory.
my $workDir = $opt->workdir;
my $tmpObject;
if (! $workDir) {
    $tmpObject = File::Temp->newdir;
    $workDir = $tmpObject->dirname;
} elsif (! -d $workDir) {
    File::Copy::Recursive::pathmk($workDir) || die "Could not create work directory: $!";
}
my $jarDir = $ENV{SEED_JARS};
if (! $jarDir) {
    die "SEED_JARS environment variable must be specified!";
}
if (! $ENV{P3API_URL}) {
    # We need the P3 data API URL in the environment.
    my $p3 = P3DataAPI->new();
    $ENV{P3API_URL} = $p3->{url};
}
# Format the parameters.
my $bufferFile = "$workDir/eval.gto";
my $reportFile = "$workDir/GenomeReport.html";
my @parms = ("gto", "--home", "BV-BRC", "--outDir", $workDir, "--output", $bufferFile, "--input", $gtoFile);
if ($opt->ref) {
    push @parms, "--ref", $opt->ref;
}
push @parms, "--format";
if ($opt->deep) {
    push @parms, "DEEP";
} else {
    push @parms, "HTML";
}
if ($opt->test) {
    push @parms, "--verbose";
}
push @parms, $opt->evaldir;
# Clear the log file.
if (-s $outFile) {
    unlink $outFile;
}
# Call the evaluator.
my $rc = system(java => "-Dlogback.configurationFile=$jarDir/dl4j.eval.logback.xml", "-Dlogfile.name=$outFile", "-jar",
        "$jarDir/dl4j.eval.jar", @parms);
if ($rc) {
    die "Error in evaluation: rc = $rc.";
}
# Copy results from work directory.
File::Copy::Recursive::fcopy($reportFile, $outHtml) || die "Could not copy html file: $!";
if (! $opt->test) {
    File::Copy::Recursive::fcopy($bufferFile, $gtoFile) || die "Could not copy updated GTO: $!";
}

