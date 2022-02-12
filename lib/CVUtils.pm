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


package CVUtils;

    use strict;
    use warnings;
    use P3Utils;
    use Stats;
    use VBin;
    use FastA;
    use VFile::Circular;
    use VFile::GenBank;
    use CGI;

=head1 CheckV Processing Utilities

This object contains useful methods for binning viruses using CheckV.  It contains the following fields.

=over 4

=item binDir

The name of the binning directory.

=item outDir

The name of the output directory.

=item dbDir

The name of the directory containing the CheckV database.

=item covgMap

Reference to a hash mapping each contig ID to its coverage.

=item stats

A L<Stats> object for tracking statistics about this run.

=item maxE

The maximum permissible error percentage for a contig to be binned.

=item minP

The minimum permissible percent of a virus that must be present to be considered a bin.

=item logH

An open file handle for output log messages.

=back

=head2 Special Methods

=head3 new

    my $checkV = CheckV->new($dbDir, $binDir, $outDir, %options);

Create a new utility object for binning viruses.

=over 4

=item dbDir

The name of the directory containing the CheckV database.

=item binDir

The name of the input binning directory.

=item outDir

The name of the directory to contain the output files.

=item options

A hash containing zero or more of the following options.

=over 8

=item stats

A L<Stats> object for keeping statistics.  The default is to create one internally.

=item maxE

The maximum permissible error percentage for a contig to be binned.  The default is C<10.0>.

=item minP

The minimum permissible percent of a virus that must be present to be considered a bin.  The default is C<10.0>.

=item logH

An open file handle for status messages.  The default is C<\*STDERR>.

=back

=back

=cut

sub new {
    my ($class, $dbDir, $binDir, $outDir, %options) = @_;
    # Get the options.
    my $stats = $options{stats} // Stats->new();
    my $maxE = $options{maxE} // 10.0;
    my $minP = $options{minP} // 10.0;
    my $logH = $options{logH} // \*STDERR;
    # Read the coverage data.
    my %covgMap;
    print $logH "Reading coverage data from $binDir.\n";
    open(my $ih, '<', "$binDir/output.contigs2reads.txt") || die "Could not open coverage file in $binDir: $!";
    my $line = <$ih>;
    while (! eof $ih) {
        $line = <$ih>;
        chomp $line;
        my ($contig, $covg) = split /\t/, $line;
        $covgMap{$contig} = $covg;
        $stats->Add(covgIn => 1);
    }
    # Create and return the object.
    my $retVal = {
        dbDir => $dbDir,
        binDir => $binDir,
        outDir => $outDir,
        stats => $stats,
        maxE => $maxE,
        minP => $minP,
        logH => $logH,
        covgMap => \%covgMap
    };
    bless $retVal, $class;
    return $retVal;
}

=head2 Processing Methods

=head3 CreateBins

    $checkV->CreateBins();

Create the bins and produce the output.

=cut

sub CreateBins {
    my ($self) = @_;
    my $logH = $self->{logH};
    # Separate the contigs into bins.
    my $binHash = $self->_computeBins();
    # Get the identifying information for each virus found.
    my $virHash = $self->_findViruses($binHash);
    # Build the bins.
    $self->_buildBins($binHash);
    # Write the report.  Note that we only output bins with a bin number.
    print $logH "Writing final report.\n";
    open(my $wh, '>', "$self->{outDir}/vbins.html") || die "Could not open virus bin web page: $!";
    open(my $oh, '>', "$self->{outDir}/vbins.tsv") || die "Could not open virus bin report: $!";
    my @cols = ("virus_id", "bin", "taxon_id", "name", "length", "completeness", "pct_error", "coverage");
    my @formats = ("text", "num", "text", "text", "num", "num", "num", "num");
    print $oh join("\t", @cols) . "\n";
    print $wh CGI::start_html(-title => "Virus Bin Summary");
    print $wh CGI::style("th, tr, td { border-style: inset; border-collapse: collapse; vertical-align: top; padding: 3px; }\nth { text-align: left; background: #EEEEEE; }\ntd.num, th.num { text-align: right; }");
    print $wh CGI::h1("Virus Bin Summary");
    print $wh CGI::p("This table lists all the known viruses found in the sample.");
    print $wh CGI::start_table();
    print $wh CGI::Tr( map { CGI::th({ class => $formats[$_] }, $cols[$_]) } 0..7);
    for my $virusID (sort { _vcmp($binHash, $a, $b) } keys %$binHash) {
        my $vBin = $binHash->{$virusID};
        if ($vBin->num > 0) {
            my ($taxonID, $name) = @{$virHash->{$virusID}};
            $name //= "Unknown virus $virusID";
            if ($taxonID) {
                @cols = ($virusID, $vBin->num, $taxonID, $name, $vBin->len, sprintf("%6.2f", $vBin->percent),
                        sprintf("%6.2f", $vBin->err), sprintf("%6.2f", $vBin->covg));
                print $oh join("\t", @cols) . "\n";
                print $wh CGI::Tr( map { CGI::td({ class => $formats[$_] }, $cols[$_]) } 0..7);
            }
        }
    }
    print $wh CGI::end_table();
    print $wh CGI::end_html();
    close $oh;
}


=head2 Internal Methods

=head3 _vcmp

    my $cmp = _vcmp($binHash, $a, $b);

Compare two bins to sort them for the final report.

=over 4

=item binHash

Reference to the master bin hash.

=item a

First virus bin ID.

=item b

Second virus bin ID.

=item RETURN

Returns a negative value if B<a> goes first, 0 if the bins are equal, and a positive value otherwise.

=back

=cut

sub _vcmp {
    my ($binHash, $a, $b) = @_;
    my ($aBin, $bBin) = map { $binHash->{$_} } ($a, $b);
    my ($aType, $bType) = ($aBin->{type}, $bBin->{type});
    my $retVal = 0;
    if ($aType eq "genbank" && $bType ne "genbank") {
        $retVal = -1;
    } elsif ($aType ne "genbank" && $bType eq "genbank") {
        $retVal = 1;
    } else {
        $retVal = (abs(100 - $aBin->percent) <=> abs(100 - $bBin->percent));
    }
    return $retVal;


}

=head3 _computeBins

    my $binHash = $checkV->_computeBins();

Process the C<completeness.tsv> file and return a reference to a hash mapping each virus ID to a L<VBin> object for that virus.

=cut

sub _computeBins {
    my ($self) = @_;
    my $stats = $self->{stats};
    my $logH = $self->{logH};
    my $binDir = $self->{binDir};
    my $maxE = $self->{maxE};
    my $covgH = $self->{covgMap};
    # This will be the return hash.
    my %retVal;
    # Open the completeness.tsv file.
    open(my $ih, '<', "$binDir/checkv/completeness.tsv") || die "Could not open virus completeness file: $!";
    my (undef, $cols) = P3Utils::find_headers($ih, completeness => 'contig_id', 'contig_length', 'aai_completeness', 'aai_error', 'aai_top_hit');
    # Read all the contig results and keep the good ones.
    my ($contigCount, $binCount) = (0, 0);
    while (! eof $ih) {
        my ($contigID, $len, $percent, $err, $virusID) = P3Utils::get_cols($ih, $cols);
        $stats->Add(contigIn => 1);
        # Skip a non-viral contig.
        if ($virusID eq 'NA') {
            $stats->Add(nonViral => 1);
        } elsif ($err eq 'NA') {
            $stats->Add(badViral => 1);
        } else {
            if ($err > $maxE) {
                $stats->Add(badHit => 1);
            } else {
                # Get the bin for this virus.
                my $vBin = $retVal{$virusID};
                if (! $vBin) {
                    $vBin = VBin->new($virusID);
                    $retVal{$virusID} = $vBin;
                    $stats->Add(virusBin => 1);
                    $binCount++;
                } else {
                    $stats->Add(virusDup => 1);
                }
                # Compute the coverage and add the contig to the bin.
                my $covg = $covgH->{$contigID} // 1.0;
                $vBin->AddContig($contigID, $len, $percent, $err, $covg);
            }
        }
        $contigCount++;
        if ($contigCount % 500 == 0) {
            print $logH "$contigCount contigs processed. $binCount bins found.\n";
        }
    }
    # Return the bin map.
    return \%retVal;
}

=head3 _findViruses

    my $virHash = $checkV->_findViruses(\%binHash);

Create a hash mapping each virus ID to its name and taxonomic ID.

=over 4

=item binHash

Hash of virus IDs to bin structures.

=item RETURN

Returns a reference to a hash mapping each virus ID to a [taxomomyID, name] pair.

=back

=cut

sub _findViruses {
    my ($self, $binHash) = @_;
    my $stats = $self->{stats};
    my $logH = $self->{logH};
    # Get the handlers for the two types of viruses.
    my $dbDir = $self->{dbDir};
    my %vHash = (circular => VFile::Circular->new($dbDir),
            genbank => VFile::GenBank->new($dbDir));
    # Get a hash for the viruses we're trying to find.
    my %idHash = map { $_ => 1 } keys %$binHash;
    # This hash will accumulate the virus IDs of each type.
    my %typeHash;
    # Open the reference file and loop through it.
    print $logH "Processing virus ID list.\n";
    open(my $ih, '<', "$dbDir/genome_db/checkv_reps.tsv") || die "Could not open virus reps file: $!";
    my $line = <$ih>;
    while (! eof $ih) {
        my ($virusID, $type) = P3Utils::get_fields($ih);
        if ($idHash{$virusID}) {
            push @{$typeHash{$type}}, $virusID;
            $binHash->{$virusID}{type} = $type;
            $stats->Add("$type-bin" => 1);
        }
    }
    # Process each virus type.
    my %retVal;
    for my $type (keys %typeHash) {
        print "Identifying $type viruses.\n";
        if (! exists $vHash{$type}) {
            print $logH "WARNING: Unknown virus type $type.\n";
        } else {
            $vHash{$type}->Process($typeHash{$type}, \%retVal);
        }
    }
    return \%retVal;
}

=head3 _buildBins

    $checkV->_buildBins($binHash);

Create the FASTA files for the virus bins.

=over 4

=item binHash

Reference to a hash mapping each virus ID to a L<VBin> object.

=back

=cut

sub _buildBins {
    my ($self, $binHash) = @_;
    my $stats = $self->{stats};
    my $outDir = $self->{outDir};
    my $logH = $self->{logH};
    my $minP = $self->{minP};
    # This will track the bin number.
    my $num = 1;
    # Load the contigs into memory.
    my $contigs = $self->_loadContigs($binHash);
    # Loop through the bins.
    my @viruses = sort { $binHash->{$b}->percent <=> $binHash->{$a}->percent } grep { $binHash->{$_}->percent >= $minP } keys %$binHash;
    for my $virusID (@viruses) {
        my $vBin = $binHash->{$virusID};
        $vBin->num($num);
        print $logH "Writing bin $num for $virusID.\n";
        open(my $oh, '>', "$outDir/vBin$num.fa") || die "Could not open output file for virus bin $virusID: $!";
        for my $contigID (@{$vBin->contigs}) {
            print $oh ">$contigID\n$contigs->{$contigID}\n";
            $stats->Add(contigOut => 1);
        }
        $stats->Add(binOut => 1);
        close $oh;
        $num++;
    }
}

=head3 _loadContigs

    my $contigHash = $checkV->_loadContigs($binHash);

Create a hash of all the contigs being put into bins.

=over 4

=item binHash

Reference to a hash of virus IDs to L<VBin> objects.

=item RETURN

Returns a hash mapping contig IDs to DNA sequences.

=back

=cut

sub _loadContigs {
    my ($self, $binHash) = @_;
    my $stats = $self->{stats};
    my $logH = $self->{logH};
    # Create a hash of the contigs we want.
    my %retVal;
    for my $virusID (keys %$binHash) {
        my $contigIDs = $binHash->{$virusID}->contigs;
        for my $contigID (@$contigIDs) {
            $retVal{$contigID} = "";
            $stats->Add(contigSet => 1);
        }
    }
    # Now load the actual contigs.
    print $logH "Loading contigs from unbinned.fasta.\n";
    my $fh = FastA->new("$self->{binDir}/unbinned.fasta");
    while ($fh->next) {
        my $contigID = $fh->id;
        if (exists $retVal{$contigID}) {
            $retVal{$contigID} = $fh->left;
            $stats->Add(contigFound => 1);
        }
    }
    return \%retVal;
}

1;


