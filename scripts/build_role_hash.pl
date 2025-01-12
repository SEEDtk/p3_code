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

=head1 Create a File of Hash Codes for Each Subsystem Role

    build_role_hash.pl [options] < variant.tbl > roleHash.tbl

This script builds an important intermediate file for the L<SubsystemProjector>. It reads the B<variant.tbl> file that
is output by the CoreSEED-processing utilities and outputs a file that maps each role name to the hash code used by
the projector.

=head2 Parameters

There are no positional parameters. The standard input should be the B<variant.tbl> file and the standard output
will be the role hash.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

=back

=cut

use strict;
use P3Utils;
use RoleParse;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::ih_options());

# Open the input file.
my $ih = P3Utils::ih($opt);
# Write the output header.
print "role_id\trole_name\tchecksum\n";
# Loop through the input until we hit the "**" marker.
print STDERR "Scanning input.\n";
my $done;
my $count = 0;
while (! eof $ih && ! $done) {
    my $line = <$ih>;
    chomp $line;
    if ($line eq "**") {
        $done = 1;
    } else {
        # Here we have a data line. The first field is the role ID and the second is the role name.
        my ($id, $name) = split /\t/, $line;
        my $checksum = RoleParse::Checksum($name);
        print "$id\t$name\t$checksum\n";
        $count++;
    }
}
print STDERR "$count lines processed.\n";
