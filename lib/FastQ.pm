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


package FastQ;

    use strict;
    use warnings;

=head1 FASTQ Reader

This package facilitates reading paired-end FASTQ files. Such files are presented as a matched set or are interlaced.
To begin using a FASTQ file, simply construct this object with a single file name (interlaced) or a pair (matched set). You can then input
pairs of reads in a single operation.

This object contains the following fields.

=over 4

=item lh

Open file handle for the left reads, or for both if we are interlaced.

=item rh

Open file handle fo the right reads, or C<undef> if we are interlaced.

=item left

Left DNA string.

=item right

Right DNA string.

=item lqual

Left quality string.

=item rqual

Right quality string.

=item id

ID of the current node.

=item singleton

If TRUE, this is a singleton file, so there are no right reads.

=item q_base

C<33> for normal fastQ, C<64> for old Illumina format.  This is the base value for the quality strings.

=item unsafe

If defined, reference to a list of saved records resulting from unpaired data.

=back

=head2 Special Methods

=head3 new

    my $fqhandle = FastQ->new($left, $right, $options);

or

    my $fqhandle = FastQ->new($interlaced, $options);

Construct a new FASTQ handler using a matched pair of files or a single interlaced file. The handler may be used to retrieve matched pairs
of reads from the file.

=over 4

=item left

Name of the file containing the left-end reads, or an open file handle for it.

=item right

Name of the file containing the right-end reads, or an open file handle for it.

=item interlaced

Name of the file containing the interlaced reads, or an open file handle for it.

=item options (optional)

If specified, a reference to a hash containing zero or more of the following options.

=over 8

=item old_illumina

If specified, the quality strings use old Illumina scoring.

=item singleton

If specified, the incoming file is a singleton file with no right strings.

=item unsafe

If specified, errors in the incoming file will be skipped.

=back

=cut

sub new {
    my ($class, $left, $right, $options) = @_;
    # Adjust for the possibility of options.
    if (ref $right eq 'HASH') {
        $options = $right;
        $right = undef;
    } elsif (! $options) {
        $options = {};
    }
    # If we are allowing unsafe reading, we need to set up a buffering queue.
    my $unsafe;
    if ($options->{unsafe} && ! $right) {
        $unsafe = [];
    }
    # This will be the new object. It starts blank.
    my $retVal = {
        left => '',  right => '',
        lqual => '', rqual => '',
        id => undef,
        unsafe => $unsafe,
        singleton => $options->{singleton},
        q_base => ($options->{old_illumina} ? 64 : 33)
    };
    # Store the handle for the left file.
    my $lh;
    if (ref($left) =~ /^(?:GLOB|IO::)/) {
        $lh = $left;
    } else {
        open($lh, "<$left") || die "Could not open FASTQ file $left: $!";
    }
    $retVal->{lh} = $lh;
    # Store the handle for the right file if we are not interlaced or singleton.
    if ($right) {
        my $rh;
        if (ref($right) =~ /^(?:GLOB|IO::)/) {
            $rh = $right;
        } else {
            open($rh, "<$right") || die "Could not open FASTQ file $right: $!";
        }
        $retVal->{rh} = $rh;
    }
    # Bless and return this object.
    bless $retVal, $class;
    return $retVal;
}

=head3 norm_id

    my $normalized = FastQ::norm_id($id);

Strip off the direction indicator from a FASTQ file ID.

=over 4

=item id

The ID to normalize.

=item RETURN

Returns the normalized sequence ID.

=back

=cut

sub norm_id {
    my ($id) = @_;
    my $retVal;
    if ($id =~ /(.+)[\/\.][12]/) {
        $retVal = $1;
    } else {
        $retVal = $id;
    }
    return $retVal;
}

=head2 OrganizeFiles

    my $filesL = FastQ::OrganizeFiles($iFlag, @fileParms);

Organize the list of files containing FastQ data into a list of FastQ object parameter lists. If the files are interlaced, this
means creating a list of singeltons; otherwise, it means creating a list of pairs. The output list can be processed sequentially
to create a sequence of FastQ objects for input.

=over 4

=item iFlag

TRUE if the input is interlaced, FALSE if it is paired.

=item fileParms

The list of files to process for input.

=item RETURN

Returns a reference to a list of specifications for FastQ constructors.

=back

=cut

sub OrganizeFiles {
    my ($iFlag, @fileParms) = @_;
    my @retVal;
    if ($iFlag) {
        @retVal = map { [$_] } @fileParms;
    } else {
        my $n = scalar @fileParms;
        if ($n & 1) {
            die "Odd number of files specified in paired mode.";
        }
        for (my $i = 0; $i < $n; $i += 2) {
            push @retVal, [$fileParms[$i], $fileParms[$i+1]];
        }
    }
    return \@retVal;
}


=head2 Public Manipulation Methods

=head3 next

    my $found = $fqhandle->next;

Move forward to the next record, returning TRUE if one was found.

=cut

sub next {
    my ($self) = @_;
    # This will be set to TRUE if everything works.
    my $retVal;
    # Get the file handles.
    my ($lh, $rh) = ($self->{lh}, $self->{rh});
    # This will hold the left ID.
    my $leftID = '';
    # Check for buffering.
    my $bufferedL = $self->{unsafe};
    if ($bufferedL && scalar @$bufferedL) {
        # We have a buffered record.
        ($self->{id}, $self->{left}, $self->{lqual}) = @$bufferedL;
        $leftID = $self->{id};
        $self->{unsafe} = [];
    } else {
        # Check for end-of-file.
        if (! eof $lh) {
            # Read the left record.
            $leftID = $self->_read_fastq($lh, 'left');
        }
    }
    if ($self->{singleton}) {
        $retVal = ($leftID ? 1 : 0);
    } else {
        # Determine from where we will get the right record. If there is no right
        # file, it will be the left file (interlaced mode).
        my $rightID = '';
        $rh //= $lh;
        if (! eof $rh) {
            # Read the right record.
            $rightID = $self->_read_fastq($rh, 'right');
            $retVal = 1;
        }
        # Insure we have valid data.
        if ($retVal && $leftID ne $rightID) {
            if ($bufferedL) {
                # Here we are unsafe.  Simulate an empty right side and save the new left.
                $self->{unsafe} = [$rightID, $self->{right}, $self->{rqual}];
                $self->{right} = '';
                $self->{rqual} = '';
            } else {
                die "Unpaired data in FASTQ files: $leftID vs $rightID.\n";
            }
        }
    }
    # Return the success indication.
    return $retVal;
}

=head3 Out

    my $fqhandle->Out($fileName);

Open output files with the given name. The file name should not have an extension: appropriate extensions will be added.

=over 4

=item fileName

The name to give to the output files, without an extension.  The extensions C<_1.fq> and C<_2.fq> will be appended.

=back

=cut

sub Out {
    my ($self, $fileName) = @_;
    require FastQ::Out;
    return FastQ::Out->new($fileName);
}


=head3 Echo

    $fqhandle->Echo($oh);

Write the current record to the specified file handle in FASTA format.

=over 4

=item oh

An open file handle onto which the current record's sequences should be written.

=back

=cut

sub Echo {
    my ($self, $oh) = @_;
    my $id = $self->id;
    print $oh ">$id/1\n$self->{left}\n";
    if (! $self->{singleton}) {
        print $oh ">$id/2\n$self->{right}\n";
    }
}

=head3 Save

    $fqhandle->Save(\@saveList);

Save the current record to a list.

=over 4

=item saveList

Reference to a list to contain saved FASTQ records.

=back

=cut

sub Save {
    my ($self, $saveList) = @_;
    my %saveRec = map { $_ => $self->{$_} } qw(id left right lqual rqual);
    push @$saveList, \%saveRec;
}

=head3 Load

    my $id = $fqhandle->Load(\@saveList);

Load a new current record from a list and return its ID.

=over 4

=item saveList

Reference to a list containing saved FASTQ records.  The last record will be popped off.

=item RETURN

Returns the ID of the record loaded.

=back

=cut

sub Load {
    my ($self, $saveList) = @_;
    my $rec = pop @$saveList;
    if ($rec) {
        for my $key (qw(id left right lqual rqual)) {
            $self->{$key} = $rec->{$key};
        }
    }
    return $self->{id};
}

=head3 Write

    $fqhandle->Write($oh);

Write the current record to the specified file handle in interlaced FASTQ format.

=over 4

=item oh

An open file handle onto which the current record's sequences should be written.

=back

=cut

sub Write {
    my ($self, $oh) = @_;
    my $id = $self->id;
    print $oh join("\n", "\@$id/1", $self->left, "+$id/1", $self->lqual, "");
    if (! $self->{singleton}) {
        print $oh join("\n", "\@$id/2", $self->right, "+$id/2", $self->rqual, "");
    }
}

=head3 WriteL

    $fqhandle->WriteL($oh);

Write the left sequence of the current record to the specified file handle in FASTQ format.

=over 4

=item oh

An open file handle onto which the current record's sequences should be written.

=back

=cut

sub WriteL {
    my ($self, $oh) = @_;
    my $id = $self->id;
    print $oh join("\n", "\@$id/1", $self->left, "+$id/1", $self->lqual, "");
}

=head3 WriteR

    $fqhandle->WriteR($oh);

Write the right sequence of the current record to the specified file handle in FASTQ format.

=over 4

=item oh

An open file handle onto which the current record's sequences should be written.

=back

=cut

sub WriteR {
    my ($self, $oh) = @_;
    my $id = $self->id;
    print $oh join("\n", "\@$id/2", $self->right, "+$id/2", $self->rqual, "");
}


=head2 Data Access Methods

=head3 id

    my $id = $fqhandle->id;

Return the current sequence ID.

=cut

sub id {
    my ($self) = @_;
    return $self->{id};
}

=head3 left

    my $dna = $fqhandle->left;

Return the left data string.

=cut

sub left {
    my ($self) = @_;
    return $self->{left};
}

=head3 lqual

    my $dna = $fqhandle->lqual;

Return the left quality string.

=cut

sub lqual {
    my ($self) = @_;
    return $self->{lqual};
}

=head3 lq_mean

    my $qNum = $fqhandle->lq_mean;

Return the mean quality of the left string.

=cut

sub lq_mean {
    my ($self) = @_;
    return $self->_qual_mean($self->{lqual});
}

=head3 right

    my $dna = $fqhandle->right;

Return the right data string.

=cut

sub right {
    my ($self) = @_;
    return $self->{right};
}

=head3 rqual

    my $dna = $fqhandle->rqual;

Return the right quality string.

=cut

sub rqual {
    my ($self) = @_;
    return $self->{rqual};
}

=head3 rq_mean

    my $qNum = $fqhandle->rq_mean;

Return the mean quality of the right string.

=cut

sub rq_mean {
    my ($self) = @_;
    return $self->_qual_mean($self->{rqual});
}

=head3 seqs

    my @seqs = $fqhandle->seqs;

Return a list of the sequences stored in the object. (There will always be two.)

=cut

sub seqs {
    my ($self) = @_;
    return ($self->{left}, $self->{right});
}

=head3 lqual_mean


=head2 Internal Utilities

=head3 _qual_mean

    my $pct = $fqhandle->_qual_mean($qualString);

Compute the mean quality from a quality string.

=over 4

=item qualString

The quality string for a sequence.

=item RETURN

Returns a fraction from 0 to 1 indicating the mean quality.

=back

=cut

use constant QBASE => 0x20;

sub _qual_mean {
    my ($self, $qualString) = @_;
    my $base = $self->{q_base};
    my $n = length $qualString;
    my $retVal = 0;
    if ($n > 0) {
        my $total = 0;
        for (my $i = 0; $i < $n; $i++) {
            my $q = ord(substr($qualString, $i, 1)) - $base;
            $q = 40 if $q > 40;
            $total += $q;
        }
        $retVal = $total / ($n * 40);
    }
    return $retVal;
}

=head3 _read_fastq

    my $id = $fqhandle->_read_fastq($ih, $dir);

Read the next record from the indicated FASTQ input stream and store its data in the specified object members.

=over 4

=item ih

Open file handle for the FASTQ file.

=item dir

C<left> for a left record and C<right> for a right record.

=item RETURN

Returns the ID for the read.

=back

=cut

sub _read_fastq {
    my ($self, $ih, $dir) = @_;
    # Compute the quality string member name.
    my $qual = substr($dir, 0, 1) . 'qual';
    # This will hold the sequence ID.
    my $id;
    # Read the ID line.
    my $line = <$ih>;
    if ($line =~ /^\@(\S+)/) {
        $id = $1;
    } else {
        die "Invalid FASTQ input: $line";
    }
    # Read the data line.
    my $data = <$ih>;
    die "Missing data line for $id in FASTQ file." if ! defined $data;
    chomp $data;
    $self->{$dir} = $data;
    # Read the quality header.
    $line = <$ih>;
    if (! $line || substr($line, 0, 1) ne '+') {
        die "Invalid quality header for $id in FASTQ file.";
    } else {
        # Read the quality data.
        $line = <$ih>;
        die "Missing quality line for $id in FASTQ file." if ! defined $line;
        chomp $line;
        if (length($line) ne length($data)) {
            die "Incorrect length for quality line belonging to $id in FASTQ file.";
        } else {
            $self->{$qual} = $line;
        }
    }
    # Normalize the ID and store it.
    my $retVal = norm_id($id);
    $self->{id} = $retVal;
    # Return the ID.
    return $retVal;

}


1;
