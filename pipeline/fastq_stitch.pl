#!/usr/bin/perl -w
#
# fastq_stitch.pl
# Use a defined spacer sequence to stitch together Read1 and Read2 data from paired read NGS data sets.

use strict;
use warnings;

sub process_fastq_block    # read in FASTQ block and return processed portions
                           # 1st argument -- filehandle for reading inputFile
                           # 2nd argument -- filename
                           # returns array containing:
                           #   [0] success flag
                           #   [1] sequence id
                           #   [2] sequence
                           #   [3] quality
{
    my ($fh, $fname) = @_; # unpacking function arguments
    my @output = ( 0, '', '', '' );
    my $line;
    my $id;
    my $seq;
    my $qual;

    if ( !defined( $line = <$fh> ) ) { return @output; }

    chomp $line;
    if ( $line =~ m|^(\@.*)|x ) { $id = $1; }
    else { die "Error: invalid FASTQ format at $line.\n"; }

    if ( !defined( $seq = <$fh> ) ) {
        die "Error: check input FASTQ format at $id\n";
    }
    chomp $seq;

    if ( !defined( $line = <$fh> ) ) {
        if ( $line != m|^\+$|x ) {
            die "Error: check input FASTQ format at $id\n";
        }
    }
    if ( !defined( $qual = <$fh> ) ) {
        die "Error: check input FASTQ format at $id\n";
    }
    chomp $qual;

    if ( (length $seq) != (length $qual) ) # sanity check: each base should have a quality call
      {     die "Error: invalid FASTQ format in $fname at $id.\n"; } # exit on error

    @output = ( 1, $id, $seq, $qual );
    return @output;
}

sub reverse_complement    # read in a sequence and return reverse-complement
                          # 1st argument -- sequence (string)
                          # returns reverse-complement (string)
{
    my ($seq) = @_;
    chomp($seq);
    $seq =~ tr /atcgATCG/tagcTAGC/;
    my $revComp = reverse($seq);
    return $revComp;
}

local $/ = "\n";

my $filename1;
my $filename2;
my $spacerSeq;

# process commandline
if ( scalar @ARGV == 2 || scalar @ARGV == 3 ) {
    $filename1 = $ARGV[0];
    $filename2 = $ARGV[1];

    if ( scalar @ARGV == 3 ) {
        $spacerSeq = $ARGV[2];
    }
    else {
        $spacerSeq = 'NNNNNN';
    }
}
else {
    die "Usage: perl fastq_stitch.pl read1.fastq read2.fastq stitch_sequence\n";
}

my $spacerQual = '#' x length($spacerSeq);

open my $infile1, '<:encoding(UTF-8)', "$filename1"    # read1 file
  or die "Cannot open $filename1.\n";

open my $infile2, '<:encoding(UTF-8)', "$filename2"    # read2 file
  or die "Cannot open $filename2.\n";

my @read1_datablock = process_fastq_block($infile1, $filename1);
my @read2_datablock = process_fastq_block($infile2, $filename2);

while ( $read1_datablock[0] ) {
    $read1_datablock[1] =~
      s|\s.*$||x;    # remove the read-specific portion of the identifiers
    $read2_datablock[1] =~ s|\s.*$||x;

    if ( $read1_datablock[1] ne $read2_datablock[1] ) {
        die
          "Error: mismatched identifiers in:\n$filename1: $read1_datablock[1]\n$filename2: $read2_datablock[1]\n";
    }

    print $read1_datablock[1], "\n";
    print $read1_datablock[2]
      . $spacerSeq
      . reverse_complement( $read2_datablock[2] ), "\n";
    print "\+\n";
    print $read1_datablock[3] . $spacerQual . reverse( $read2_datablock[3] ),
      "\n";
    @read1_datablock = process_fastq_block($infile1, $filename1);
    @read2_datablock = process_fastq_block($infile2, $filename2);
}

close $infile1;
close $infile2;
