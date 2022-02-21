#!/usr/bin/perl -w
#
# fastq_asym_barcode_transfer.pl
# Transfer UMI (barcode) from Read2 to Read1 in a paired-read NGS library data set
#  that has been sequenced asymmetrically (Read2 data will ultimately be discarded).

use strict;
use warnings;

local $/ = "\n";
my $usage =
    "Usage: perl fastq_asym_barcode_transfer.pl adapter_annotated_read1.fastq adapter_annotated_read2.fastq \'fwd_adapter_name\' \'rev_adapter_name\' 'preamble\' \'barcode\' \'post\'.\n";

my $filename1;
my $filename2;
my $fwd_adapter_name;
my $rev_adapter_name;

my $preamble_pattern;
my $barcode_pattern;
my $post_pattern;

# processing commandline arguments
if ( scalar @ARGV == 7 ) {
    $filename1          = $ARGV[0];
    $filename2          = $ARGV[1];
    $fwd_adapter_name   = $ARGV[2];
    $rev_adapter_name   = $ARGV[3];
    $preamble_pattern   = $ARGV[4];
    $barcode_pattern    = $ARGV[5];
    $post_pattern       = $ARGV[6];
}
else { die $usage; }

my $umi = qr<($preamble_pattern)($barcode_pattern)($post_pattern)>;

my $read1_fh;
my $read2_fh;

my $header1;
my $header2;

my $label;

my $sequence1;
my $sequence2;

my $quality1;
my $quality2;

my $adapter1;
my $adapter2;

my $preamble_seq;
my $barcode_seq;
my $post_seq;
my $trim_seq;
my $trim_seq_len;
my $line;

my $orientation;
open $read1_fh, '<:encoding(UTF-8)', "$filename1"    # annotated read1
  or die "Cannot open $filename1.\n";

open $read2_fh, '<:encoding(UTF-8)', "$filename2"    # annotated read2
  or die "Cannot open $filename2.\n";

while ( defined( $header1 = <$read1_fh> ) ) {
    ############## FASTQ1 harvest block ##############
    if ( !defined( $sequence1 = <$read1_fh> ) ) {
        die "Error: couldn't retrieve sequence from $header1\n";
    }
    if ( !defined( $line = <$read1_fh> ) ) {
        die "Error: FASTQ format error for $header1\n";
    }
    if ( !defined( $quality1 = <$read1_fh> ) ) {
        die "Error: couldn't retrieve quality scores from $header1\n";
    }
    ############## FASTQ1 harvest block ##############

    ############## FASTQ2 harvest block ##############
    if ( !defined( $header2 = <$read2_fh> ) ) {
        die "Error: reached the end of $filename1 at\n$header1\n\tbefore\n$filename2\n";
    }
    if ( !defined( $sequence2 = <$read2_fh> ) ) {
        die "Error: couldn't retrieve sequence from $header2\n";
    }
    if ( !defined( $line = <$read2_fh> ) ) {
        die "Error: FASTQ format error for $header2\n";
    }
    if ( !defined( $quality2 = <$read2_fh> ) ) {
        die "Error: couldn't retrieve quality scores from $header2\n";
    }
    ############## FASTQ2 harvest block ##############

    ############### FASTQ label check ################
    if ( $header1 =~ m|^(@[\S]+)|x )    # harvest the read1 label for validation
    { $label = $1; }
    else { die "Error reading label:\n$header1\n"; }

    if ( $header2 =~ m|^(@[\S]+)|x ) {
        if ( $1 ne $label ) {
            die "Error: label mismatch:\n$header1\n$header2\n";
        }
    }
    ############### FASTQ label check ################

    if ( $header1 =~
        m|^@.+;([\w\_\-]+)_adapter$|x )    # harvest the read1 adapter sequence
    { $adapter1 = $1; }
    else {
        die
          "Error with adapter annotation in the following pair:\n$header1\n$header2\n";
    }

    if ( $header2 =~
        m|^@.+;([\w\_\-]+)_adapter$|x )    # harvest the read2 adapter sequence
    { $adapter2 = $1; }
    else {
        die
          "Error with adapter annotation in the following pair:\n$header1\n$header2\n";
    }

    if ( $adapter1 eq $fwd_adapter_name )    # looks like a forward read
    {
        if ( $adapter2 eq $rev_adapter_name )    # forward-reverse
        { $orientation = "orient_fwd"; }
        else { $orientation = "orient_unk"; }

        if ( $sequence1 =~ m|^$umi|x )          # harvest the UMI barcode
        {
            $preamble_seq = $1;
            $barcode_seq  = $2;
            $post_seq     = $3;

            $trim_seq       = $preamble_seq . $barcode_seq;
            $trim_seq_len = length($trim_seq);
            $sequence1 =~ s|^$trim_seq||x; # trim the barcode from the sequence
            $quality1  =~ s|^.{$trim_seq_len}||x
              ;    # trim the barcode-length string from the quality string
        }
        else { $barcode_seq = "unknown"; }
    }
    elsif (
        $adapter2 eq $fwd_adapter_name )    # check if read2 has a forward adapter
    {
        if ( $adapter1 eq $rev_adapter_name )    # reverse-forward
        { $orientation = "orient_rev"; }
        else { $orientation = "orient_unk"; }

        if ( $sequence2 =~ m|^$umi|x )          # harvest the UMI barcode
        {
            $preamble_seq = $1;
            $barcode_seq  = $2;
            $post_seq     = $3;
        }
        else { $barcode_seq = "unknown"; }
    }
    else    # all failed
    {
        $orientation = "orient_unk";
        $barcode_seq     = "unknown";
    }

    # write output by adding informartion to the read1 header
    chomp $header1;
    print $header1, "-$adapter2\_adapter;$orientation;barcode=$barcode_seq\n";
    print $sequence1;
    print "+\n";
    print $quality1;

}

close $read1_fh;
close $read2_fh;
