#!/usr/bin/perl -w
#
# fastx_asym_orientation_fix.pl
# Use the forward/reverse primer information to uniformly orient annotated sequences in an NGS data set

use strict;
use warnings;

local $/ = "\n";

my $infilename;
my $fwd_primer_name;
my $rev_primer_name;
my $DEBUG;

if ( scalar @ARGV == 3 || scalar @ARGV == 4 ) {
    $infilename    = $ARGV[0];
    $fwd_primer_name = $ARGV[1];
    $rev_primer_name = $ARGV[2];

    if ( scalar @ARGV == 4 ) {
        $DEBUG = $ARGV[3];
    }
}
else {
    die "Usage: perl fastq_asym_orientation_fix.pl trimmed_annotated.fastq/a FWD_primer_name REV_primer_name [DEBUG]\n";
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

my $line;
my $header;
my $sequence;
my $quality;

my $readLabel;
my $primer1;
my $primer2;
my $fastq_flag;

my $orientation;
my $count = 0;

if    ( $infilename =~ m|\.fastq$|x ) { $fastq_flag = 1; }
elsif ( $infilename =~ m|\.fasta$|x ) { $fastq_flag = 0; }
else { die "Error: invalid filetype (neither FASTQ not FASTA)."; }

open my $infile, '<:encoding(UTF-8)', "$infilename"    # annotated read1
  or die "Cannot open $infilename.\n";

while ( defined( $header = <$infile> ) ) {
    ############## FASTX harvest block ##############
    if ( !defined( $sequence = <$infile> ) ) {
        die "Error: couldn't retrieve sequence from $header\n";
    }
    if ($fastq_flag) {
        if ( !defined( $line = <$infile> ) ) {
            die "Error: FASTQ format error for $header\n";
        }
        if ( !defined( $quality = <$infile> ) ) {
            die "Error: couldn't retrieve quality scores from $header\n";
        }
    }
    ############## FASTX harvest block ##############

    chomp($header);
    chomp($sequence);
    if ($fastq_flag) { chomp($quality); }
    $orientation = "orient_unk";

    if ( $header !~ m|^(.*);([^;]*);([^;]*)$|x ) {
        die "Error: unrecognized file format. Make sure the primer annotation is present.\n";
    }
    else {
        $readLabel = $1;
        $primer1   = $2;
        $primer2   = $3;
    }

    if ( $primer1 =~ m|$fwd_primer_name|x ) {
        if ( $primer1 eq $fwd_primer_name . "_adapter" ) {
            if ( $primer2 eq "$rev_primer_name" . "_revcomp_adapter" ) {
                $orientation = 'orient_fwd';
            }
        }
    }
    elsif ( $primer1 =~ m|$rev_primer_name|x ) {
        if ( $primer1 eq $rev_primer_name . "_adapter" ) {
            if ( $primer2 eq "$fwd_primer_name" . "_revcomp_adapter" ) {
                $orientation = 'orient_rev';
            }
        }
    }
    else { $orientation = 'orient_unk'; }

    if ( $orientation eq 'orient_rev' ) {
        $sequence    = reverse_complement($sequence);
        $orientation = "orient_fwd";
        $count++;
        if ($fastq_flag) { $quality = reverse($quality); }
    }

    print "$readLabel;$primer1-$primer2;$orientation\n$sequence\n";
      if ($fastq_flag) { print "\+\n$quality\n"; }
}
close $infile;

if ( defined $DEBUG && $DEBUG eq 'DEBUG' ) {
    print ">DEBUG: reversed $count reads.\nNNNNN\n";
}
