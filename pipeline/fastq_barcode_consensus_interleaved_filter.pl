#!/usr/bin/perl -w
#
# fastq_barcode_consensus_interleaved_filter.pl
# Use UMI(barcode)-annotated data set to prepare strictly interleaved input file for FLASH processing

use strict;
use warnings;

my $filename;

sub process_fastq_block # read in FASTQ block and return processed portions
                  # 1st argument -- filehandle for reading inputFile
                  # 2nd argument -- filename
                  # returns array containing:
                  #   [0] success flag
                  #   [1] sequence id
                  #   [2] orientation
                  #   [3] retention report
                  #   [4] number of sequences retained
                  #   [5] sequence
                  #   [6] quality
{
  my ($fh, $fname) = @_; # unpacking function arguments
  my @output = ( 0, '', '', '', '', '', '');
  my $line;
  my $id;
  my $orientation;
  my $retention_report;
  my $num_retained;
  my $seq;
  my $qual;

  if(!defined( $line = <$fh> ))
  { return @output; }

  chomp $line;
  if ( $line =~ m|^(\@MIG\d+;)[^;]+;orient_(.{3});(barcode=[ATCG]+);valid;size=(\d+);retained=(\d+)|x ) {
      $id = $1 . $3;
      $orientation = $2;
      $retention_report = $4 . ':' . $5;
      $num_retained = $5;
  }
  else { die "Error: the following ID is incorrectly formatted:\n$line\n";}

  if ( !defined( $seq = <$fh> ) ) {
      die "Error: check input FASTQ format at $id\n";
  }

  if ( defined( $line = <$fh> ) ) {
      if ( $line ne "+\n" ) {
          die "Error: check input FASTQ format at $id\n";
      }
  }
  else { die "Error: reached an abrupt end of FASTQ input at $id\n";}

  if ( !defined( $qual = <$fh> ) ) {
      die "Error: reached an abrupt end of FASTQ input at $id\n";
  }

  if ( (length $seq) != (length $qual) ) # sanity check: each base should have a quality call
    {     die "Error: invalid FASTQ format in $fname at $id.\n"; } # exit on error

  @output = (1, $id, $orientation, $retention_report, $num_retained, $seq, $qual);
  return @output;
}

# process commandline
if ( scalar @ARGV == 1 ) {
    $filename = $ARGV[0];
}
else {
    die
      "Usage: perl fastq_barcode_consensus_interleaved_filter.pl source.fastq\n";
}

my @current_datablock = (1,'','','','','','');
my @fwd_data          = @current_datablock;
my @rev_data          = @current_datablock;
my $identifier;

# my $file_read_successful = $current_datablock[0];
open my $file, '<:encoding(UTF-8)', "$filename"
  or die "Cannot open $filename.\n";

while ($current_datablock[0]) {
  @current_datablock = process_fastq_block ($file, $filename);

  if ($current_datablock[2] eq 'fwd'){
    @fwd_data = @current_datablock;
  }
  elsif($current_datablock[2] eq 'rev'){
    @rev_data = @current_datablock;
  }

  if($fwd_data[1] eq $rev_data[1]){
    $identifier = $fwd_data[1] . ';fwd:' . $fwd_data[3] . ':rev:' . $rev_data[3];
    $identifier .= ';retained=' . ($fwd_data[4] + $rev_data[4]);
    print "$identifier\n$fwd_data[5]+\n$fwd_data[6]";
    print "$identifier\n$rev_data[5]+\n$rev_data[6]";
  }
}
close $file;
