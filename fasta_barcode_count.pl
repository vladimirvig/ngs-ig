#!/usr/bin/perl -w
#
# barcode_count.pl
#
# This script uses the known adapter sequences from UMI-tagged NGS library data
#  to extract the UMI barcodes for downstream processing.
use strict;
use warnings;

# process commandline
if ( scalar @ARGV != 4 ) {
    die
      "Usage: perl barcode_count.pl source.fasta \'preamble\' \'barcode\' \'post\'\n";
}

my $source_filename = $ARGV[0];
my $preamble        = $ARGV[1];
my $barcode         = $ARGV[2];
my $post            = $ARGV[3];

my %barcodes;
my %unknowns;
my @keys;

my $line;
my $reference;
my $id = 1;
my $index;

my $umi = qr<$preamble($barcode)$post>;

open my $source_file, '<:encoding(UTF-8)', "$source_filename"
  or die "Cannot open $source_filename.\n";

while ( defined( $line = <$source_file> ) ) {
    if ( $line =~ m|^\>|x ) { $reference = tell $source_file; }
    else {
        if ( $line =~ m|^$preamble|x ) {
            if ( $line =~ m|^$umi|x ) {
                $barcodes{$1}[0]++;
                push @{ $barcodes{$1} }, $reference;
            }
            else {
                $unknowns{'unknown'}[0]++;
                push @{ $unknowns{'unknown'} }, $reference;
            }
        }
        else {
            $unknowns{'garbage'}[0]++;
            push @{ $unknowns{'garbage'} }, $reference;
        }
    }
}

@keys = sort { $barcodes{$b}[0] <=> $barcodes{$a}[0] } ( keys %barcodes );

# The Main Set
foreach (@keys) {
    $index = $barcodes{$_}[0];

    while ($index) {
        print ">MIG", $id, ";barcode=", $_, ";size=", $barcodes{$_}[0],
          ";element=", $index, "\n";
        $reference = pop @{ $barcodes{$_} };
        seek( $source_file, $reference, 0 );
        if ( !defined( $line = <$source_file> ) ) {
            die "Error at barcode $barcodes{$_}[0].\n";
        }
        $line =~ s|^$umi||x;
        print $line;
        $index--;
    }
    $id++;
}

@keys = keys %unknowns;

# The unknowns (trash)
foreach (@keys) {
    $index = $unknowns{$_}[0];

    while ($index) {
        print ">MIG", $id, ";barcode=", $_, ";size=", $unknowns{$_}[0],
          ";element=", $index, "\n";
        $reference = pop @{ $unknowns{$_} };
        seek( $source_file, $reference, 0 );
        if ( !defined( $line = <$source_file> ) ) {
            die "Error at barcode $unknowns{$_}[0].\n";
        }
        print $line;
        $index--;
    }
    $id++;
}

close $source_file;
