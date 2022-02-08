#!/usr/bin/perl -w
#
# fastq_asym_barcode_order.pl
# Order sequences in an NGS data set according to decreasing UMI (barcode) abundance
#  Also, assign element number to the sequences with the same UMI.
#  Write the invalid sequences (with unknown barcodes) in the end.

use strict;
use warnings;

local $/ = "\n";

my $filename;
my $debug_flag;

# process commandline
if ( scalar @ARGV == 1 || scalar @ARGV == 2 ) {
    $filename = $ARGV[0];
    if ( scalar @ARGV == 2 ) {
            $debug_flag = $ARGV[1];
    }
}
else {
    die
      "Usage: perl fastq_asym_barcode_order.pl barcode_annotated.fastq [DEBUG]\n";
}

my $line;
my $reference;
my $identifier;
my $fwd_adapter_name;
my $rev_adapter_name;

my $barcode;
my %barcodes;
my @keys;

my %fwd_reads;
my %rev_reads;
my @fwd_keys;
my @rev_keys;
my %unknowns;

my $garbage_counter = 0;
my $total_counter   = 0;
my $orphan_barcodes = 0;
my $counter        = 0;
my $id             = 1;

open my $infile, '<:encoding(UTF-8)', "$filename"    # annotated read1
  or die "Cannot open $filename.\n";

while ( defined( $line = <$infile> ) ) {

    # processing the first line in the FASTQ block
    if ( $line =~ m|^\@.+(;.+;orient_\w{3};barcode=\w+)$|x ) {
        $identifier = $1;
        $reference  = tell $infile;

        if ( $identifier =~ m|^;(.+)_adapter-(.+)_adapter|x ) {
            $fwd_adapter_name = $1;
            $rev_adapter_name = $2;
        }
        else { die "Error: unrecognized adapter annotation in $identifier\n"; }

        if (   ( $fwd_adapter_name eq 'no' )
            or ( $rev_adapter_name eq 'no' )
            or ( $fwd_adapter_name eq $rev_adapter_name ) )
        {
            $garbage_counter++;
            $unknowns{$identifier}[0]++;
            push @{ $unknowns{$identifier} }, $reference;
        }
        elsif ( ( $identifier =~ m|;orient_fwd;barcode=([ATCG]{16})$|x ) ) {
            $barcode = $1;
            $barcodes{$barcode}++;
            $fwd_reads{$barcode}[0]++;

            if ( $fwd_reads{$barcode}[0] == 1 ) {
                push @{ $fwd_reads{$barcode} }, $identifier;
            }
            elsif ( $fwd_reads{$barcode}[1] !~ m|^$identifier$|x ) {
                die
                  "Error: the $identifier doesn't match the storage target $fwd_reads{$barcode}[1]\n";
            }

            push @{ $fwd_reads{$barcode} }, $reference;
        }
        elsif ( ( $identifier =~ m|;orient_rev;barcode=([ATCG]{16})$|x ) ) {
            $barcode = $1;
            $barcodes{$barcode}++;
            $rev_reads{$barcode}[0]++;

            if ( $rev_reads{$barcode}[0] == 1 ) {
                push @{ $rev_reads{$barcode} }, $identifier;
            }
            elsif ( $rev_reads{$barcode}[1] !~ m|^$identifier$|x ) {
                die
                  "Error: the $identifier doesn't match the storage target $rev_reads{$barcode}[1]\n";
            }

            push @{ $rev_reads{$barcode} }, $reference;
        }
        else {
            $garbage_counter++;
            $unknowns{$identifier}[0]++;
            push @{ $unknowns{$identifier} }, $reference;
        }

        $total_counter++;
    }
    else {
        die
          "Error: did not find the expected annotations in FASTQ format at\n$line";
    }

    # skipping the rest of lines in the FASTQ block: sequence, +, quality
    $counter = 3;
    while ( $counter > 0 ) {
        if ( !defined( $line = <$infile> ) ) {
            die
              "Error: problem with FASTQ format below header for: $identifier\n";
        }
        $counter--;
    }
}

# order the barcodes by abundance
@keys = sort { $barcodes{$b} <=> $barcodes{$a} } ( keys %barcodes );

@fwd_keys = keys %fwd_reads;
@rev_keys = keys %rev_reads;

foreach (@keys) {
    if ( defined $fwd_reads{$_} ) {
        $counter    = $fwd_reads{$_}[0];
        $identifier = $fwd_reads{$_}[1];

        while ($counter) {
            print "@", $id, $identifier, ";valid;size=", $fwd_reads{$_}[0],
              ";element=", $counter, "\n";
            $reference = pop @{ $fwd_reads{$_} };
            seek( $infile, $reference, 0 );

            if ( !defined( $line = <$infile> ) ) {
                die "Error at barcode $fwd_reads{$_}[0].\n";
            }
            print $line;
            if ( !defined( $line = <$infile> ) ) {
                die "Error at barcode $fwd_reads{$_}[0].\n";
            }
            print $line;
            if ( !defined( $line = <$infile> ) ) {
                die "Error at barcode $fwd_reads{$_}[0].\n";
            }
            print $line;

            $counter--;
        }
    }

    if ( defined $rev_reads{$_} ) {
        $counter    = $rev_reads{$_}[0];
        $identifier = $rev_reads{$_}[1];

        while ($counter) {
            print "@", $id, $identifier, ";valid;size=", $rev_reads{$_}[0],
              ";element=", $counter, "\n";
            $reference = pop @{ $rev_reads{$_} };
            seek( $infile, $reference, 0 );

            if ( !defined( $line = <$infile> ) ) {
                die "Error at barcode $rev_reads{$_}[0].\n";
            }
            print $line;
            if ( !defined( $line = <$infile> ) ) {
                die "Error at barcode $rev_reads{$_}[0].\n";
            }
            print $line;
            if ( !defined( $line = <$infile> ) ) {
                die "Error at barcode $rev_reads{$_}[0].\n";
            }
            print $line;

            $counter--;
        }
    }
    $id++;
}

@keys = keys %unknowns;

foreach (@keys) {
    $counter    = $unknowns{$_}[0];
    $identifier = $_;

    while ($counter) {
        print "@", $id, $identifier, ";invalid;size=", $unknowns{$_}[0],
          ";element=", $counter, "\n";
        $reference = pop @{ $unknowns{$_} };
        seek( $infile, $reference, 0 );

        if ( !defined( $line = <$infile> ) ) {
            die "Error at barcode $unknowns{$_}[0].\n";
        }
        print $line;
        if ( !defined( $line = <$infile> ) ) {
            die "Error at barcode $unknowns{$_}[0].\n";
        }
        print $line;
        if ( !defined( $line = <$infile> ) ) {
            die "Error at barcode $unknowns{$_}[0].\n";
        }
        print $line;

        $counter--;
    }
    $id++;
}

close $infile;

# Accounting:
if ( defined $debug_flag && $debug_flag eq 'DEBUG' ) {
    foreach (@keys) {
        if ( ( !defined $fwd_reads{$_} ) or ( !defined $rev_reads{$_} ) ) {
            $orphan_barcodes++;
        }
    }

    print "Total: $total_counter\n";
    @fwd_keys = keys %fwd_reads;
    print "Forward: ", scalar @fwd_keys, "\n";
    @rev_keys = keys %rev_reads;
    print "Reverse: ",        scalar @rev_keys, "\n";
    print "Total Barcodes: ", scalar @keys,    "\n";
    print "Orphans: $orphan_barcodes\n";
    print "Garbage: $garbage_counter\n";
}
