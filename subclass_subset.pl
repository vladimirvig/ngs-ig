#!/usr/bin/perl -w
#
# subclass_subset.pl
#
# This script generates a list of different subclass representatives from each cluster

use strict;
use warnings;

local $/ = "\n";

my $filename;
# process commandline
if ( scalar @ARGV == 1 ) {
    $filename = $ARGV[0];
}
else {
    die
      "Usage: perl subclass_subset.pl subclass_annotated.fasta\n";
}

my $line;
my %bins;
my $cluster_id;
my $subclass_gene;

# need top 10 representatives for each subclass from a given cluster
# the assumption is that the sequence set is ordered by abundance
open my $infile, '<:encoding(UTF-8)', "$filename"    # blast_out file
  or die "Cannot open $filename.\n";

while ( defined( $line = <$infile> ) ) {
    if ( $line =~ m|^>\d+\-\d+;(\w+)\-.+subclass:(IGHG[^\*]+)|x ) {
        $cluster_id    = $1;
        $subclass_gene = $2;

        ++$bins{"$cluster_id:$subclass_gene"};

        if (   ( $bins{"$cluster_id:$subclass_gene"} < 11 )
            || ( $subclass_gene =~ m|^IGHG3$|x ) )
        {
            print $line;
            $line = <$infile>;
            print $line;
        }
    }
}

close $infile;
