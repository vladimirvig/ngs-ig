#!/usr/bin/perl -w
#
# clonotype_annotate.pl
# take a clonotype dictionary and annotate its source FASTA file with clonotype labels

use strict;
use warnings;

if ( scalar @ARGV != 2 ) {
    die
"Usage: perl clonotype_annotate.pl dictionary.clonotype_dict source.fasta\n";
}


my $dictfilename = $ARGV[0];
my $srcfilename = $ARGV[1];

my %dictionary;
my $current;
my $description;

open my $dictfile, '<:encoding(UTF-8)', "$dictfilename"
  or die "Cannot open $dictfilename.\n";

# load dictionary
while ( defined( $current = <$dictfile> ) ) {
    if (
        $current =~ m|^(\w{5})\t  # clonotype ID
                    (\w+)\t       # CDR3aa
                    \d+\t         # number of V/J pairs on the line
                    (\d+)\t       # number of sequences in the set (excluding chimeras)
                    V(\d+)\t      # V family
                    [jJ](\d+)     # J family
                  |x
      )
    {
        $dictionary{$2}[0] = $1;    # store clonotype ID
        $dictionary{$2}[1] = $4;    # store dominant V family
        $dictionary{$2}[2] = $5;    # store dominant J family
        $dictionary{$2}[3] = $3;    # store clonotype population size
    }
    elsif (
        $current =~ m|^(\w{5})\t  # clonotype ID
                        (\w+)\t   # CDR3aa
                        \d+\t     # number of V/J pairs on the line
                        (\d+)     # number of sequences in the set (excluding chimeras)
                        .+\wambig # V/Jambig flag
                      |x
      )    # deal with ambiguous (empty) clonotypes
    {
        $dictionary{$2}[0] = $1;
        $dictionary{$2}[1] = 'ambig';
        $dictionary{$2}[2] = 'ambig';
        $dictionary{$2}[3] = $3;        # store clonotype population size
    }
    elsif ( $current =~ m|^#\s|x )      # skip comments in input
    { }
    else { die "!!!Error: Unknown $current\n"; }
}
close $dictfile;

open my $srcfile, '<:encoding(UTF-8)', "$srcfilename"
  or die "Cannot open $srcfilename.\n";

# annotate the dataset
while ( defined( $current = <$srcfile> ) ) {
    if ( $current =~ m|^\>(\S+)\t(.+)$|x ) {
        print ">", $1;
        $description = $2;

        if (
            $description =~ m|V[[:alpha:]]?(\d+)            # V-segment assignment
                                \S+\t?\S*\t                 # possible D-segment assignment
                                \w+[jJ][[:alpha:]]?(\d+).+  # J-segment assignment
                                CDR3aa\:\w\w(\w+)\w         # capture the CDR3 without
                                                            # the leading and trailing aa's
                             |x
          )
        {
            if ( !defined $dictionary{$3} ) {
                die "Error: failed to look up $3.\n";
            }

            print ";",
              $dictionary{$3}[0];    # label the sequence ID with clonotype ID
            print "-", $dictionary{$3}[3]
              ;    # append the clonotype population size to clonotype ID
            if (   ( $dictionary{$3}[1] eq 'ambig' )
                or ( $dictionary{$3}[1] != $1 )
                or ( $dictionary{$3}[2] != $2 ) )
            {
                print "-chimera";
            }      # mark as chimeric if the major V and J don't match
        }
        else { die "Error parsing $description.\n"; }

        print "\t", $description, "\n";

        #fetch and print the sequence
        $current = <$srcfile>;
        print $current;
    }
}

close $srcfile;
