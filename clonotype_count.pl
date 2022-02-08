#!/usr/bin/perl -w
#
# clonotype_count.pl
#  take an annotated FASTA file from immunoglobulin sequences, carry out accounting,
#  and generate a clonotype dictionary

use strict;
use warnings;

srand 12345;    # Set random seed for output reproducibility.

#random string generator
sub rndStr {
    my $size = shift;
    return join '', @_[ map { rand @_ } 1 .. $size ];
}

if ( scalar @ARGV != 1 ) {
    die "Usage: perl clonotype_count.pl source.fasta\n";
}
my $src_filename = $ARGV[0];

my %clonotype;
my @keys;
my %dictionary;
my $current;
my $clust_id;
my %clust_ids;
my $index;
my $vj_assignment;
my $v_fam;
my $j_fam;
my $rearr;
my $vj_rearr_code;

my $vj_tally;

open my $src_file, '<:encoding(UTF-8)', "$src_filename"
  or die "Cannot open $src_filename.\n";

while ( defined( $current = <$src_file> ) ) {
    if (
        $current =~ m|^>[^\t]+[-=](\d+)\t
                       ([^\s\,]*V[^\s\,]+)\S*\t
                       ([^\s\,]*[DN][^\s\,]+\S+\t)?
                       ([^\s\,]*[jJ][^\s\,]+).+CDR3aa\:(\w{2})(\w+)(\w){1}
                    |x
      )
    {
        if ( ( $5 . $6 . $7 ) ne '0null0' ) {
            $vj_rearr_code = $2 . ';' . $4 . ';' . $6;
        }    # collecting Vfam, Jfam and CDR3aa
        else {
            $vj_rearr_code = $2 . ';' . $4 . ';' . '0null0';
        }    # keep the invalid CDR3 for future annotation

        if ( !defined( $clonotype{$vj_rearr_code} ) ) {
            $clonotype{$vj_rearr_code}[0] = 1;
            $clonotype{$vj_rearr_code}[1] = $1;

            #$clonotype{$2}[2] = $3;
        }
        else {
            $clonotype{$vj_rearr_code}[0]++;
            $clonotype{$vj_rearr_code}[1] += $1;
        }

        $current = <$src_file>;
    }
    else { die "Problem parsing $current.\n"; }
}

@keys = sort { $clonotype{$b}[1] <=> $clonotype{$a}[1] } ( keys %clonotype );

print "# There were ", scalar @keys, " clonotypes detected.\n";

# load the dictionary with rearrangements from the clonotypes hash.
foreach (@keys) {
    $current = $_;
    if ( $current =~ m|([^\s;]*;[^\s;]*);([^\s;]*)|x ) {
        $vj_assignment = $1;
        $rearr        = $2;

        if (
            $vj_assignment =~ m|V[[:alpha:]]?(\d+).*\;.*[jJ][[:alpha:]]?(\d+)|x )
        {
            $v_fam = $1;
            $j_fam = $2;
        }
        else { die "Error interpreting $vj_assignment.\n"; }
    }
    else { die "Error with: $current\n"; }

    if ( !defined( $dictionary{$rearr} )
      ) # the first one is expected to be dominant due to the clonotype sort; check later.
    {
        # assign a unique identifier
        $clust_id =
          rndStr( 5, 'a' .. 'z', 'A' .. 'Z', 0 .. 9 );    # generate identifier
        $clust_ids{$clust_id}++;    # initialize counter for the identifier

        $index = 10
          ;  # in cases when new random string is not unique, try up to 10 times
        if ( $clust_ids{$clust_id} > 1 )    # if not unique ...
        {
            while ( $clust_ids{$clust_id} and $index-- ) {
                $clust_id = rndStr( 5, 'a' .. 'z', 'A' .. 'Z', 0 .. 9 );
            }                             #try to generate another identifier
            if ( !($index) ) {
                die "# Error! failed to get a unique cluster ID.";
            }
        }

        $dictionary{$rearr}[0] =
          $clust_id;    #assign identifier that has been verified unique
        $dictionary{$rearr}[1] = $rearr;
        $dictionary{$rearr}[2] =
          1;           # store the number of VJ pairs for this rearrangement
        $dictionary{$rearr}[3] = $clonotype{$_}[1]
          ;    # store the total number of sequences for this rearrangement
        $dictionary{$rearr}[4] =
          $v_fam;    # store dominant V-family for this rearrangement
        $dictionary{$rearr}[5] =
          $j_fam;    # store dominant J-family for this rearrangement

        $dictionary{$rearr}[6] = $vj_assignment;      # store current VJ pair
        $dictionary{$rearr}[7] = $clonotype{$_}[1]
          ;    # store the total number of sequences for current VJ pair
    }
    else {
        $index = $dictionary{$rearr}[2]
          ;    # figure out where the new clonotype will be added
        $dictionary{$rearr}[2]++
          ;    # bump up the total number of clonotypes with this rearrangement
        $dictionary{$rearr}[3] += $clonotype{$_}[1]
          ;  # bump up the number of sequences represented by this rearrangement

        $dictionary{$rearr}[ 6 + $index * 2 ] =
          $vj_assignment;    # store current VJ pair
        $dictionary{$rearr}[ 7 + $index * 2 ] = $clonotype{$_}[1]
          ;    # store the total number of sequences for current VJ pair
    }
}

@keys = keys %dictionary;

# check dominant clonotype and mark the ambiguous V/J assignments
foreach (@keys) {
    if ( $dictionary{$_}[7] <= ( $dictionary{$_}[3] / 2 )
      ) # if the first clonotype is not the majority, check the V/J majority assignment
    {
        $index   = 0;
        $vj_tally = 0;

        while ( $index < $dictionary{$_}[2] ) {
            if (
                $dictionary{$_}[ $index * 2 + 6 ] =~ m|
                   V[[:alpha:]]?(\d+).*\;.*[jJ][[:alpha:]]?(\d+)
                 |x
              )
            {
                $v_fam = $1;
                $j_fam = $2;
            }
            else { die "Error interpreting $dictionary{$_}[$index*2+6].\n"; }

            if (    ( $v_fam == $dictionary{$_}[4] )
                and ( $j_fam == $dictionary{$_}[5] ) )
            {
                $vj_tally += $dictionary{$_}[ $index * 2 + 7 ];
            }
            $index++;
        }

        if ( $vj_tally <= ( $dictionary{$_}[3] / 2 ) ) {
            $dictionary{$_}[4] = 'ambig';
            $dictionary{$_}[5] = 'ambig';
        }
    }

    $index = 0;
    while ( $index < $dictionary{$_}[2] ) {
        if (
            $dictionary{$_}[ $index * 2 + 6 ] =~ m|
            V[[:alpha:]]?(\d+).*\;.*[jJ][[:alpha:]]?(\d+)
            |x
          )
        {
            $v_fam = $1;
            $j_fam = $2;
        }
        else { die "Error interpreting $dictionary{$_}[$index*2+6].\n"; }

        if (   ( $dictionary{$_}[4] =~ m|ambig| )
            or ( $v_fam != $dictionary{$_}[4] )
            or ( $j_fam != $dictionary{$_}[5] ) )
        {
            $dictionary{$_}[ $index * 2 + 6 ] .= ';chimera';
            $dictionary{$_}[3] -= $dictionary{$_}[ $index * 2 + 7 ]
              ;    # subtract the chimeric sequence counts from the total
        }

        $index++;
    }

    $dictionary{$_}[4] = 'V' . $dictionary{$_}[4];
    $dictionary{$_}[5] = 'J' . $dictionary{$_}[5];

}

print "# These collapsed down to ", scalar @keys, " rearrangements.\n";
print "# Counts:\n";

@keys = sort { $dictionary{$b}[3] <=> $dictionary{$a}[3] } ( keys %dictionary );

foreach (@keys) {
    print join( "\t", @{ $dictionary{$_} } ), "\n";
}

close $src_file;
