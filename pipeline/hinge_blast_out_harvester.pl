#!/usr/bin/perl -w
#
# harvest hinge-blast-result data and output an annotated fasta file

use strict;
use warnings;

local $/ = "\n";

sub verify_ids # check whether the query id matched the header of the target sequence
               # 1st argument -- query id
               # 2nd argument -- header
               # returns nothing
{
    my ( $test_query_id, $test_header ) = @_;    # unpacking function arguments
    my $sequence_id;

    if ( !( defined $test_query_id and defined $test_header ) ) {
        die "Code error in the verification subroutine.\n";
    }

    if ( $test_header =~ m|^>([\S]*)|x ) {
        $sequence_id = $1;
        if ( $test_query_id ne $sequence_id ) {
            die
              "Error: mismatched entries ...\n$test_query_id\n$test_header\n";
        }
        else { undef $sequence_id; }
    }
    else { die "Error parsing: $test_header"; }
    return;
}

sub order_hits    # order the hits loaded into hash according to top bitscore
                  #  and return the string to be appended as "subclass"
                  # 1st argument -- hash of hits/bitscores
                  # returns the top hits string composed of allele names
{
    my (%hits_bits) = @_;           # unpacking function arguments
    my @hit_ids;
    my @top_hit_ids;
    my $top_hits_string;

    if ( scalar %hits_bits ) {

# sort hits in reverse numerical order according to bit score, keeping to ASCII order for ties
        @hit_ids = sort { $hits_bits{$b} <=> $hits_bits{$a} or $b cmp $a }
          ( keys %hits_bits );

        foreach (@hit_ids) {
            if ( scalar @top_hit_ids ) {
                if ( $hits_bits{$_} == $hits_bits{ $top_hit_ids[0] } ) {
                    push @top_hit_ids, $_;
                }
            }
            else { push @top_hit_ids, $hit_ids[0]; }
        }

        # harvest just the allele information from the assignments
        foreach (@top_hit_ids) { s/.*\|(IGHG[^\s\|]*)\|.*/$1/x; }
        $top_hits_string = join( ',', @top_hit_ids );
    }
    else { $top_hits_string = "0null0"; }
    return $top_hits_string;
}

my $filename1;
my $filename2;

# process commandline
if ( scalar @ARGV == 2 ) {
    $filename1 = $ARGV[0];    # result.blast_out
    $filename2 = $ARGV[1];    # source.fasta
}
else {
    die
"Usage: perl hinge_blast_out_harvester.pl result.blast_out source.fasta\n";
}

open my $infile1, '<:encoding(UTF-8)', "$filename1"    # blast_out file
  or die "Cannot open $filename1.\n";
open my $infile2, '<:encoding(UTF-8)', "$filename2"    # input fasta
  or die "Cannot open $filename2.\n";

my $line;     # reading from the input blast_out file
my $line2;    # reading from the input fasta file
my $end_section = 0;
my $expected_format = "# Fields: query id, subject id, % identity, query length, subject length, alignment length, % query coverage per subject, bit score, evalue";
my $known_format;
my @query;

my %hits_bitscores;

$line = <$infile1>;    # eat the first "BLASTN"-labeled line to initialize

while ( defined( $line = <$infile1> ) ) {
    chomp $line;
    # lines marked as comments need special treatment
    if ( $line =~ m|^\#\s|x ) {
        # catch the start of the next section and set a flag
        if    ( $line =~ m|^\#\sBLAST[N\s]|x ) { $end_section = 1; }
        elsif ( $line =~ m|^\#\sQuery:\s([\S]+)|x ) {
            $query[0] = $1;
        }    # assign id to the section
        elsif ( $line eq $expected_format ) { $known_format = 1; }

        if    ($end_section)    # the main processing section
        {
            if ( defined( $line2 = <$infile2> ) ) {
                chomp $line2;
                verify_ids( $query[0], $line2 );
                print $line2, "\tsubclass:", order_hits(%hits_bitscores), "\n";
            }
            else { die "Error reading from the input FASTA." }

            # print the corresponding sequence
            if ( defined( $line2 = <$infile2> ) ) { print $line2; }
            else { die "Error reading from the input FASTA." }

            $end_section  = 0;        # starting new section of hits
            $known_format = 0;        # don't assume known format
            undef %hits_bitscores;    # clear the hits hash
            undef @query;             # clear the main query container array
        }
    }
    # data lines need to get parsed according to the predefined BLAST outfmt 7 set out in $expected_format
    else {
        if ( $known_format
            and $line =~ m|^([\S]+)\t # field1
                           ([\S]+)\t # field2
                           ([\S]+)\t # field3
                           ([\S]+)\t # field4
                           ([\S]+)\t # field5
                           ([\S]+)\t # field6
                           ([\S]+)\t # field7
                           ([\S]+)\t # field8
                           ([\S]+)   # field9
                        |x
          )
        {
            $query[1] = $1;    # query id
            $query[2] = $2;    # subject id
            $query[3] = $3;    # percent identity
            $query[4] = $4;    # query length
            $query[5] = $5;    # subject length
            $query[6] = $6;    # align length
            $query[7] = $7;    # percent coverage per subject
            $query[8] = $8;    # bitscore   <== these are being tracked
            $query[9] = $9;    # e-value

            if ( $query[0] ne $query[1] ) {
                die "\"$query[0]\" does not match \"$1\".\n";
            }

        }
        else { die "Encountered unknown blast output format. \n$line\n"; }

        if (
            !(
                defined $hits_bitscores{ $query[2] }
                and ( $query[8] <= $hits_bitscores{ $query[2] } )
            )
          )
        {
            $hits_bitscores{ $query[2] } = $query[8];
        }
    }

}

close $infile1;
close $infile2;
