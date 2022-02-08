#!/usr/bin/perl -w
#
# harvest igblast-result data and output an annotated fasta file

use strict;
use warnings;

local $/ = "\n";

sub codon2aa #adapted from http://www.techcuriosity.com/resources/bioinformatics/dna2protein.php
{
    my ($codon) = @_;
    $codon = uc $codon;
    my %gcode = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','AGC'=>'S','AGT'=>'S',
                 'TTC'=>'F','TTT'=>'F',
                 'TTA'=>'L','TTG'=>'L','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
                 'TAC'=>'Y','TAT'=>'Y',
                 'TAA'=>'*','TAG'=>'*','TGA'=>'*',
                 'TGC'=>'C','TGT'=>'C',
                 'TGG'=>'W',
                 'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
                 'CAC'=>'H','CAT'=>'H',
                 'CAA'=>'Q','CAG'=>'Q',
                 'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','AGA'=>'R','AGG'=>'R',
                 'ATA'=>'I','ATC'=>'I','ATT'=>'I',
                 'ATG'=>'M',
                 'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
                 'AAC'=>'N','AAT'=>'N',
                 'AAA'=>'K','AAG'=>'K',
                 'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
                 'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
                 'GAC'=>'D','GAT'=>'D',
                 'GAA'=>'E','GAG'=>'E',
                 'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G',
                 'NAA'=>'X','NAT'=>'X','NAG'=>'X','NAC'=>'X','NTA'=>'X','NTT'=>'X','NTG'=>'X','NTC'=>'X',
                 'NGA'=>'X','NGT'=>'X','NGG'=>'X','NGC'=>'X','NCA'=>'X','NCT'=>'X','NCG'=>'X','NCC'=>'X',
                 'ANA'=>'X','ANT'=>'X','ANG'=>'X','ANC'=>'X','TNA'=>'X','TNT'=>'X','TNG'=>'X','TNC'=>'X',
                 'GNA'=>'X','GNT'=>'X','GNG'=>'X','GNC'=>'X','CNA'=>'X','CNT'=>'X','CNG'=>'X','CNC'=>'X',
                 'AAN'=>'X','ATN'=>'X','AGN'=>'X','ACN'=>'X','TAN'=>'X','TTN'=>'X','TGN'=>'X','TCN'=>'X',
                 'GAN'=>'X','GTN'=>'X','GGN'=>'X','GCN'=>'X','CAN'=>'X','CTN'=>'X','CGN'=>'X','CCN'=>'X',
                 'ANN'=>'X','TNN'=>'X','CNN'=>'X','GNN'=>'X',
                 'NAN'=>'X','NTN'=>'X','NCN'=>'X','NGN'=>'X',
                 'NNA'=>'X','NNT'=>'X','NNC'=>'X','NNG'=>'X',
                 'NNN'=>'X'
    );

    if ( exists $gcode{$codon} ) {
        return $gcode{$codon};
    }
    else {
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}

sub translate #adapted from http://www.techcuriosity.com/resources/bioinformatics/dna2protein.php
{
    my ($nuc) = @_;
    my $aa    = '';
    my $codon = '';
    my $index = 0;

    for ( $index = 0 ; $index < ( length($nuc) - 2 ) ; $index += 3 ) {
        $codon = substr( $nuc, $index, 3 );
        $aa .= codon2aa($codon);
    }

    return $aa;
}

sub rev_comp    # determine the reverse-complement sequence
                # 1st argument -- input sequence
                # return reverse-complement

{
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return $seq;
}

sub parse_igblast_block   # read from the IgBLAST output file and store data
                          # 1st argument -- filehandle for reading inputFile
                          # 2nd argument -- first line in the block ("Query...")
                          # returns hash containing (keys):
                          #   ['query']
                          #   ['q_length']
                          #   ['gene_usage']
                          #   ['rearr']
                          #   ['cdr3_nt']
                          #   ['cdr3_aa']
                          #   ['perc_ident']
                          #   ['cov']
                          #   ['trunc_flags'] array
                          #   ['fwk_bounds']  array
                          #   ['q_rev_flag']

{
    my ( $fh, $line ) = @_;

    my $end_table_flag = 0;
    my $end_block_flag = 0;
    my $aa_set = "[ACDEFGHIKLMNPQRSTVWXY\*]";
    my $nt_set = "[ACGTN]";
    my %alignment_tbl_data;
    my %result;

    $result{'q_length'}   = 0;
    $result{'gene_usage'} = '';
    $result{'rearr'}      = '';
    $result{'cdr3_aa'}    = '0null0';
    $result{'cdr3_nt'}    = '0null0';
    $result{'perc_ident'} = 0;
    $result{'cov'}        = 0;

    # positions for Fr1, CDR1, Fr2, CDR2, Fr3, CDR3, J-family assignment
    $result{'trunc_flags'} = [ 1, 1, 1, 1, 1, 1, 1 ];
    $result{'fwk_bounds'}  = [];
    $result{'q_rev_flag'}  = 0;

    $line =~ m|^Query=\s(\S*)|xs
      or die "Error: invald start of IgBLAST output section at $line.";
    $result{'query'} = $1;

    while ( !$end_block_flag && defined( $line = <$fh> ) ) {
        for ($line) {
            if (/^Length=(\d+)/x) { $result{'q_length'} = $1; }
            elsif (/Effective search space used/) { $end_block_flag = 1; }
            elsif (/rearrangement summary for query sequence/) {
                $line = <$fh>;    # grab the next line for processing
                chomp $line;
                $result{'gene_usage'} = $line;
                if ( $line =~ m|\t\-$|x ) { $result{'q_rev_flag'} = 1; }
            }
            elsif (/junction details based on top germline gene/) {
                $line = <$fh>;    # grab the next line for processing
                chomp $line;
                $result{'rearr'} = $line;

                # the last field denotes presence of J chain
                if ( $result{'rearr'} !~ m|\tN\/A\s*$|x ) {
                    $result{'trunc_flags'}[6] = 0;
                }

                $result{'rearr'} =~ s|N\/A||gx;

                # This may indicate overlapping sequences; proceed with caution.
                $result{'rearr'} =~ s|[\(,\)]||gx;
                $result{'rearr'} =~ s|\t||gx;
            }
            elsif (/region sequence details/) {
                $line = <$fh>;    # grab the next line for processing
                $line =~ m|^CDR3\s+($nt_set+)\s+($aa_set*)\s.*|x
                  or die
                    "Error: invalid IgBLAST output format (CDR3 region); check data for $result{'query'}\n";
                $result{'cdr3_nt'} = $1;
                $result{'cdr3_aa'} = length($2) ? $2 : "0null0";
            }
            elsif (/^Alignment\ssummary/x)
            # Framework alignment summary table; expected to end with the "Total" line
            {
                while ( !$end_table_flag && defined( $line = <$fh> ) ) {
                    # Expected fields (from, to, length, matches, mismatches, gaps, percent identity)
                    if ($line =~ /^([^\t]+)\t(\S+)\t(\S+)\t(\d+)\t\S+\t\S+\t\S+\t(\d+\.?\d*)/x) {
                        $alignment_tbl_data{'rowname'} = $1;
                        $alignment_tbl_data{'from'} = $2;
                        $alignment_tbl_data{'to'} = $3;
                        $alignment_tbl_data{'length'} = $4;
                        $alignment_tbl_data{'pct_id'} = $5;

                        for ($alignment_tbl_data{'rowname'}) {
                            if    (/FR1/) { $result{'trunc_flags'}[0] = 0; }
                            elsif (/CDR1/) {
                              $result{'trunc_flags'}[1] = 0;
                            }
                            elsif (/FR2/) { $result{'trunc_flags'}[2] = 0; }
                            elsif (/CDR2/) {
                              $result{'trunc_flags'}[3] = 0;
                            }
                            elsif (/FR3/) { $result{'trunc_flags'}[4] = 0; }
                            elsif (/CDR3/) {
                              $result{'trunc_flags'}[5] = 0;
                            }
                            elsif (/Total/) {
                              $result{'perc_ident'} = $alignment_tbl_data{'pct_id'};
                              $result{'cov'}        = sprintf( "%.1f",
                                  $alignment_tbl_data{'length'} / $result{'q_length'} * 100 );
                              $end_table_flag       = 1;
                            }
                            else {
                              die
                                  "Error interpreting alignment at $result{'query'} summary.";
                            }
                        }
                        # collect framework bounds for reading frame determination
                        if ( $alignment_tbl_data{'from'} ne 'N/A' ) {
                            push @{ $result{'fwk_bounds'} }, $alignment_tbl_data{'from'};
                            push @{ $result{'fwk_bounds'} }, $alignment_tbl_data{'to'};
                        }
                    }
                }
            }
            elsif (/\Q***** No hits found *****\E/x) {
                $result{'gene_usage'} = "invalid_query_seq";
                $result{'rearr'}      = '';
                $result{'cov'}        = 0;
                $result{'perc_ident'} = 0;
            }
        }
    }
    return %result;
}

#### main section
my $filename1;
my $filename2;

# process commandline
if ( scalar @ARGV == 2 ) {
    $filename1 = $ARGV[0];
    $filename2 = $ARGV[1];
}
else {
    die
      "Usage: perl igblast-out_harvester.pl source.igblast_out source.fastq\n";
}

my $line;

my $query_id;
my $query_seq;
my $working_seq;
my $readframe;
my %igblast_data;
my $rearr_aa;
my $translation;
my $idx;

open my $infile1, '<:encoding(UTF-8)', "$filename1"    # igblast_out file
  or die "Cannot open $filename1.\n";
open my $infile2, '<:encoding(UTF-8)',, "$filename2"    # input fasta
  or die "Cannot open $filename1.\n";

while ( defined( $line = <$infile1> ) ) {
    if ( $line =~ m|^Query=\s|xs ) {

        # initialization block
        $rearr_aa    = '0null0';
        $readframe   = 0;
        $translation = '';         #clear the translation string

        %igblast_data = parse_igblast_block( $infile1, $line );

        # grab the sequence from the FASTA file
        if ( defined( $line = <$infile2> )
            && ( $line =~ m|^>$igblast_data{'query'}.*|x ) )
        {
            $query_id = $line; # if the line is long, IgBLAST truncates the header; store the real one
            chomp $query_id;
            $query_id =~ s|^>||x;

            if ( !defined( $query_seq = <$infile2> ) ) {
                die "Error: cannot retrieve\n\t$query_id\n\tfrom $filename2\n";
            }
            chomp $query_seq;
        }
        else { die "Error at $igblast_data{'query'} FASTA entry retrieval."; }

        # resolve the reading frame from the framework boundaries
        # based on IgBLAST-determined framework boundaries
        if (
               ( defined $igblast_data{'fwk_bounds'}[2] )
            && ( defined $igblast_data{'fwk_bounds'}[1] )
            && (
                ( $igblast_data{'fwk_bounds'}[2] -
                        $igblast_data{'fwk_bounds'}[1] ) == 1 )
               )
        {
            $readframe = ( $igblast_data{'fwk_bounds'}[1] ) % 3 + 1;
        }
        else { $readframe = 0; }

        # if the sequence is determined to be "reversed", determine the revcomp
        if ( $igblast_data{'q_rev_flag'} ) { $working_seq = rev_comp($query_seq); }
        else { $working_seq = $query_seq; }

        # translate the sequence, possibly reverse-complement
        if ($readframe) {
            $translation =
              translate( substr( $working_seq, ( $readframe - 1 ) ) );

            # fix the translation for cases when only a portion of the sequence
            #   was used in the igblast annotation
            if (  !( $igblast_data{'cdr3_aa'} eq '0null0' )
                && ( $translation !~ m|\Q$igblast_data{'cdr3_aa'}\E|x ) )
            {
                $readframe = 0; # this will be reset based on the correct translation
                while (( $readframe < 3 )
                    && ( $translation !~ m|\Q$igblast_data{'cdr3_aa'}\E|x ) )
                {
                    $translation =
                      translate( substr( $working_seq, ( $readframe++ ) ) );
                }
            }
        }

       # obtain context residues for the cdr3_aa:
       #    \Q and \E ensure that special characters are escaped (see quotemeta)
        if ( !( $igblast_data{'cdr3_aa'} eq '0null0' ) &&
              ( $translation =~ m| [\w\*]+
                                  ([\w\*]{3}\Q$igblast_data{'cdr3_aa'}\E[\w\*]{2})
                                   [\w\*]+
                                 |x )
           )
        {
              $igblast_data{'cdr3_aa'} = $1;
        }
        else { $igblast_data{'cdr3_aa'} = '0null0'; }

        # generate the translation for the junction sequence TODO: disagreements with old code
        if ( length( $igblast_data{'rearr'} ) && !( $igblast_data{'cdr3_aa'} eq '0null0' ) ) {
            $idx = 0;
            while (( $idx < 3 )
                && ( $igblast_data{'cdr3_aa'} !~ m|\Q$rearr_aa\E|x ) )
            {
                $rearr_aa = translate( substr( $igblast_data{'rearr'}, $idx ) );
                $idx++;    # increment frame index
            }

            if ( $translation !~ m|\Q$rearr_aa\E|x ) { $rearr_aa = "0null0"; }
        }

        ### the output block starts here
        print "\>$query_id\t$igblast_data{'gene_usage'}\t";

        if (   $igblast_data{'trunc_flags'}[0]
            or $igblast_data{'trunc_flags'}[1]
            or $igblast_data{'trunc_flags'}[2]
            or $igblast_data{'trunc_flags'}[3]
            or $igblast_data{'trunc_flags'}[4] )
        {
            print "Vtruncated.", @{ $igblast_data{'trunc_flags'} }, "\t";
        }
        elsif ( $igblast_data{'trunc_flags'}[6] ) {
            print "Jtruncated.", @{ $igblast_data{'trunc_flags'} }, "\t";
        }
        else { print "Vintact\t"; }

        if ( $igblast_data{'q_rev_flag'} ) { $readframe = -$readframe; }

        print "junctnn:$igblast_data{'rearr'}\t";
        print "junctaa:$rearr_aa\t";

        # print   "CDR3nt:$igblast_data{'cdr3_nt'}\t";
        print "CDR3aa:$igblast_data{'cdr3_aa'}\t";
        print "frame:$readframe\t";
        print "pcov:$igblast_data{'cov'}\t";
        print "pid:$igblast_data{'perc_ident'}\t";
        print "transl:$translation\n";
        print "$query_seq\n";
        ### end output block
    }
}

# handle premature end of IgBLAST processing output
if ( defined( $line = <$infile2> ) ){
  close $infile1;

  print STDERR "!!! Warning !!! There are more sequenced in the input file than have been annotated.\n";
  print STDERR "Processing stopped at the following id:\n";
  print STDERR $line;

  my $outfilename = time() . ".unprocessed.fasta";
  open my $outfile, '>:encoding(UTF-8)', "$outfilename" # unprocessed FASTA reamainder
    or die "Cannot open $outfilename.\n";

  print $outfile $line;
  while ( $line = <$infile2> ){
    print $outfile $line;
  }

  close $outfile;
  close $infile2;
  die "!!! Warning !!! IgBLAST output processing has stopped prematurely.\n";
}

close $infile1;
close $infile2;
