"""
igblast-out_harvester.py
  harvest igblast-result data and output an annotated fasta file
"""

import sys
import argparse
import re
from os.path import exists


def codon2aa(codon):
    '''
        Look up the amino acid from a codon.
        1st argument--codon
        returns an amino acid
        adapted from http://www.techcuriosity.com/resources/bioinformatics/dna2protein.php
    '''

    gcode = {'TCA':'S','TCC':'S','TCG':'S','TCT':'S','AGC':'S','AGT':'S',
             'TTC':'F','TTT':'F',
             'TTA':'L','TTG':'L','CTA':'L','CTC':'L','CTG':'L','CTT':'L',
             'TAC':'Y','TAT':'Y',
             'TAA':'*','TAG':'*','TGA':'*',
             'TGC':'C','TGT':'C',
             'TGG':'W',
             'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
             'CAC':'H','CAT':'H',
             'CAA':'Q','CAG':'Q',
             'CGA':'R','CGC':'R','CGG':'R','CGT':'R','AGA':'R','AGG':'R',
             'ATA':'I','ATC':'I','ATT':'I',
             'ATG':'M',
             'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
             'AAC':'N','AAT':'N',
             'AAA':'K','AAG':'K',
             'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
             'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
             'GAC':'D','GAT':'D',
             'GAA':'E','GAG':'E',
             'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
             'NAA':'X','NAT':'X','NAG':'X','NAC':'X','NTA':'X','NTT':'X','NTG':'X','NTC':'X',
             'NGA':'X','NGT':'X','NGG':'X','NGC':'X','NCA':'X','NCT':'X','NCG':'X','NCC':'X',
             'ANA':'X','ANT':'X','ANG':'X','ANC':'X','TNA':'X','TNT':'X','TNG':'X','TNC':'X',
             'GNA':'X','GNT':'X','GNG':'X','GNC':'X','CNA':'X','CNT':'X','CNG':'X','CNC':'X',
             'AAN':'X','ATN':'X','AGN':'X','ACN':'X','TAN':'X','TTN':'X','TGN':'X','TCN':'X',
             'GAN':'X','GTN':'X','GGN':'X','GCN':'X','CAN':'X','CTN':'X','CGN':'X','CCN':'X',
             'ANN':'X','TNN':'X','CNN':'X','GNN':'X',
             'NAN':'X','NTN':'X','NCN':'X','NGN':'X',
             'NNA':'X','NNT':'X','NNC':'X','NNG':'X',
             'NNN':'X'
            }

    if codon in gcode:
        return gcode[codon]

    sys.exit('Bad codon \"' + codon + '\"!!\n\n')

def translate(seq):
    '''
    translate input sequence
      adapted from http://www.techcuriosity.com/resources/bioinformatics/dna2protein.php
    '''

    codon = ''
    transl = ''
    ind = 0

    while ind < (len(seq) - 2):
        codon = seq[ind:ind+3]
        transl = transl + codon2aa(codon)
        ind = ind + 3

    return transl

def rev_comp(seq):
    '''
    returns reverse-complement
    '''
    seq = seq[::-1]
    complement = {"A":"T", "T":"A",
                  "C":"G", "G":"C",
                  "a":"t", "t":"a",
                  "c":"g", "g":"c"}
    table = seq.maketrans(complement)
    return seq.translate(table)

def parse_igblast_block(file, line):
    '''
    read from the IgBLAST output file and store data
    1st argument -- filehandle for reading inputFile
    2nd argument -- first line in the block ("Query...")
    returns hash containing (keys):
      ['query']
      ['q_length']
      ['gene_usage']
      ['rearr']
      ['cdr3_nt']      ... CDR3 nucleotide sequence
      ['cdr3_bounds']  array
      ['cdr3_aa']      ... CDR3 amino acid sequence
      ['perc_ident']
      ['cov']
      ['trunc_flags'] array
      ['fwk_bounds']  array
      ['q_rev_flag']
    '''

    end_table_flag        = 0
    end_block_flag        = 0
    aa_set                = r'[ACDEFGHIKLMNPQRSTVWXY\*]'
    nt_set                = r'[ACGTN]'
    alignment_tbl_data    = {}
    result                = {}

    result['q_length']    = 0
    result['gene_usage']  = ''
    result['rearr']       = ''
    result['cdr3_aa']     = '0null0'
    result['cdr3_nt']     = '0null0'
    result['perc_ident']  = 0
    result['cov']         = 0

    # positions for Fr1, CDR1, Fr2, CDR2, Fr3, CDR3, J-family assignment
    result['trunc_flags'] = [1, 1, 1, 1, 1, 1, 1 ]
    result['cdr3_bounds'] = [0, 0]
    result['fwk_bounds']  = []
    result['q_rev_flag']  = 0

    match_result = re.search(r'^Query=\s(\S*)', line)

    if match_result is None:
        sys.exit('Error: invald start of IgBLAST output section at ' + line + '.')

    result['query'] = match_result.group(1)

    line = file.readline() # initial read

    while not end_block_flag and line:
        # skip all empty lines
        while re.search(r'^\s*\n', line):
            line = file.readline()

        match_result = re.search(r'Length=(\d+)', line)
        if match_result:
            result['q_length'] = int(match_result.group(1))

            line = file.readline()
            # skip all empty lines
            while re.search(r'^\s*\n', line):
                line = file.readline()

        match_result = re.search(r'rearrangement summary for query sequence', line)
        if match_result:
            line = file.readline()
            result['gene_usage'] = line.strip()
            if re.search(r'\t\-\t?[^\t]*$', line):
                result['q_rev_flag'] = 1

            line = file.readline()
            # skip all empty lines
            while re.search(r'^\s*\n', line):
                line = file.readline()

        match_result = re.search(r'junction details based on top germline gene', line)
        if match_result:
            line = file.readline()
            # the last field denotes presence of J chain
            if re.search(r'\tN\/A\s*$', line) is None:
                result['trunc_flags'][6] = 0

            # This may indicate overlapping sequences; proceed with caution.
            line = re.sub(r'N\/A|\(|\)|\t', '', line)
            result['rearr'] = line.strip()

            line = file.readline()
            # skip all empty lines
            while re.search(r'^\s*\n', line):
                line = file.readline()

        match_result = re.search(r'region sequence details', line)
        if match_result:
            line = file.readline()
            match_result = re.search(r'^CDR3\s+(' + nt_set \
                            + r'+)\s+(' + aa_set \
                            + r'*)\t(\d+)\t(\d+)', line)
            if match_result is None:
                sys.exit('Error: invalid IgBLAST output format (CDR3 region); check data for '\
                            + result['query'])

            result['cdr3_nt'] = match_result.group(1)
            result['cdr3_aa']  = match_result.group(2) if match_result.group(2) else '0null0'
            result['cdr3_bounds'][0] = int(match_result.group(3))
            result['cdr3_bounds'][1] = int(match_result.group(4))

            line = file.readline()
            # skip all empty lines
            while re.search(r'^\s*\n', line):
                line = file.readline()

        match_result = re.search(r'^Alignment summary', line)
        # Framework alignment summary table; expected to end with the "Total" line
        if match_result:
            while not end_table_flag and line:
                line = file.readline()
                match_result = re.search(r'^([^\t]+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t(\S+)', line)
                if match_result:
                    alignment_tbl_data['rowname'] = match_result.group(1)
                    alignment_tbl_data['from'] = match_result.group(2)
                    alignment_tbl_data['to'] = match_result.group(3)

                    if re.match(r'\d+',match_result.group(4)):
                        alignment_tbl_data['length'] = int(match_result.group(4))
                    else:
                        alignment_tbl_data['length'] = 0

                    if re.match(r'\d+\.?\d*',match_result.group(5)):
                        alignment_tbl_data['pct_id'] = float(match_result.group(5))
                    else:
                        alignment_tbl_data['pct_id'] = 0

                    if re.search(r'^FR1', alignment_tbl_data['rowname']):
                        result['trunc_flags'][0] = 0
                    elif re.search(r'^CDR1', alignment_tbl_data['rowname']):
                        result['trunc_flags'][1] = 0
                    elif re.search(r'^FR2', alignment_tbl_data['rowname']):
                        result['trunc_flags'][2] = 0
                    elif re.search(r'^CDR2', alignment_tbl_data['rowname']):
                        result['trunc_flags'][3] = 0
                    elif re.search(r'^FR3', alignment_tbl_data['rowname']):
                        result['trunc_flags'][4] = 0
                    elif re.search(r'^CDR3', alignment_tbl_data['rowname']):
                        result['trunc_flags'][5] = 0
                    elif alignment_tbl_data['rowname'] == 'Total':
                        result['perc_ident'] = alignment_tbl_data['pct_id']
                        if result['q_length']:
                            result['cov'] = 100 * alignment_tbl_data['length'] / result['q_length']
                        else:
                            result['cov'] = 0
                        end_table_flag = 1
                    else:
                        sys.exit('Error interpreting alignment at ' + result['query'] + ' summary.')

                # collect framework bounds for reading frame determination
                if alignment_tbl_data['rowname'] != 'Total':
                    result['fwk_bounds'].append(int(alignment_tbl_data['from']))
                    result['fwk_bounds'].append(int(alignment_tbl_data['to']))

        match_result = re.search(re.escape(r'***** No hits found *****'), line)
        if match_result:
            result['gene_usage'] = 'invalid_query_seq'
            result['rearr']      = ''
            result['cov']        = 0
            result['perc_ident'] = 0

        match_result = re.search(r'Effective search space used', line)
        if match_result:
            end_block_flag = 1

        line = file.readline()

    return result

def compose_fasta_block(file, filename, line, igblast_data_dict):
    '''
    read a FASTA entry, determine the reading frame from igblast_data_dict,
        compose the new description line, and return the FASTA block
    1st argument -- file handle for input FASTA file
    2nd argument -- filename for input FASTA file
    3rd argument -- first line of input FASTA file
    4th argument -- igblast_data_dict
    '''
    # initialize
    result      = {}
    readframe   = 0
    rearr_aa    = '0null0'
    translation = ''
    aa_set      = r'[ACDEFGHIKLMNPQRSTVWXY\*]'

    if re.search(r'^>' + igblast_data_dict['query'] + '.*', line) is None:
        sys.exit('Error at ' + igblast_data_dict['query'] + ' FASTA entry retrieval.')

    result['query_id'] = line.strip()
    result['query_id'] = re.sub(r'^>', '', result['query_id'])

    # grab the sequence from the FASTA file
    line = file.readline()
    if line:
        result['query_seq'] = line.strip()
    else:
        sys.exit('Error: cannot retrieve the sequence for\n\t' + \
          result['query_id'] + '\n\tfrom ' + filename)

    if len(igblast_data_dict['fwk_bounds']) > 1:
        readframe = igblast_data_dict['fwk_bounds'][1] % 3 + 1

    # another possibility that isn't always available:
    # readframe = str((result['cdr3_bounds'][0]-1)%3+1)

    # if the sequence is determined to be "reversed", determine the revcomp
    working_seq = rev_comp(result['query_seq']) \
        if igblast_data_dict['q_rev_flag'] else result['query_seq']

    # translate the sequence, possibly reverse-complement
    if readframe:
        translation = translate(working_seq[readframe - 1:])

        # fix the translation for cases when only a portion of the sequence
        #   was used in the igblast annotation
        if igblast_data_dict['cdr3_aa'] != '0null0' \
          and re.search(re.escape(igblast_data_dict['cdr3_aa']), translation) is None:
            readframe = 3
            while readframe >= 1:
                translation = translate(working_seq[readframe - 1 :])
                if re.search(re.escape(igblast_data_dict['cdr3_aa']), translation):
                    break
                readframe = readframe - 1

    # obtain context residues for the cdr3_aa
    if igblast_data_dict['cdr3_aa'] != '0null0':
        match_result = re.search(aa_set + r'+(' + aa_set + r'{3}' + \
                       re.escape(igblast_data_dict['cdr3_aa']) + \
                       aa_set + r'{2})', translation)
        if match_result:
            igblast_data_dict['cdr3_aa'] = match_result.group(1)
        else:
            igblast_data_dict['cdr3_aa'] = '0null0'  # TODO, see below
            # sys.exit('Error: cannot obtain CDR3 amino acid sequence context.')

    if igblast_data_dict['q_rev_flag']:
        readframe = -readframe

    # generate the translation for the junction sequence
    if len(igblast_data_dict['rearr']) and igblast_data_dict['cdr3_aa'] != '0null0':
        for idx in range(0,3):
            rearr_aa = translate(igblast_data_dict['rearr'][idx:])
            if re.search(re.escape(rearr_aa), igblast_data_dict['cdr3_aa']):
                break

        if re.search(re.escape(rearr_aa), translation) is None:
            rearr_aa = '0null0'

    ### result construction start
    result['query_id'] = '>' + result['query_id'] + '\t' \
                             + igblast_data_dict['gene_usage'] + '\t'
    if   igblast_data_dict['trunc_flags'][0] \
      or igblast_data_dict['trunc_flags'][1] \
      or igblast_data_dict['trunc_flags'][2] \
      or igblast_data_dict['trunc_flags'][3] \
      or igblast_data_dict['trunc_flags'][4]:
        result['query_id'] = result['query_id'] + 'Vtruncated.' + \
          ''.join(map(str,igblast_data_dict['trunc_flags']))
    elif igblast_data_dict['trunc_flags'][6]:
        result['query_id'] = result['query_id'] + 'Jtruncated.' + \
          ''.join(map(str,igblast_data_dict['trunc_flags']))
    else:
        result['query_id'] = result['query_id'] + 'Vintact'

    result['query_id'] = result['query_id'] + '\t' + \
      'junctnn:' + igblast_data_dict['rearr'] + '\t' + \
      'junctaa:' + rearr_aa + '\t' + \
      'CDR3aa:' + igblast_data_dict['cdr3_aa'] + '\t' + \
      'frame:' + str(readframe) + '\t' + \
      'pcov:' + f"{igblast_data_dict['cov']:.1f}" + '\t' + \
      'pid:' + f"{igblast_data_dict['perc_ident']:.1f}" + '\t' + \
      'transl:' + translation

    return result

#-------------------------------------------------------------------------------
#### main section
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('igblastOut_name', \
        help='Filename for the IgBLAST output (e.g., "source.igblast_out")')
    parser.add_argument('fasta_name', \
        help='Filename for the FASTA data set (e.g., "source.fasta")')
    #parser.add_argument('--debug', help='output debug information', action='store_true')
    args = parser.parse_args()

try:
    with open(args.igblastOut_name, encoding="utf8") as igblast, \
         open(args.fasta_name, encoding="utf8") as fasta:
        in_line = igblast.readline()

        while in_line:
            # print("0")
            if re.search(r'^Query=\s', in_line):
                igblast_data = parse_igblast_block(igblast, in_line)
                in_line = fasta.readline()
                annotated_fasta = compose_fasta_block(fasta, \
                  args.fasta_name, in_line, igblast_data)
                print(annotated_fasta['query_id'])
                print(annotated_fasta['query_seq'])
            in_line = igblast.readline()

except FileNotFoundError:
    if not exists(args.igblastOut_name):
        sys.exit('File ' + args.igblastOut_name + ' was not found!')
    else:
        sys.exit('File ' + args.fasta_name + ' was not found!')
