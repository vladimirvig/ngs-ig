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

def compose_fasta_block(keys, data):
    '''
    read a FASTA entry, determine the reading frame from igblast_data_dict,
        compose the new description line, and return the FASTA block
    1st argument -- header from the AIRR table
    2nd argument -- data line from the AIRR table\
    '''
    # initialize
    result      = {}
    readframe   = 0

    # positions for Fr1, CDR1, Fr2, CDR2, Fr3, CDR3, J-segment annotation (if any)
    trunc_flags = [1,1,1,1,1,1,1]
    fwr_ends = [0,0,0,0,0,0,0]
    
    translation = ''
    aa_set      = r'[ACDEFGHIKLMNPQRSTVWXY\*]'
    
    if len(keys) != len(data):
        sys.exit('Error reading the AIRR table due to field number mismatch at: ' + data[0])
        
    entry_dict = dict(zip(keys, data))
    
    # store the sequence for output
    result['query_seq'] = entry_dict['sequence']
    
    # compose the identifier line from annotations stored in entry_dict
    result['query_id'] = entry_dict['sequence_id'] + '\t'
    
    # V/D/J assignments
    if entry_dict['v_call'] != '':
        result['query_id'] = result['query_id'] + entry_dict['v_call'] + '\t'
    else:
        result['query_id'] = result['query_id'] + 'N/A' + '\t'
    if entry_dict['d_call'] != '':
        result['query_id'] = result['query_id'] + entry_dict['d_call'] + '\t'
    else:
        result['query_id'] = result['query_id'] + 'N/A' + '\t'
    if entry_dict['j_call'] != '':
        result['query_id'] = result['query_id'] + entry_dict['j_call'] + '\t'
    else:
        result['query_id'] = result['query_id'] + 'N/A' + '\t'
    if entry_dict['locus'] != '':
        result['query_id'] = result['query_id'] + entry_dict['locus'] + '\t'
    else:
        result['query_id'] = result['query_id'] + 'invalid_query_seq' + '\t'

    # miscelaneous annotations (stop codons, productive status, frame-shifts, etc.)
    if entry_dict['stop_codon'] == 'F':
        result['query_id'] = result['query_id'] + 'No' + '\t'
    elif entry_dict['stop_codon'] == 'T':
        result['query_id'] = result['query_id'] + 'Yes' + '\t'
    else:
        result['query_id'] = result['query_id'] + 'Unknown' + '\t'
    if entry_dict['vj_in_frame'] == 'F':
        result['query_id'] = result['query_id'] + 'Out-of-frame' + '\t'
    elif entry_dict['vj_in_frame'] == 'T':
        result['query_id'] = result['query_id'] + 'In-frame' + '\t'
    else:
        result['query_id'] = result['query_id'] + 'Unknown' + '\t'
    if entry_dict['productive'] == 'F':
        result['query_id'] = result['query_id'] + 'No' + '\t'
    elif entry_dict['productive'] == 'T':
        result['query_id'] = result['query_id'] + 'Yes' + '\t'
    else:
        result['query_id'] = result['query_id'] + 'Unknown' + '\t'
    if entry_dict['rev_comp'] == 'F':
        result['query_id'] = result['query_id'] + '+' + '\t'
    elif entry_dict['rev_comp'] == 'T':
        result['query_id'] = result['query_id'] + '-' + '\t'
    else:
        result['query_id'] = result['query_id'] + 'Unknown' + '\t'
    if entry_dict['v_frameshift'] == 'F':
        result['query_id'] = result['query_id'] + 'No' + '\t'
    elif entry_dict['v_frameshift'] == 'T':
        result['query_id'] = result['query_id'] + 'Yes' + '\t'
    else:
        result['query_id'] = result['query_id'] + 'Unknown' + '\t'

    # Determine which regions are present for the truncation annotation
    if entry_dict['fwr1'] != '':
        trunc_flags[0] = 0
        fwr_ends[0] = int(entry_dict['fwr1_end'])
    if entry_dict['cdr1'] != '':
        trunc_flags[1] = 0
        fwr_ends[1] = int(entry_dict['cdr1_end'])
    if entry_dict['fwr2'] != '':
        trunc_flags[2] = 0
        fwr_ends[2] = int(entry_dict['fwr2_end'])
    if entry_dict['cdr2'] != '':
        trunc_flags[3] = 0
        fwr_ends[3] = int(entry_dict['cdr2_end'])
    if entry_dict['fwr3'] != '':
        trunc_flags[4] = 0
        fwr_ends[4] = int(entry_dict['fwr3_end'])
    if entry_dict['cdr3'] != '':
        trunc_flags[5] = 0
        fwr_ends[5] = int(entry_dict['cdr3_end'])
    if entry_dict['fwr4'] != '':
        trunc_flags[6] = 0
        fwr_ends[6] = int(entry_dict['fwr4_end'])
    
    if   trunc_flags[0] \
      or trunc_flags[1] \
      or trunc_flags[2] \
      or trunc_flags[3] \
      or trunc_flags[4]:
        result['query_id'] = result['query_id'] + 'Vtruncated.' + \
          ''.join(map(str,trunc_flags)) + '\t'
    # elif trunc_flags[6]:
    #     result['query_id'] = result['query_id'] + 'Jtruncated.'  + '\t' \
    #       ''.join(map(str,trunc_flags))
    else:
        result['query_id'] = result['query_id'] + 'Vintact' + '\t'
    
    if entry_dict['junction'] != '':        
        result['query_id'] = result['query_id'] + 'junctnn:' + entry_dict['junction'] + '\t'
    else:
        result['query_id'] = result['query_id'] + 'junctnn:' + '0null0' + '\t'

    if entry_dict['junction_aa'] != '':
        result['query_id'] = result['query_id'] + 'junctaa:' + entry_dict['junction_aa'] + '\t'
    else:
        result['query_id'] = result['query_id'] + 'junctaa:' + '0null0' + '\t'
    
    # determine reading frame
    for value in fwr_ends:
        if value:
            # if readframe and readframe != value % 3 + 1:
            #     sys.exit('Conflict among readframes in :' + result['query_id'] + ' ' + ' '.join(map(str,fwr_ends)))
            readframe = value % 3 + 1
 
    # determine translation
    if readframe != 0:
        translation = translate(entry_dict['sequence'][readframe-1:])
    else:
        translation = '0null0'
    
    # obtain context residues for the cdr3_aa
    if entry_dict['cdr3_aa'] != '':
        match_result = re.search(aa_set + r'+(' + aa_set + r'{3}' + \
                       re.escape(entry_dict['cdr3_aa']) + \
                       aa_set + r'{2})', translation)
        if match_result:
            cdr3_aa = match_result.group(1)
        else:
            cdr3_aa = '0null0'
    else: 
        cdr3_aa = '0null0'
    result['query_id'] = result['query_id'] + 'CDR3aa:' + cdr3_aa + '\t'

    # declare reading frame
    result['query_id'] = result['query_id'] + 'frame:' + str(readframe) + '\t'
    
    # declare percent covered by the alignment
    if entry_dict['sequence_alignment'] != '':
        seq_aln = entry_dict['sequence_alignment'].replace('-','')
        pcov = 100 * len(seq_aln) / len(entry_dict['sequence'])
        result['query_id'] = result['query_id'] + 'pcov:' + str(round(pcov,1)) + '\t'
    else:
        result['query_id'] = result['query_id'] + 'pcov:' + '0null0' + '\t'
        
    # declare percent identity for the V-segment
    if entry_dict['v_identity'] != '':
        identity = round(float(entry_dict['v_identity']),1)
    else:
        identity = '0null0'
    result['query_id'] = result['query_id'] + 'pid:' + str(identity) + '\t'
 
    # declare the tranlation
    result['query_id'] = result['query_id'] + 'transl:' + translation

    ## resulting convention:
    # 1: 'sequence_id'
    # 2: 'v_call'
    # 3: 'd_call'
    # 4: 'j_call'
    # 5: 'locus'
    # 6: 'stop_codon'
    # 7: 'vj_in_frame'
    # 8: 'productive'
    # 9: 'rev_comp'
    # 10: 'v_frameshift'
    # 11: 'Vintact'/'Vtruncated'
    # 12: 'junction'
    # 13: 'junction_aa'
    # 14: 'CDR3aa'
    # 15: 'readframe'
    # 16: 'pcov'
    # 17: 'pid'
    # 18: 'transl'

    return result

#-------------------------------------------------------------------------------
#### main section
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('igblastOut_name', \
        help='Filename for the IgBLAST output (e.g., "source.igblast_out.airr.tsv")')
    #parser.add_argument('--debug', help='output debug information', action='store_true')
    args = parser.parse_args()
    keys = []

try:
    with open(args.igblastOut_name, encoding="utf8") as igblast:
        in_line = igblast.readline()

        while in_line:
            if re.search(r'^sequence_id', in_line):
                keys = in_line.split("\t")
            else:
                data = in_line.split("\t")
                
                # igblast_data = parse_igblast_block(igblast, in_line)
                # in_line = fasta.readline()
                annotated_fasta = compose_fasta_block(keys, data)
                print('>' + annotated_fasta['query_id'])
                print(annotated_fasta['query_seq'])
            in_line = igblast.readline()

except FileNotFoundError:
    if not exists(args.igblastOut_name):
        sys.exit('File ' + args.igblastOut_name + ' was not found!')