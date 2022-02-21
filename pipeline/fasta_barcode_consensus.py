#!/usr/bin/python3
'''
 fasta_barcode_consensus.py
   Determine consensus sequences for (previously organized) reads with the same
   barcodes. Implemented using a position frequency matrix; quality scores derived as
   Cumulative Quality Score. This avoids introduction of gaps and significantly
   improves runtime (real alignment takes a long time).
'''
# Q-score encoding conventions, taken from Wikipedia:
# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
# ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
# ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
# .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ.....................
# LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
# |                         |    |        |                              |                     |
# 33                        59   64       73                            104                   126
# 0........................26...31.......40
#                          -5....0........9.............................40
#                                0........9.............................40
#                                   3.....9..............................41
# 0.2......................26...31........41
#
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 41)
#    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
#    (Note: See discussion above).
# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

import sys
import argparse
import re

#-------------------------------------------------------------------------------
def get_seed_middle (seq, half_seed_len, offset):
    '''
    returns string containing seed sequence
    1st argument--sequence
    2nd argument--half-seed length (e.g., 10)
    3rd argument--offset with respect to middle
    '''
    # start in the center of the sequence
    start = int(len(seq) / 2) - half_seed_len + offset - 1
    end = start + half_seed_len * 2
    return seq[start:end]

#-------------------------------------------------------------------------------
def consensus_generator (input_seqs, half_seed_len, offset_rng, max_mismatch_cnt, debug_flag = 0):
    '''
    generate consensus sequence/quality score pair from an array of arrays of strings
    1st argument--arrays of sequences
    2nd argument--half-seed length (e.g., 10)
    3rd argument--offset range (e.g., 5)
    4th argument--maxMismatch (e.g., 3)
    5th argument--debug-flag (e.g., 0 or 1)
    returns array of consensus sequence, quality, number of sequences used
    adapted from the MIGEC code (PMID: 24793455)
    https://github.com/mikessh/migec/blob/master/src/main/groovy/com/milaboratory/migec/Assemble.groovy
    '''
    nt_codes = { 'N' : 0, 'A' : 1, 'T' : 2, 'C' : 3, 'G' : 4 }
    code_nts = ['N','A','T','C','G']
    seed_dict = {}
    valid_seq_refs = []
    seqs_with_offsets = []
    best_seed = ''
    best_seed_data = [0,0]
    max_left_arm = 0
    max_right_arm = 0
    removed_elements =[]

    # find all possible seed sequences and their offsets, check lengths and find maximum
    #  storing the indices for later reference from the input array
    for ind, seq in enumerate(input_seqs):
        # make sure the sequence is long enough
        if len(seq) > (half_seed_len * 2 + offset_rng + 1):
            for offset in range(-offset_rng, offset_rng + 1):
                seed = get_seed_middle(seq, half_seed_len, offset)
                # initialize seed_dict entry, if necessary
                if seed not in seed_dict:
                    seed_dict[seed] = [0,0]
                seed_dict[seed][0] += 1      # increment count
                seed_dict[seed][1] += offset # add to the cumulative offset
            valid_seq_refs.append(ind)
        else:
            removed_elements.append(len(input_seqs)-ind)

    # identify the most abundant seed sequence
    for seed_seq, seed_data in seed_dict.items():
        if seed_data[0] > best_seed_data[0] or \
                (seed_data[0] == best_seed_data[0] and\
                 seed_data[1] < best_seed_data[1]):
        ### TODO: and what do we do when (seed_dict[seed][0] == best_seed_data[0]
        ###          and seed_dict[seed][1] == best_seed_data[1])?
        ### Perl makes a random choice, Python3 doesn't.

        # true if seed is more abundant, or if it's a tie and the cumulative offset is smaller
            best_seed = seed_seq
            best_seed_data = seed_data.copy()

    # identify the best match to seed for each sequence, load into the position frequency matrix
    for ind in valid_seq_refs:
        best_offset = 0
        best_mismatch_cnt = half_seed_len * 2 + 1 # max out the mismatch_cnt

        for offset in range(-offset_rng, offset_rng + 1):
            mismatch_cnt = 0
            seed = get_seed_middle(input_seqs[ind], half_seed_len, offset)

            if seed == best_seed:
                best_offset = offset
                best_mismatch_cnt = 0
                break

            mismatch_cnt = sum(c1!=c2 for c1,c2 in zip(seed,best_seed))
            if mismatch_cnt < best_mismatch_cnt:
                best_mismatch_cnt = mismatch_cnt
                best_offset = offset

        if best_mismatch_cnt <= max_mismatch_cnt:
            left_arm = int(len(input_seqs[ind])/2) + best_offset
            right_arm = len(input_seqs[ind]) - left_arm
            max_left_arm = left_arm if max_left_arm < left_arm else max_left_arm
            max_right_arm = right_arm if max_right_arm < right_arm else max_right_arm
            # load the sequence set with information for padding and the reference to the sequence
            seqs_with_offsets.append([left_arm,right_arm,ind])
        elif debug_flag:
            seed = get_seed_middle(input_seqs[ind], half_seed_len, best_offset)
            print ("#### Tossing element", len(input_seqs)-ind, "with", best_mismatch_cnt,\
                   "mismatches for", seed, "with respect to", best_seed, ":")
            print (input_seqs[ind])
            for offset in range(-offset_rng, offset_rng + 1):
                seed = get_seed_middle(input_seqs[ind], half_seed_len, offset)
                mismatch_cnt = sum(c1!=c2 for c1,c2 in zip(seed,best_seed))
                print(offset, ":", seed, "--", mismatch_cnt)


    if debug_flag:
        print('#### Array size:', len(seqs_with_offsets), '.####')
        print('#### Original sequence set size:', len(input_seqs), '.####')
        if removed_elements is not None:
            print('#### Removed the following elements (too short)', removed_elements)

    if len(seqs_with_offsets) > 1:
        # position_freq_matrix length should be the sum of maximum arms
        #    from the sequences in the alignment dataset
        pfm_length = max_left_arm + max_right_arm
        # initialize the position_freq_matrix array
        position_freq_matrix = []
        pfm_row = []
        for nt_code in range(len(code_nts)):
            pfm_row.append(0)
        for pos in range(pfm_length):
            position_freq_matrix.append(pfm_row.copy())

        if debug_flag:
            print("#### Alignment for the current MIG. ####")

        # load the position frequency matrix
        for ind, seq_plus_offsets in enumerate(seqs_with_offsets):
            seq = 'N'*(max_left_arm - seq_plus_offsets[0]) +\
                input_seqs[seq_plus_offsets[2]] + 'N'*(max_right_arm -\
                                                          seq_plus_offsets[1])
            for pos in range(pfm_length):
                position_freq_matrix[pos][nt_codes[seq[pos]]] += 1

            if debug_flag:
                print(">element" + str(len(input_seqs)-seq_plus_offsets[2]),\
                      "leftOffset =", seq_plus_offsets[0], "rightOffset =",\
                          seq_plus_offsets[1])
                print(seq)

        if debug_flag:
            print("########################################")
            print ("#### Contents of the positional frequency matrix for the current MIG.")
            print ("######## N\tA\tT\tC\tG")
            for row in position_freq_matrix:
                print("# values ", end='')
                for element in row:
                    print(element, end='\t')
                print()
            print('.....\n\n')

        consensus_seq = '' # initialize consensus sequence
        consensus_qual = '' # initialize consensus quality
        for pos in range(pfm_length):
            best_base = 0
            max_count = 0
            for nt_code in range(1,len(nt_codes)):
                if position_freq_matrix[pos][0]:
                    position_freq_matrix[pos][nt_code] += position_freq_matrix[pos][0]/4
                if max_count <  position_freq_matrix[pos][nt_code]:
                    max_count = position_freq_matrix[pos][nt_code]
                    best_base = code_nts[nt_code]
            consensus_seq += best_base

            if best_base != 'N':
                # cumulative quality score (CQS) from MIGEC #####
                best_base_qual = int(((max_count/len(seqs_with_offsets) - 0.25)) / 0.75 * 40 + 33)
                consensus_qual += chr(35 if best_base_qual < 35 else best_base_qual)
            else:
                consensus_qual += chr(35)
    else:
        # return a null value if dealing with a singlet (after tossing the bad sequences)
        return None

    # Prepare the return object
    consensus = [len(seqs_with_offsets), consensus_seq, consensus_qual]
    return consensus

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sourceName', help='Filename for the data set (e.g., "source.fasta")')
    parser.add_argument('-H', '--half_seed_length', nargs='?', type=int, default=10,\
                        help='Number of bases in the half-seed (default is 10 for a 21-base seed)')
    parser.add_argument('-M', '--max_mismatch_count', nargs='?', type=int, default=3,\
                        help='Number of mismatches to tolerate (default is 3)')
    parser.add_argument('-O', '--offset_range', nargs='?', type=int, default=5,\
                        help='Number of offsets to check in both directions (default is 5)')
    parser.add_argument('--debug', help='output debug information', action='store_true')
    args = parser.parse_args()

try:
    with open(args.sourceName, encoding="utf8") as sourceFile:
        clusterSeqAlignment = []
        consensusArray = []
        half_seed_length = args.half_seed_length
        offset_range = args.offset_range
        max_mismatch_count = args.max_mismatch_count

        # sys.exit('half_seed_length:' + str(half_seed_length) + ' offset_range:' +\
        #     str(offset_range) + ' max_mismatch_count:' + str(max_mismatch_count))

        header = sourceFile.readline()
        while header:
            header = header.rstrip()

            # only work with identifiable barcodes; dump the rest
            if re.search('barcode=unknown', header) is None:
                m = re.search(r'size=(\d+)[:;]element=(\d+)', header)
                if m is not None:
                    clusterSize = int(m.group(1))
                    elementID = int(m.group(2))
                    ############## FASTA harvest block ##############
                    sequence = sourceFile.readline()
                    if sequence is None:
                        sys.exit('!!!Error at:\n' + header)
                    sequence = sequence.rstrip()
                    #################################################
                    # start of a new cluster; initialize
                    if clusterSize == elementID:
                        counter = clusterSize
                        clusterID = re.sub('^>MIG','',header)
                        clusterID = re.sub(r'[:;]element=\d+','',clusterID)

                    if clusterSize >= 2:
                        clusterSeqAlignment.append(sequence)
                        counter -= 1
                    else:
                        print('@MIG' + clusterID + ";retained=1")
                        print(sequence + '\n' + '+')
                        print('#' * len(sequence))

                    # at the end of cluster: for non-singlets,
                    #   store the padded sequences in the alignment
                    if counter == 0:
                        consensusArray = consensus_generator(clusterSeqAlignment,\
                                                             half_seed_length, offset_range,\
                                                                 max_mismatch_count, args.debug)
                        clusterSeqAlignment = [] # clear the collection when done

                        if consensusArray is not None:
                            print('@MIG' + clusterID + ';retained=' + str(consensusArray[0]))
                            print(consensusArray[1] + '\n+')
                            print(consensusArray[2])
            header = sourceFile.readline()

except FileNotFoundError:
    sys.exit('File ' + args.sourceName + ' was not found!')
