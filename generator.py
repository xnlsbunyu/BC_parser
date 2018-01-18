# Author: Xianan Liu
import re
import os
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import sys
forward_file = sys.argv[1]
reverse_file = sys.argv[2]
f_file = FastqGeneralIterator(open(forward_file, "r"))
r_file = FastqGeneralIterator(open(reverse_file, "r"))

#Define a function to generate all possible regex with allowed number of mismatches
def mismatch_re(s, mis):
    import itertools
    import re
    all_list = []
    mismatch_pos = itertools.combinations(range(len(s)), mis)
    for i in mismatch_pos:
        temp_s = s
        for pos in i:
            temp_s = "".join((temp_s[:pos], ".", temp_s[pos+1:]))
        all_list.append(temp_s)
    return re.compile('|'.join(all_list))

multitag = ['GCTTGACTTCTC','CAAGAAAGGGTT', 'ATGCGAGGCGTCC', 'CATAGGGTATTG', 'TTTTTGGCGCGC', 'TGGAAGGCACTC', 'CCGGGCTGTTTG', 'CTACTCCCTGCA', 'AGTGATGATCTG', 'TTCGGTGCTTAA']

# generate the regex with mismatch for each multitag
re_dict = {}
for k, v in enumerate(multitag):
    re_dict[k] = [mismatch_re(v[0:6], 2), mismatch_re(v[6:], 2)]
multitag_dict = {}
for k, m in enumerate(multitag):
    multitag_dict[k] = m

# open file to write multitags
for i in multitag:
    vars()[i + '_multitag'] = open(i + '_multitag.txt', 'w')
    vars()[i + '_lintag1'] = open(i + '_lintag1.txt', 'w')
    vars()[i + '_lintag2'] = open(i + '_lintag2.txt', 'w')
    vars()[i + '_seqtag'] = open(i + '_seqtag.txt', "w")
    vars()[i + '_lintag1_umi'] = open(i + '_lintag1_umi.txt', "w")
    vars()[i + '_lintag2_umi'] = open(i + '_lintag2_umi.txt', "w")
#f_multitag = mismatch_re(multitag[0:6], 1)
#r_multitag = mismatch_re(multitag[6:], 1)
lintag1_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*')
lintag2_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAC|T.AC|TT.C|TTA.)\D*')
lintag1_f_clipper = re.compile('(.ACC|T.CC|TA.C|TAC.)')
lintag2_f_clipper = re.compile('(.ACC|T.CC|TA.C|TAC.)')
lintag1_r_clipper = re.compile('(.TAC|T.AC|TT.C|TTA.)')
lintag2_r_clipper = re.compile('(.TAC|T.AC|TT.C|TTA.)')

#define some boundaries
seqtag_pos = 0
f_multitag_pos = 8
f_barcode_pos = 59
r_multitag_pos = 8
r_barcode_pos = 45
barcode_length = 38

# Minimum quality score
min_qs = 30

#define several important counters
bar_pattern_counts = 0
non_bar_pattern_counts = 0
bar_pattern_quality_counts = 0
bar_pattern_non_quality_counts = 0
bar_pattern_quality_multitag_counts = 0
bar_pattern_quality_non_multitags_counts = 0
total_counts = 0


# Use FastqGeneralIterator to parse each fastq as tuples Each tuple contains the title, the read sequence and the quality scores as strings
for f_record, r_record in zip(f_file, r_file):
    # read the sequence by indexing the tuples
    fr = f_record[1]
    rr = r_record[1]
    # read the quality scores and convert it into a list of integers
    fq = [ord(i) - 33 for i in list(f_record[2])]
    rq = [ord(i) - 33 for i in list(r_record[2])]
    lintag1_grep = lintag1_re.match(str(fr)[f_barcode_pos:f_barcode_pos + barcode_length])
    lintag2_grep = lintag2_re.match(str(rr)[r_barcode_pos:r_barcode_pos + barcode_length])
    total_counts += 1
    #After trial, checking barcode pattern is preferred
    if lintag1_grep is not None and lintag2_grep is not None:
        bar_pattern_counts += 1
        if np.mean(fq[lintag1_grep.start() : lintag1_grep.end()]) >= min_qs and np.mean(rq[lintag2_grep.start() : lintag2_grep.end()]) > min_qs:
            bar_pattern_quality_counts += 1
            for k, v in re_dict.items():
                f_grep = v[0].match(fr[f_multitag_pos:])
                r_grep = v[1].match(rr[r_multitag_pos:])
                if f_grep is not None and r_grep is not None:
                    bar_pattern_quality_multitag_counts += 1

                    vars()[multitag_dict[k] + '_lintag1'].write(lintag1_grep.group() + '\n')
                    vars()[multitag_dict[k] + '_lintag2'].write(lintag2_grep.group() + '\n')
                    vars()[multitag_dict[k] + '_multitag'].write(f_grep.group() +  r_grep.group() + '\n')
                    vars()[multitag_dict[k] + '_lintag1_umi'].write(lintag1_grep.group() + "," + fr[0:f_multitag_pos] + rr[0:r_multitag_pos] + '\n')
                    vars()[multitag_dict[k] + '_lintag2_umi'].write(lintag1_grep.group() + "," + fr[0:f_multitag_pos] + rr[0:r_multitag_pos] + '\n')
                    vars()[multitag_dict[k] + '_seqtag'].write(fr[0:f_multitag_pos] + rr[0:r_multitag_pos] + '\n')
                    break
        else:
            bar_pattern_non_quality_counts += 1
    else:
        non_bar_pattern_counts += 1

print ('The total counts of reads is ' + str(total_counts))
print ('The barcode pattern matching reads is ' + str(bar_pattern_counts))
print ('The non barcode pattern matching reads is ' + str(non_bar_pattern_counts))
print ('The barcode pattern matching and quality reads is ' + str(bar_pattern_quality_counts))
print ('The barcode pattern matching and non quality reads is ' + str(bar_pattern_non_quality_counts))
print ('The barcode pattern matching, quality, and multitag matching reads is ' + str(bar_pattern_quality_multitag_counts))


#Close the files
for i in multitag:
    vars()[i + '_multitag'].close()
    vars()[i + '_lintag1'].close()
    vars()[i + '_lintag2'].close()

f_file.close()
r_file.close()

# run bartender to cluster the barcodes
#for i in multitag:
#    os.system("bartender_single_com -f " + i + "_lintag1.txt -c 2 -t 8 -d 2 -o " + i + "_bartender_lintag1_stdout" )
#    os.system("bartender_single_com -f " + i + "_lintag2.txt -c 2 -t 8 -d 2 -o " + i + "_bartender_lintag2_stdout")




