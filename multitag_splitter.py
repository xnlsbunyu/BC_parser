import re
from regex_tools import mismatch_re
multitag = ['GCTTGACTTCTC','CAAGAAAGGGTT', 'ATGCGAGGCGTC', 'CATAGGGTATTG', 'TTTTTGGCGCGC', 'TGGAAGGCACTC', 'CCGGGCTGTTTG', 'CTACTCCCTGCA', 'AGTGATGATCTG', 'TTCGGTGCTTAA']
multitag_f_dict = {}
multitag_r_dict = {}
for i in multitag:
    multitag_f_dict[i[0:6]], multitag_r_dict[i[6:]] = re.compile(mismatch_re(i[0:6], 2, "^@", "")), re.compile(mismatch_re(i[6:], 2, "", "@$"))

for i in multitag:
    vars()[i + "_lintag1"] = open(i + "_lintag1.txt", "w")
    vars()[i + "_lintag2"] = open(i + "_lintag2.txt", "w")
    vars()[i + "_seqtag"] = open(i + "_seqtag.txt", "w")
    vars()[i + "_lintag1_list"] = []
    vars()[i + "_lintag2_list"] = []
    vars()[i + "_seqtag_list"] = []

with open("everything") as f:
    while True:
        line = f.readline()
        if not line: break
        for (k1, v1), (k2, v2) in zip(multitag_f_dict.items(), multitag_r_dict.items()):
            if v1.match(line) is not None and v2.search(line) is not None:
                line_split = line.split(',')
                if k1 + k2 in multitag:
                    vars()[k1+k2+"_lintag1_list"].append(line_split[1] + "\n")
                    vars()[k1+k2+"_lintag2_list"].append(line_split[2] + "\n")
                    vars()[k1+k2+"_seqtag_list" ].append(line_split[3] + "\n")
                else:
                    continue

for i in multitag:
    vars()[i + "_lintag1"].writelines(vars()[i + "_lintag1_list"])
    vars()[i + "_lintag1"].close()
    vars()[i + "_lintag2"].writelines(vars()[i + "_lintag2_list"])
    vars()[i + "_lintag2"].close()
    vars()[i + "_seqtag"].writelines(vars()[i + "_seqtag_list"])
    vars()[i + "_seqtag"].close()




