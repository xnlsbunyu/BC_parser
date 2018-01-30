#Define a function to generate all possible regex with allowed number of mismatches
def mismatch_re(s, mis):
    import itertools
    import re
    all_list_f = []
    all_list_r = []
    mismatch_pos = itertools.combinations(range(6), mis)
    for i in mismatch_pos:
        temp_s_f, temp_s_r = s[0:6], s[6:]
        for pos in i:
            temp_s_f = "".join((temp_s_f[:pos], ".", temp_s_f[pos+1:]))
            temp_s_r = "".join((temp_s_r[:pos], ".", temp_s_r[pos+1:]))
            all_list_f.append("^@" + temp_s_f)
            all_list_r.append(temp_s_r + "@$")
    return (re.compile('|'.join(all_list_f)),re.compile('|'.join(all_list_r)))
#    with open(s[0:6] + "_re.txt", "w") as f:
#        f.writelines(''.join(all_list_f))
#    f.close()
#    with open(s[6:] + "_re.txt", "w") as r:
#        r.writelines(''.join(all_list_r))
#    r.close()

multitag = ['GCTTGACTTCTC','CAAGAAAGGGTT', 'ATGCGAGGCGTC', 'CATAGGGTATTG', 'TTTTTGGCGCGC', 'TGGAAGGCACTC', 'CCGGGCTGTTTG', 'CTACTCCCTGCA', 'AGTGATGATCTG', 'TTCGGTGCTTAA']
multitag_f_dict = {}
multitag_r_dict = {}
for i in multitag:
    multitag_f_dict[i[0:6]], multitag_r_dict[i[6:]] = mismatch_re(i, 2)

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




