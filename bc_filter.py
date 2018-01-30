# generate all possible combinations for a string with given number of mismatchs
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
    return '|'.join(all_list)

# get full list of multitag used in one run and join them together to generate one big regex to filter multitag
def get_full_re(multitag):
    import re
    multitag_list = []
    for i in multitag:
        multitag_list.append(mismatch_re(i, 2))
    return re.compile('|'.join(multitag_list))



# Follow biopython FastqGeneralIterator generator function
def barcode_filter_generator(handle, read_direction):
    import re
    import numpy as np
    lintag1_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*')
    lintag2_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAC|T.AC|TT.C|TTA.)\D*')
    f_clipper = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)') 
    bc_len = 38
    multitag_pos = 8
    multitag = ['GCTTGACTTCTC','CAAGAAAGGGTT', 'ATGCGAGGCGTC', 'CATAGGGTATTG', 'TTTTTGGCGCGC', 'TGGAAGGCACTC', 'CCGGGCTGTTTG', 'CTACTCCCTGCA', 'AGTGATGATCTG', 'TTCGGTGCTTAA']
    
    # forward multitag big regex
    f_multi = [i[0:6] for i in multitag]
    # reverse multitag big regex
    r_multi = [i[6:] for i in multitag]
    if read_direction == "f":
        bc_filter = lintag1_re
        r_clipper = re.compile('\D*?(AAT.|AA.A|A.TA|.ATA)')
        bc_pos = 57
        multitag_filter = get_full_re(f_multi)
    elif read_direction == "r":
        bc_filter = lintag2_re
        r_clipper = re.compile('\D*?(CAT.|CA.T|C.TT|.ATT)')
        bc_pos = 43
        multitag_filter = get_full_re(r_multi)
    count = 0
    handle_readline = handle.readline
    while True:
        line = handle_readline()
        if not line:
            return
        if line[0] == "@":
            break
        if isinstance(line[0], int):
            raise ValueError ("Is this handle in binary mode not text mode")

    while line:
        count += 1
        if line[0] != "@":
            raise ValueError(
                    "Records in Fastq files should start with '@' character"
                    )
        title_line = line[1:].rstrip()
        seq_string = handle_readline().rstrip()
        while True:
            line = handle_readline()
            if not line:
                raise ValueError("End of file without quality information.")
            if line[0] == "+":
                second_title = line[1:].rstrip()
                if second_title and second_title != title_line:
                    raise ValueError("Sequence and quality captions differ.")
                break
            seq_string += line.rstrip()
        if " " in seq_string or "\t" in seq_string:
            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)

        quality_string = handle_readline().rstrip()
        while True:
            line = handle_readline()
            if not line:
                break # EOF
            if line[0] == "@":
                if len(quality_string) >= seq_len:
                    break
            quality_string += line.rstrip()
        if seq_len != len(quality_string):
            raise ValueError("Lengths of sequence and quality values differs "
                                " for %s (%i and %i)."
                                % (title_line, seq_len, len(quality_string)))
        bc_grep = bc_filter.match(seq_string[bc_pos:bc_pos+bc_len])
        multitag_grep = multitag_filter.match(seq_string[multitag_pos:])
        if bc_grep != None and multitag_grep != None:
            quality_score = np.fromstring(quality_string[bc_grep.start() + bc_pos: bc_grep.end() + bc_pos], np.int8)-33
            if np.mean(quality_score) >= 29:
                raw_bc = bc_grep.group()
                bc_start = f_clipper.match(raw_bc).end()
                bc_end = r_clipper.search(raw_bc[::-1]).end()*(-1)
                yield (multitag_grep.group(), bc_grep.group()[bc_start:bc_end], seq_string[0:8], count)
            else:
                continue
        else:
            continue
    print (count)
    raise StopIteration



# Build two separate dictionaries with line number as key and quality matching read as value
# What if the fastq file is very big, not sure the size limit for dictionary
# I need to optimize the dictionary, empty the dictionaries maybe every 100000 reads
bc1_dict = {}
bc2_dict = {}
multitag = ['GCTTGACTTCTC','CAAGAAAGGGTT', 'ATGCGAGGCGTC', 'CATAGGGTATTG', 'TTTTTGGCGCGC', 'TGGAAGGCACTC', 'CCGGGCTGTTTG', 'CTACTCCCTGCA', 'AGTGATGATCTG', 'TTCGGTGCTTAA']
multitag_f_dict = {}
multitag_r_dict = {}
for i in multitag:
    multitag_f_dict[i[0:6]] = get_full_re(i[0:6])
    multitag_r_dict[i[6:]] = get_full_re(i[6:])
iteration = 1

#open the file for writing everything 
vars()["everything_test"] = open("everything_test.txt", "w")
for (m1, bc1, umi1, c1), (m2, bc2, umi2, c2) in zip(barcode_filter_generator(open("small_f.fastq", "r"), "f"), 
        barcode_filter_generator(open("small_r.fastq", "r"), "r")):
    bc1_dict[c1], bc2_dict[c2] = [m1, bc1, umi1], [m2, bc2, umi2]
    # set the buffer size when the dictionary size is too big, compare two dicts and write into everything file
    # use the line number as the filter, maybe use the length of the dictionary will be more reasonable
    if max(c1, c2) >10000*iteration: # every 10000 lines in the fastq parser
        #open everything file and write as append mode
        vars()["everything_test"] = open("everything_test.txt", "a")
        # get the common keys in these two dicts
        dbc_keys = sorted(list(bc1_dict.keys() & bc2_dict.keys()))
        # format for each line
        # Start with "@", following forward multitag, "," lintag1, ",", lintag2, "," , contatenated seqtag, reverse multitag ending with "@" and newline
        every_list = ["@" + bc1_dict[i][0] + "," + bc1_dict[i][1] + "," + bc2_dict[i][1] + "," + bc1_dict[i][2]+bc2_dict[i][2] + "," + bc2_dict[i][0] + "@\n" for i in dbc_keys]
        # write the whole list of lines at once into everything, avoid calling write function millions of times
        vars()["everything_test"].writelines(every_list)
        # This is avoiding missing records due to the fact that the keys of two dictionaries are not always paired
        if c1 > c2:
            bc1_dict = {k:v for k, v in bc1_dict.items() if k > c2}
            bc2_dict = {}
        elif c1 < c2:
            bc1_dict = {}
            bc2_dict = {k:v for k, v in bc2_dict.items() if k > c1}
        else:
            bc1_dict = {}
            bc2_dict = {}
        iteration += 1
    else:
        continue

# write the remaining contents of dictionaries into the file
dbc_keys = sorted(list(bc1_dict.keys() & bc2_dict.keys()))
#
every_list = ["@" + bc1_dict[i][0] + "," + bc1_dict[i][1] + "," + bc2_dict[i][1] + "," + bc1_dict[i][2]+bc2_dict[i][2] + "," + bc2_dict[i][0] + "@\n" for i in dbc_keys]
vars()["everything_test"].writelines(every_list)
vars()["everything_test"].close()

