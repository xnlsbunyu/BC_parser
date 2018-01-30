from regex_tools import mismatch_re, get_full_re
from barcode_filter_generator import barcode_filter_generator


# Build two separate dictionaries with line number as key and quality matching read as value
# What if the fastq file is very big, not sure the size limit for dictionary
# I need to optimize the dictionary, empty the dictionaries maybe every 100000 reads
def barcode_recorder(forward, reverse, multitag):
    import gzip
    bc1_dict = {}
    bc2_dict = {}
   # multitag = ['GCTTGACTTCTC','CAAGAAAGGGTT', 'ATGCGAGGCGTC', 'CATAGGGTATTG', 'TTTTTGGCGCGC', 'TGGAAGGCACTC', 'CCGGGCTGTTTG', 'CTACTCCCTGCA', 'AGTGATGATCTG', 'TTCGGTGCTTAA']
    multitag_f_dict = {}
    multitag_r_dict = {}
    for i in multitag:
        multitag_f_dict[i[0:6]] = get_full_re(i[0:6])
        multitag_r_dict[i[6:]] = get_full_re(i[6:])
    iteration = 1

    #open the file for writing everything 
    vars()["everything_test"] = open("everything_test.txt", "w")
    for (m1, bc1, umi1, c1), (m2, bc2, umi2, c2) in zip(barcode_filter_generator(gzip.open("small_f.fastq.gz", "rt"), "f", multitag), 
            barcode_filter_generator(gzip.open("small_r.fastq.gz", "rt"), "r", multitag)):
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

