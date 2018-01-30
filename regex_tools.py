# small functions for generating regex patterns with mismatches
def mismatch_re(s, mis, prefix, suffix):
    import itertools
    import re
    all_list = []
    mismatch_pos = itertools.combinations(range(len(s)), mis)
    for i in mismatch_pos:
        temp_s = s
        for pos in i:
            temp_s = "".join((temp_s[:pos], ".", temp_s[pos+1:]))
        all_list.append(prefix+ temp_s + suffix)
    return '|'.join(all_list)

# get full list of multitag used in one run and join them together to generate one big regex to filter multitag
def get_full_re(multitag):
    import re
    multitag_list = []
    for i in multitag:
        multitag_list.append(mismatch_re(i, 2, "", ""))
    return '|'.join(multitag_list)
