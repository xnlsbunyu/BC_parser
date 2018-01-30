# Follow biopython FastqGeneralIterator generator function
def barcode_filter_generator(handle, read_direction, multitag):  
    from regex_tools import mismatch_re, get_full_re
    import re
    import numpy as np
    lintag1_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*')
    lintag2_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAC|T.AC|TT.C|TTA.)\D*')
    f_clipper = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)') 
    bc_len = 38
    multitag_pos = 8
    
    # forward multitag big regex
    f_multi = [i[0:6] for i in multitag]
    # reverse multitag big regex
    r_multi = [i[6:] for i in multitag]
    if read_direction == "f":
        bc_filter = lintag1_re
        r_clipper = re.compile('\D*?(AAT.|AA.A|A.TA|.ATA)')
        bc_pos = 57
        multitag_filter = re.compile(get_full_re(f_multi))
    elif read_direction == "r":
        bc_filter = lintag2_re
        r_clipper = re.compile('\D*?(CAT.|CA.T|C.TT|.ATT)')
        bc_pos = 43
        multitag_filter = re.compile(get_full_re(r_multi))
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
