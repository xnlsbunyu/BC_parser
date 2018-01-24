# Follow biopython FastqGeneralIterator generator function
def barcode_filter_generator(handle, read_direction):
    import re
    import numpy as np
    lintag1_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*')
    lintag2_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAC|T.AC|TT.C|TTA.)\D*')
    bc_len = 38
    if read_direction == "f":
        bc_filter = lintag1_re
        bc_pos = 57
    elif read_direction == "r":
        bc_filter = lintag2_re
        bc_pos = 43
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

        if bc_grep != None:
            quality_score = np.fromstring(quality_string[bc_grep.start() + bc_pos: bc_grep.end() + bc_pos], np.int8)-33
            if np.mean(quality_score) >= 30:
#                count += 1
                yield (title_line, seq_string, quality_string, count)
            else:
                continue
        else:
            continue
    print (count)
    raise StopIteration
bc1_dict = {}
bc2_dict = {}
for t, s, q, c in barcode_filter_generator(open("small_f.fastq", "r"), "f"):
    bc1_dict[c] = s
for t, s, q, c in barcode_filter_generator(open("small_r.fastq", "r"), "r"):
    bc2_dict[c] = s
dbc_keys = sorted(list(bc1_dict.keys() & bc2_dict.keys()))
vars()["f_bc.txt"] = open("f_bc.txt", "w")
vars()["r_bc.txt"] = open("r_bc.txt", "w")
for i in dbc_keys:
    vars()["f_bc.txt"].write(bc1_dict[i] + '\n')
    vars()["r_bc.txt"].write(bc2_dict[i] + '\n')
vars()["f_bc.txt"].close()
vars()["r_bc.txt"].close()
#bc1 = barcode_filter_generator(open("small_f.fastq", "r"), "f")
#bc2 = barcode_filter_generator(open("small_r.fastq", "r"), "r")
#for (t1, s1, q1), (t2, s2, q2) in zip(bc1, bc2):

