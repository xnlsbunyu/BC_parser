# Follow biopython FastqGeneralIterator generator function
def barcode_filter_generator(handle, read_direction):
    import re
    import numpy as np
    lintag1_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAA|A.AA|AT.A|ATA.)\D*')
    lintag2_re = re.compile('\D*?(.ACC|T.CC|TA.C|TAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAC|T.AC|TT.C|TTA.)\D*')
    if read_direction == "f":
        bc_filter = lintag1_re
    elif read_direction == "r":
        bc_filter = lintag2_re
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
        bc_grep = bc_filter.match(seq_string[57:57+38])

        if bc_grep != None:
            quality_score = np.fromstring(quality_string[bc_grep.start() + 57: bc_grep.end() + 57], np.int8)-33
            if np.mean(quality_score) >= 30:
                count += 1
                yield (title_line, seq_string, quality_string)
            else:
                continue
        else:
            continue
    print (count)
    raise StopIteration

for t, s, q in barcode_filter_generator(open("small_f.fastq", "r"), "f"):
    print (s)
