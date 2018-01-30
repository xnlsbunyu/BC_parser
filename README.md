# BC_parser
1. modified the fastq parser from biopython (FastqGeneralIterator), which is a generator function
2. Add quality score, and barcode regex, and multitag regex filter
3. Output a big file containing each record per line with defined format
4. Demultiplex the big file using multitag_splitter.py script.


