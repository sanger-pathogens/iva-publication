#!/usr/bin/env python3

import argparse
import fastaq

parser = argparse.ArgumentParser(
    description = 'Gets N (default 200) reads equally spaced from input fastq file',
    usage = '%(prog)s <infile> <outfile>')
parser.add_argument('infile', help='Input fastq file')
parser.add_argument('outfile', help='Name of output file')
parser.add_argument('--reads', type=int, help='Number of reads wanted [%(default)s]', default=200, metavar='INT')
options = parser.parse_args()

wc = int(fastaq.utils.syscall_get_stdout('wc -l ' + options.infile)[0].split()[0])
total_reads = wc / 4
step = int(total_reads / (options.reads + 2))
f = fastaq.sequences.file_reader(options.infile)
n = 1
fout = fastaq.utils.open_file_write(options.outfile)
found = 0

for seq in f:
    if n % step == 0 and n != step:
        print(seq, file=fout)
        found += 1
    n += 1
    if found == options.reads:
        break

fastaq.utils.close(fout)
