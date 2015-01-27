#!/usr/bin/env python3

import iva
import argparse

parser = argparse.ArgumentParser(
    description = 'Wrapper to pick closest ref from reads',
    usage = '%(prog)s <db> <reads> <outprefix>')
parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1)
parser.add_argument('db', help='database dir, made by iva_qc_make_db')
parser.add_argument('reads', help='Reads fastq filename')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

db = iva.kraken.Database(options.db, preload=True, threads=options.threads)
ref = db.choose_reference(options.reads, options.outprefix)
with open(options.outprefix + '.ref.info', 'w') as f:
    print(ref, file=f)
f.close()
