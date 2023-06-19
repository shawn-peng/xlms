#!/bin/env python
import sys
import Bio.SeqIO
from Bio.Seq import Seq
from copy import copy

def gen_decoy_db(input_db):
    def get_output_list():
        yield from Bio.SeqIO.parse(input_db, 'fasta')

        for rec in Bio.SeqIO.parse(input_db, 'fasta'):
            rev_rec = copy(rec)
            rev_rec.id = 'reverse_' + rec.id
            rev_rec.seq = Seq(''.join(reversed(rec.seq)))
            yield rev_rec
    
    pos = input_db.rfind('.')
    output_db = input_db[:pos] + '_decoy' + input_db[pos:]
    print(output_db)
    Bio.SeqIO.write(get_output_list(), output_db, 'fasta')


if __name__ == '__main__':
    assert(len(sys.argv) >= 2)

    input_db = sys.argv[1]
    gen_decoy_db(input_db)


