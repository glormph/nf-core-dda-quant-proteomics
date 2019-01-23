#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq


def insilico_trypsinized(seq) :
    segments = []
    seg = []
    for i in range(len(seq)) :
        if seq[i] in ('K','R') :
            if i == len(seq)-1 :
                seg.append(seq[i])
            elif seq[i+1] == 'P' :
                seg.append(seq[i])
            else :
            #found first tryptic site
                if len(seg) :
                    segments.append(seg)
                segments.append( [seq[i]] )
                seg = []
        else :
            seg.append(seq[i])
    if len(seg) :
        segments.append(seg)
    segs_len = sum([len(x) for x in segments])
    try :
        assert(segs_len == len(seq))
    except Exception as e :
        segged_seq = []
        for s in segments :
            segged_seq.extend(s)
        print >> sys.stderr , "lens:" , len(seq), len(segged_seq)
        print >> sys.stderr , "original_seq:"
        print >> sys.stderr , "".join(seq)
        print >> sys.stderr , "new_seq:"
        print >> sys.stderr , "".join(segged_seq)
        raise(e)
    return segments


def tryp_rev(seq):
    segments = insilico_trypsinized(seq)
    final_seq = []
    for s in segments :
        if len(s) > 1 :
            if s[-1] in ['R', 'K']:
                new_s = s[:-1]
                new_s.reverse()
                new_s.append(s[-1])
            else:
                new_s = s
                new_s.reverse()
        else :
            new_s = s
        final_seq.extend(new_s)
    seq.seq = Seq(''.join(final_seq))
    seq.id = 'decoy_{}'.format(seq.name)
    return seq


def main():
    fa = sys.argv[1]
    out = os.path.join(os.path.split(fa)[0], 'decoy_{}'.format(os.path.basename(fa)))
    print(out)
    with open(fa) as fp, open(out, 'w') as wfp:
        seqs = SeqIO.parse(fp, 'fasta')
        SeqIO.write((tryp_rev(x) for x in seqs), wfp, 'fasta')

if __name__ == '__main__':
    main()
