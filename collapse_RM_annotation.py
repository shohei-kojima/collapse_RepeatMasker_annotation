#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Usage: python %prog -i in.fa.out -o out
Requirement: Python 3.7 or later
Licence: see LICENSE file
'''

import os,sys,gzip,datetime,collections,argparse,errno

version='2021/09/20'

description='''
    This script collapses TE annotations from RepeatMasker (i.e. genome.fa.out file).
    The TE annotations from RepeatMasker are frequently overlap each other. If so, when counting
    NGS reads mapping to TEs, reads mapping to two or more TE annotations may not be counted.
    In such case, it is preferable to use non-overlapping TE annotations. This script makes
    non-overlapping TE annotations by removing TEs with lower bit score. Processing of one file
    (e.g. hg38, mm10) requires ~3 min and ~100 MB RAM. This script only uses single thread.
    Requirement: Python 3.7 or later
'''

# args
parser=argparse.ArgumentParser(description=description)
parser.add_argument('-i', metavar='str', type=str, help='Specify input genome.fa.out file.', required=True)
parser.add_argument('-o', metavar='str', type=str, help='Specify output basename. [basename].gtf.gz and [basename].bed.gz will be generated.', required=True)
parser.add_argument('-gap', metavar='int', type=int, help='Optional. Specify gap distance to connect two TEs belonging to the same repeat. Default: 0', default=0)
parser.add_argument('-keep_simple_repeat', help='Optional. Specify if you want to keep "Simple_repeat" and "Low_complexity".', action='store_true')
parser.add_argument('-quiet', help='Optional. Specify if you do not want to output processing status.', action='store_true')
parser.add_argument('-v', '--version', action='version', version='Version: %s %s' % (os.path.basename(__file__), version))
parser.add_argument('-testrun', action='store_true', help=argparse.SUPPRESS)
args=parser.parse_args()

# set up
ogtf='%s.gtf.gz' % args.o
obed='%s.bed.gz' % args.o
gap=args.gap
annots=('gene', 'transcript', 'exon')
remove_simple_low_complex={'Simple_repeat', 'Low_complexity'}
if args.keep_simple_repeat is True:
    remove_simple_low_complex={}
_date=datetime.datetime.now()

# check input and output
if os.path.exists(args.i) is False:
    print(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.i))
    exit(1)
outdir=os.path.dirname(args.o)
if '/' in outdir:
    os.makedirs(outdir, exist_ok=True)



class IntervalTree:
    def __init__(self):
        self.root=None
        
    def insert(self, start, end, score, info):
        if self.root is None:
            self.root=IntervalNode(start, end, score, info)
        else:
            self.root=self.root.insert(start, end, score, info)
    
    def find_top_score(self, start, end):
        return self.root.find_top_score(start, end)


class IntervalNode:
    def __init__(self, start, end, score, info):
        self.start=start
        self.end=end
        self.score=score
        self.info=info
        self.maxend=self.end
        self.minend=self.end
        self.left=None
        self.right=None

    def insert(self, start, end, score, info):
        root=self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right=self.right.insert(start, end, score, info)
            else:
                self.right=IntervalNode(start, end, score, info)
            # rebalance tree
            if self.score < self.right.score:
                root=self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left=self.left.insert(start, end, score, info)
            else:
                self.left=IntervalNode(start, end, score, info)
            # rebalance tree
            if self.score < self.left.score:
                root=self.rotateright()
        if root.right and root.left:
            root.maxend=max(root.end, root.right.maxend, root.left.maxend)
            root.minend=min(root.end, root.right.minend, root.left.minend)
        elif root.right:
            root.maxend=max(root.end, root.right.maxend)
            root.minend=min(root.end, root.right.minend)
        elif root.left:
            root.maxend=max(root.end, root.left.maxend)
            root.minend=min(root.end, root.left.minend)
        return root
    
    def rotateright(self):
        root=self.left
        self.left=self.left.right
        root.right=self
        if self.right and self.left:
            self.maxend=max(self.end, self.right.maxend, self.left.maxend)
            self.minend=min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend=max(self.end, self.right.maxend)
            self.minend=min(self.end, self.right.minend)
        elif self.left:
            self.maxend=max(self.end, self.left.maxend)
            self.minend=min(self.end, self.left.minend)
        return root

    def rotateleft(self):
        root=self.right
        self.right=self.right.left
        root.left=self
        if self.right and self.left:
            self.maxend=max(self.end, self.right.maxend, self.left.maxend)
            self.minend=min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend=max(self.end, self.right.maxend)
            self.minend=min(self.end, self.right.minend)
        elif self.left:
            self.maxend=max(self.end, self.left.maxend)
            self.minend=min(self.end, self.left.minend)
        return root
    
    def find_top_score(self, start, end):
        if start < self.end and end > self.start:
            return (self.score, self.info)
        if self.left and start < self.left.maxend:
            return self.left.find_top_score(start, end)
        if self.right and end > self.start:
            return self.right.find_top_score(start, end)


def connect(l):
    connected=[]
    prev_e= (gap * -1) - 1
    prev_i=None
    for _rep in l:
        dist= _rep.start - prev_e
        if _rep.info == prev_i and dist <= gap:
            prev_e=_rep.end
        else:
            if not prev_i is None:
                connected.append(Rep(prev_s, prev_e, prev_r, prev_i))
            prev_s,prev_e,prev_r,prev_i=_rep
    connected.append(Rep(prev_s, prev_e, prev_r, prev_i))
    return connected


def collapse(tmp, chr, gft_id_n):
    if len(tmp) == 1:
        connected=tmp
    else:
        tree=IntervalTree()
        poss=set()
        for _rep in tmp:
            tree.insert(*_rep)
            poss |= {*_rep[:2]}
        poss=sorted(list(poss))
        results=[]
        for n in range(len(poss) - 1):
            se=poss[n:n + 2]
            result=tree.find_top_score(*se)
            if not result is None:
                results.append(Rep(*se, *result))
        connected=connect(results)
    gtf_lines=[]
    bed_lines=[]
    for _rep in connected:  # start, end, score, info
        strand,rep_name=_rep.info
        gft_id='RM_%s.%s' % (gft_id_n, rep_name)
        gene_attr='gene_id "%s"; gene_name "%s"; score "%s";' % (gft_id, gft_id, _rep.score)
        tran_attr='gene_id "%s"; transcript_id "t_%s"; gene_name "%s"; score "%s";' % (gft_id, gft_id, gft_id, _rep.score)
        exon_attr='gene_id "%s"; transcript_id "t_%s"; gene_name "%s"; exon_number 1; exon_id "e_%s"; score "%s";' % (gft_id, gft_id, gft_id, gft_id, _rep.score)
        attrs=(gene_attr, tran_attr, exon_attr)
        for annot,attr in zip(annots, attrs):
            l=[chr, 'RepeatMasker', annot, str(_rep.start), str(_rep.end), '.', strand, '.', attr]
            gtf_lines.append('\t'.join(l) +'\n')
        l=[chr, str(_rep.start), str(_rep.end), gft_id, str(_rep.score), strand]
        bed_lines.append('\t'.join(l) +'\n')
        gft_id_n += 1
    return gtf_lines, bed_lines, gft_id_n


def parse_line(ls):
    start= int(ls[5]) - 1
    end= int(ls[6])
    score= int(ls[0])
    strand=ls[8]
    if strand == 'C':
        strand='-'
    repname='%s.%s.%s' % (ls[9], ls[10], ls[14])
    info=(strand, repname)
    return Rep(start, end, score, info)


def per_chr(reps, gft_id_n):
    if args.quiet is False:
        print('processing %s...' % prev_chr)
    prev_end= (gap * -1) - 1
    reps=sorted(reps)
    for _rep in reps:
        dist= _rep[0] - prev_end
        if dist <= gap:
            tmp.append(_rep)
        else:
            if prev_end > 0:
                gtf_lines,bed_lines,gft_id_n=collapse(tmp, prev_chr, gft_id_n)
                outgtf.write(''.join(gtf_lines))
                outbed.write(''.join(bed_lines))
            tmp=[_rep]
        prev_end=_rep[1]
    gtf_lines,bed_lines,gft_id_n=collapse(tmp, prev_chr, gft_id_n)
    outgtf.write(''.join(gtf_lines))
    outbed.write(''.join(bed_lines))
    return gft_id_n



# main
if args.testrun is False:
    outgtf=gzip.open(ogtf, 'wt')
    outbed=gzip.open(obed, 'wt')
    
    outgtf.write('##format: gtf\n')
    outgtf.write('##date: %s\n' % _date)
    outgtf.write('##version: %s %s\n' % (os.path.basename(__file__), version))
    outgtf.write('##original file: %s\n' % args.i)
    outbed.write('#date: %s\n' % _date)
    outbed.write('#version: %s %s\n' % (os.path.basename(__file__), version))
    outbed.write('#original file: %s\n' % args.i)
    outbed.write('#format: chr start end name bit_score strand\n')
    
    prev_chr='dummy'
    gft_id_n=0
    Rep=collections.namedtuple('Rep', ['start', 'end', 'score', 'info'])
    with open(args.i) as infile:
        for _ in range(3):
            next(infile)
        for line in infile:
            ls=line.split()
            if ls[10] in remove_simple_low_complex:
                continue
            _rep=parse_line(ls)
            # per chr processing
            if not ls[4] == prev_chr:
                if not prev_chr == 'dummy':
                    gft_id_n=per_chr(reps, gft_id_n)
                prev_chr=ls[4]
                reps=[_rep]
            else:
                reps.append(_rep)
    gft_id_n=per_chr(reps, gft_id_n)
    
    outgtf.close()
    outbed.close()
    print('\n%s repeats were retained.\nThank you for using this script!\n' % gft_id_n)

else:
    # for debug
    Rep=collections.namedtuple('Rep', ['start', 'end', 'score', 'info'])
    tree=IntervalTree()
    tree.insert(100, 200, 55, ('+', 'LTR'))
    tree.insert(110, 120, 40, ('-', 'SINE'))
    tree.insert(190, 290, 60, ('-', 'LINE'))
    tree.insert(195, 300, 80, ('+', 'DNA'))
    poss=[100, 150, 190, 195, 200, 290, 300]
    results=[]
    for n in range(len(poss) - 1):
        se=poss[n:n + 2]
        result=tree.find_top_score(*se)
        results.append(Rep(*se, *result))
    final=connect(results)
    print(final)
