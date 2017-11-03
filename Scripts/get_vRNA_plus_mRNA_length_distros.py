#!/usr/bin/env python

import sys

mRNA_lengths_file=sys.argv[1]
vRNA_lengths_file=sys.argv[2]

mRNA_seqlen_count_dict={}
vRNA_seqlen_count_dict={}

for i in range(25,51):
    mRNA_seqlen_count_dict[i]=0
    vRNA_seqlen_count_dict[i]=0

mRNA_total_count=0
vRNA_total_count=0

with open(mRNA_lengths_file) as mRNA_content:
    for line in mRNA_content:
        line=line.split()
        no_reads=int(line[0].strip())
        seqlen=int(line[1].strip())
        if seqlen in mRNA_seqlen_count_dict.keys():
            mRNA_seqlen_count_dict[seqlen]+=no_reads
            mRNA_total_count+=no_reads

with open(vRNA_lengths_file) as vRNA_content:
    for line in vRNA_content:
        line=line.split()
        no_reads=int(line[0].strip())
        seqlen=int(line[1].strip())
        if seqlen in vRNA_seqlen_count_dict.keys():
            vRNA_seqlen_count_dict[seqlen]+=no_reads
            vRNA_total_count+=no_reads

print """Length_nt\tmRNA\tvRNA"""
for i in range(25,51):
    mRNA_count=mRNA_seqlen_count_dict[i]
    if mRNA_count>0:
        mRNA_percentage=(float(mRNA_count)/float(mRNA_total_count))*100
    else:
        mRNA_percentage=0
    vRNA_count=vRNA_seqlen_count_dict[i]
    if vRNA_count>0:
        vRNA_percentage=(float(vRNA_count)/float(vRNA_total_count))*100
    else:
        vRNA_percentage=0
    print """{0}\t{1}\t{2}""".format(i,mRNA_percentage,vRNA_percentage)
