#!/usr/bin/env python

import sys

host_framing_file=sys.argv[1]
virus_framing_file=sys.argv[2]

libname=host_framing_file.split(".framing.txt")[0]

host_frames_dict={}
host_frames_dict[1]=0
host_frames_dict[2]=0
host_frames_dict[3]=0

virus_frames_dict={}
virus_frames_dict[1]=0
virus_frames_dict[2]=0
virus_frames_dict[3]=0

total_host_reads=0
total_virus_reads=0

with open(host_framing_file) as framing_content:
    for line in framing_content:
        line=line.split()
        length=int(line[0].strip())
        no_reads=int(line[1].strip())
        frame=int(line[2].strip())
        total_host_reads+=no_reads
        # changing frame names from 0,1,2 to 1,2,3
        if frame==0:
            frame=1
        elif frame==1:
            frame=2
        elif frame==2:
            frame=3
        host_frames_dict[frame]+=no_reads

with open(virus_framing_file) as framing_content:
    for line in framing_content:
        line=line.split()
        length=int(line[0].strip())
        no_reads=int(line[1].strip())
        frame=int(line[2].strip())
        total_virus_reads+=no_reads
        # changing frame names from 0,1,2 to 1,2,3
        if frame==0:
            frame=1
        elif frame==1:
            frame=2
        elif frame==2:
            frame=3
        virus_frames_dict[frame]+=no_reads

print libname+"\t",

frames_list=[1,2,3]

	# virus_1	virus_2	virus_3	host_1	host_2	host_3

for frame in frames_list:
    no_virus_reads=virus_frames_dict[frame]
    virus_percent_in_frame=(float(no_virus_reads)/float(total_virus_reads))*100
    print str(virus_percent_in_frame)+"\t",

for frame in frames_list:
    no_host_reads=host_frames_dict[frame]
    host_percent_in_frame=(float(no_host_reads)/float(total_host_reads))*100
    print str(host_percent_in_frame)+"\t",


print "\n",
