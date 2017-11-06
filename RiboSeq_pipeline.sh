#!/bin/bash

# Set parameters:

workingdir=$(pwd)

# Directories
databasedir="X" # directory containing bowtie databases
stardbdir="X" # directory containing STAR databases
scriptsdir="X" # directory containing analysis scripts (python scripts, etc.)
plotsdir="X" # directory containing plot scripts (R scripts)

# Database names
databases1="rRNA:rRNA/rRNA"
databases2="mRNA:mRNA/mRNA ncRNA2:ncRNA_other/ncRNA_other gDNA:genome/genome"

# colours used for R plots:
dbdatacol1="rRNA:616 vRNA:552 mRNA:494 ncRNA_other:78"
#blue, red, light green, brown
dbdatacol2="gDNA:142"
#yellow

# Adaptor sequence to be trimmed from reads + min length of trimmed reads
adaptor="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCA"
trimlen=25
export adaptor
export trimlen

#-----------------------------------------------------------------------

# Place all gzipped Fastq files in the current working directory.
#    (files should have the extension .fq.gz)

# gunzip files in parallel
parallel 'gunzip {}' ::: *.fq.gz

#Check number of cycles (most abundant read length and number of reads of this
#  length in the first 100 reads in each file):

for library in $(awk '{print $1}' libraries.txt)
do
   head -400 $library.fq | \
    awk '{if (NR%4==2) print length($1)}' | sort -n | uniq -c | sort -nr | \
    head -1 | awk '{printf "%s: %s x %s nt\n","'"$library"'",$1,$2}'
done

#-------------------------------------------------------------------------------

# trim reads in parallel
parallel 'fastx_clipper -Q33 -l "$trimlen" -a "$adaptor" -c -n \
    -v -i {} > {.}.trimmed.fq 2>> {.}.log.txt' ::: *.fq

#The read names line in trimmed fastq files (e.g.
#  "@M00964:54:000000000-A3D68:1:1101:16462:1518 1:N:0:1") should contain two
#  fields. So if it is >2, just save the first 2 fields, and if it is only 1,
#  then add the dummy field "1:N:0:1" to each read name line.

for library in $(awk '{print $1}' libraries.txt)
do
  test=$(head -1 $library.trimmed.fq | awk '{print NF}') #check no. fields in 1st line
  if [ $test != 2 ]; then
    awk '{if (NR%4==1) {print $1,"1:N:0:1"} else {print $0}}' \
      $library.trimmed.fq | \
    awk '{if (NR%4==1) {print $1,$2} else {print $0}}' > $library.temp1
    mv $library.temp1 $library.trimmed.fq
  fi
done

#-----------------------------------------------------------------------

#Check that all bowtie databases are present for mapping

for host in $(awk '{print $3}' libraries.txt | sort | uniq)
do
  for dbline in $(echo $databases1 $databases2)
  do
    dbbowtie=$(echo $dbline | awk -F: '{print $2}')
    dbdir=$(ls $databasedir/$host/$dbbowtie.*ebwt | wc -l | awk '{print $1}')
    if [ -z $dbdir ]; then
      echo "Can't find bowtie database $databasedir/$host/$dbbowtie"
    fi
  done
done

for virus in $(awk '{print $2}' libraries.txt | sort | uniq)
do
  test=$(echo $virus | awk '{if ($1=="NA") {print 0} else {print 1}}')
  if [ -n $test ]; then
    dbdir=$(ls $databasedir/$virus/$virus.*ebwt | wc -l | awk '{print $1}')
    if [ -z $dbdir ]; then
      echo "Can't find bowtie database $databasedir/$virus/$virus"
    fi
  fi
done


#-----------------------------------------------------------------------

# map to ribosomal RNA

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "rRNA" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonrRNA.fq \
       $databasedir/$hostname/rRNA/rRNA \
       -q ./$library.trimmed.fq \
        > $library.rRNA.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done

#-----------------------------------------------------------------------

# map remaining reads to vRNA

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "vRNA" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonvRNA.fq \
       $databasedir/$virus/$virus \
       -q $library.nonrRNA.fq \
        > $library.vRNA.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done


#-----------------------------------------------------------------------

# map remaining reads to mRNA
for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "mRNA" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonmRNA.fq \
       $databasedir/$hostname/mRNA/mRNA \
       -q $library.nonvRNA.fq \
        > $library.mRNA.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done

#-----------------------------------------------------------------

# Map to other ncRNAs

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    echo "ncRNA_other" >> $library.log.txt
    bowtie -p 8 -v 2 --best --un $library.nonncRNA.fq \
       $databasedir/$hostname/ncRNA_other/ncRNA_other \
       -q $library.nonmRNA.fq \
        > ./$library.ncRNA_other.bowtie 2>> $library.log.txt
    echo >> $library.log.txt
done

#-----------------------------------------------------------------------

# map to the host genome with STAR

# NOTE: for reads with multiple mappings, we randomly select a single
# alignment using the flags: --outMultimapperOrder Random --outSAMmultNmax 1

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    STAR --runMode alignReads \
    --runThreadN 8  --outFileNamePrefix $library. --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 8 --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 2  \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outMultimapperOrder Random --outSAMmultNmax 1 --genomeLoad LoadAndKeep \
    --limitBAMsortRAM 60000000000 \
    --readFilesIn $library.nonncRNA.fq \
    --genomeDir $stardbdir/$hostname
done

#------------------------------------------------------------------------------

#Mapping summary plots

rm -f readcounts.summary

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    description="$library-$condition-$libtype"
    total=$(grep "Input:" ./$library.log.txt | head -1 | awk '{print $2}')
    tooshort=$(grep "too-short" ./$library.log.txt | head -1 | awk '{print $2}')
    adapters=$(grep "adapter-only" ./$library.log.txt | head -1 | awk '{print $2}')
    nonclipped=$(grep "non-clipped" ./$library.log.txt | head -1 | awk '{print $2}')
    echo -n "$description $total $tooshort $adapters $nonclipped " >> readcounts.summary
    for dbline in $(echo $dbdatacol1) #everything except gDNA
    do
        db=$(echo $dbline | awk -F: '{print $1}')
        fwdcounts=$(awk '{print $3}' $library.$db.bowtie | grep "+" | wc -l | awk '{print $1}')
        revcounts=$(awk '{print $3}' $library.$db.bowtie | grep "-" | wc -l | awk '{print $1}')
        echo -n "$db $fwdcounts $revcounts " >> readcounts.summary
    done
    for dbline in $(echo $dbdatacol2) #gDNA
    do
        db=$(echo $dbline | awk -F: '{print $1}')
        totalcounts=$(samtools view -F 4 -c $library.Aligned.sortedByCoord.out.bam | awk '{print $1}')
        echo -n "$db $totalcounts 0" >> readcounts.summary #NOTE: inserting 0 for reverse-strand count, since it doesn't make sense to distiguish strands for gDNA
    done
    echo >> readcounts.summary
done

ndb1=$(echo $dbdatacol1 | wc | awk '{print $2}')
ndb2=$(echo $dbdatacol2 | wc | awk '{print $2}')
ndb=$((ndb1+ndb2))
cat $plotsdir/readcounts_1.R | sed 's/ggg/'virus'/' | sed 's/nnn/'$ndb'/' > readcounts.summary.R
fwdcol=7 #forward mapping reads are specified in the 4th column of readcounts file
revcol=8
count=0
for dbline in $(echo $dbdatacol1)
do
    db=$(echo $dbline | awk -F: '{print $1}')
    col=$(echo $dbline | awk -F: '{print $2}')
    cat $plotsdir/readcounts_2.R | sed 's/ddd/'$db'/' | \
        sed 's/aaa/'$fwdcol'/' | sed 's/bbb/'$revcol'/' | \
        sed 's/ccc/'$col'/' | sed 's/xxx/'$count'/g' >> readcounts.summary.R
    fwdcol=$((fwdcol+3))
    revcol=$((revcol+3))
    count=$((count+1))
done
for dbline in $(echo $dbdatacol2)
do
    db=$(echo $dbline | awk -F: '{print $1}')
    col=$(echo $dbline | awk -F: '{print $2}')
    cat $plotsdir/readcounts_2.R | sed 's/ddd/'$db'/' | \
        sed 's/aaa/'$fwdcol'/' | sed 's/bbb/'$revcol'/' | \
        sed 's/ccc/'$col'/' | sed 's/xxx/'$count'/g' >> readcounts.summary.R
    fwdcol=$((fwdcol+3))
    revcol=$((revcol+3))
    count=$((count+1))
done

# make R plots
echo "dev.off()" >> readcounts.summary.R
sed 's/^.*://' readcounts.summary > readcounts.summary_pl
R --no-save --slave < readcounts.summary.R #argument --slave makes R run as 'quietly' as possible.


#-----------------------------------------------------------------------

# sort + count unmapped reads

for line in $(awk '{printf "%s:%s\n", $1,$2}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    hostname=$(echo $line | awk -F: '{print $2}')
    awk '{if (NR%4==2) print $0}' $library.Unmapped.out.mate1 \
      | sort | uniq -c | sort -nr > $library.unmapped.counted
done


#-----------------------------------------------------------------------------

# length distribution and phasing of reads on host mRNAs

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    description="$library-$condition-$libtype"
    awk '{if ($3=="+") print $4,$5+1,length($6)}' $library.mRNA.bowtie | \
      egrep -v "<|>|join" > $library.mRNAhits.total
    cat $library.mRNAhits.total | sed 's/_/ /g' | sed 's/\.\./ /' |  \
     sed 's/NM /NM_/' | sed 's/XM /XM_/' | \
    awk '{if (0+$5>=15+$2&&0+$5<=$3-15) print $6}'| sort -n | uniq -c \
       > $library.lenhist.mRNA.total
     sed 's/ttt/'$description'/' $plotdir/lengthdist.R | \
      sed 's/lll/'$library'/' | sed 's/:/ /g' > temp1
     R --no-save --slave < temp1
     rm -f temp1
done

#Framing of reads on host mRNAs.

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    description="$library-$condition-$libtype"
    echo "Frame of 5' end of reads that map within host mRNA coding sequences" \
      > $library.framing_by_length.txt
    echo >> $library.framing_by_length.txt
    rm -f $library.framing.txt #remove this file at the start of each loop
    for len in $(seq 25 35)
    do
        rm -f temp1 #remove these files at the start of each loop
        echo "$len nt reads:" >> $library.framing_by_length.txt
        echo "RefSeq mRNA reads" >> $library.framing_by_length.txt
        cat $library.mRNAhits.total | awk '{if ($3=="'"$len"'") print $0}' | \
          sed 's/_/ /g' | sed 's/\.\./ /' |  sed 's/NM /NM_/' | \
           sed 's/XM /XM_/' |  awk '{if (0+$5>=0+$2&&0+$5<=$3-30) print ($5-$2)%3}' | \
           awk '{s[$1]+=1}END{printf "%s 0\n%s 1\n%s 2\n",0+s[0],0+s[1],0+s[2]}' >> temp1
        paste temp1 >> $library.framing_by_length.txt
        paste temp1 | awk '{print "'"$len"'",$0}' >> $library.framing.txt
        echo >> $library.framing_by_length.txt
        rm -f temp1
    done
    cat $plotdir/framing.R | sed 's/ttt/'$description'/' | \
    sed 's/lll/'$library'/' | sed 's/@/ /g' > $library.framing.R
    R --no-save --slave < $library.framing.R
done


#Plot histograms to show positions of reads relative to start and stop codons

#Select transcripts with > 60 nt 5' and 3' UTRs and CDS >= 150 codons.
#Offset so that the A of AUG is at 0 ($5-$2), and similarly relative to
#  the Z of UXZ stop codon ($5-$3).
#Exclude reads with 5' end mapping more than 60 nt 5' of the AUG, or 3' end
#Find + plot 6 most abundant lengths.

for line in $(awk '{printf "%s@%s@%s\n", $1,$2,$4}' libraries.txt | grep -v RNASeq)
do
    library=$(echo $line | awk -F"@" '{print $1}')
    cat $library.mRNAhits.total | sed 's/_/ /g' | sed 's/\.\./ /' | \
        sed 's/NM /NM_/' | sed 's/XM /XM_/' | \
        awk '{if ($2-1>=0+60&&$4-$3>=0+60&&$3-$2+1>=0+450) \
        print $5-$2,$5-$3,$6}' | awk '{if (0+$1>=0-60&&$2+$3-1<=0+60) \
        print $0}' > $library.relative_to_starts_stops_length.total
        #as per above, $5=5' read position (on the transcript; NOT the genomic coordinate!), $2=start of gene (on the transcript; NOT the genomic coordinate!), $3=end of gene (on the transcript; NOT the genomic coordinate!), $4=length of gene.
        #the printed output here is in the format: 5' position relative to start, 5' position relative to stop, length of read
    sizes=$(awk '{print $3}' \
        $library.relative_to_starts_stops_length.total | \
        sort -n | uniq -c | sort -nr | head -6 | awk '{print $2}' | sort -n | \
        awk '{printf "%s,",$1}' | sed 's/,$//') #"head -6" to find the 6 most abundant read lengths
    cat $plotdir/relative_to_starts_and_stops.R | \
        sed 's/ttt/'$line'/' | sed 's/lll/'$library'/' | sed 's/@/ /g' | \
        sed 's/sss/'$sizes'/' | sed 's/mmm/'total'/' \
        > $library.relative_to_starts_and_stops.total.R
    R --no-save --slave < $library.relative_to_starts_and_stops.total.R
    cat $plotdir/relative_to_starts_and_stops_zm.R | \
        sed 's/ttt/'$line'/' | sed 's/lll/'$library'/' | sed 's/@/ /g' | \
        sed 's/sss/'$sizes'/' | sed 's/mmm/'total'/' \
        > $library.relative_to_starts_and_stops_zm.total.R
    R --no-save --slave < $library.relative_to_starts_and_stops_zm.total.R
done


#-------------------------------------------------------------------------

# Assess contamination by RNPs; i.e. proteins binding to RNAs and forming
#  complexes which happen to sediment with ribosomes and produce refseq_mRNA
#  footprints. These will have a different read length size distribution from
#  bona fide ribosome footprints. We can't recognize them on an individual
#  basis, but we can quantify the level of contamination by comparing the read
#  length size distributions of footprints in CDSs (mostly real ribosome
#  footprints) to footprints in regions expected to have very few bona fide
#  footprints (e.g. 3'UTRs).

#Mostly relevant to virus work (nucleocapsids non-specifically binding RNA) and
#  can gauge RNP contamination by comparing infected vs mock.

#Select refseq_mRNAs with > 100 codons 3'UTR and >150 codons CDS. Calculate density of
#  ribosomes in first 100 codons of 3'UTR relative to last 100 codons of CDS,
#  avoiding 10 codons around the stop codon. Also compare footprint size
#  distributions. Do for fwd and rev.

offset=12 #approximate P-site location
minlen=25
maxlen=50
rm -f RNPcontamination.txt
for line in $(awk '{printf "%s@%s@%s\n", $1,$4,$5}' libraries.txt)
do
  library=$(echo $line | awk -F"@" '{print $1}')
  awk '{if ($3=="+") print $4,$5+1+"'"$offset"'",length($6)}' \
    $library.mRNA.bowtie | egrep -v "<|>|join" | sed 's/_/ /g' | \
    sed 's/\.\./ /' | sed 's/NM /NM_/' | sed 's/XM /XM_/' | \
    awk '{if ($4-$3>=0+300&&1+$3-$2>=0+450) print $0}' > temp1.fwd
  awk '{print $5-$3,$6}' temp1.fwd > temp1.fwd.total
  sort temp1.fwd | uniq | awk '{print $5-$3,$6}' > temp1.fwd.uniq
  awk '{if ($3=="-") print $4,$5+1+"'"$offset"'",length($6)}' \
    $library.mRNA.bowtie | egrep -v "<|>|join" | sed 's/_/ /g' | \
    sed 's/\.\./ /' | sed 's/NM /NM_/' | sed 's/XM /XM_/' | \
    awk '{if ($4-$3>=0+300&&1+$3-$2>=0+450) print $0}' > temp1.rev
  awk '{print $5-$3,$6}' temp1.rev > temp1.rev.total
  sort temp1.rev | uniq | awk '{print $5-$3,$6}' > temp1.rev.uniq
  for mode in uniq total
  do
    for dir in fwd rev
    do
      awk '{if (0+$1>=0-300&&0+$1<=0+300) print $0}' temp1.$dir.$mode | \
        awk '{if (0+$1<=0-30||0+$1>=0+30) print $0}' > temp2
      awk '{if (0+$1<0+0) print $0}' temp2 > temp2.$dir.cds #if the mapping pos is upstream of the start codon
      awk '{if (0+$1>0+0) print $0}' temp2 > temp2.$dir.utr
      rm -f temp2
      for region in cds utr
      do
        rm -f $library.RNP.$dir.$region
        rm -f $library.RNP.$mode.$dir.$region
        for len in $(seq $minlen $maxlen)
        do
          num=$(awk '{if ($2=="'"$len"'") s+=1}END{print 0+s}'\
            temp2.$dir.$region)
          echo $len $num  >> $library.RNP.$mode.$dir.$region
        done
      done
    done
    nCDSfwd=$(wc -l temp2.fwd.cds | awk '{print $1}') #number of line in this file=number of reads being considered
    nUTRfwd=$(wc -l temp2.fwd.utr | awk '{print $1}')
    nCDSrev=$(wc -l temp2.rev.cds | awk '{print $1}')
    nUTRrev=$(wc -l temp2.rev.utr | awk '{print $1}')
    echo $line $mode $nCDSfwd $nUTRfwd $nCDSrev $nUTRrev | sed 's/@/ /g' \
      >> RNPcontamination.txt
    cat $plotsdir/RNPcontamination.R | sed 's/ttt/'$line'/' | \
      sed 's/lll/'$library'/' | sed 's/mmm/'$mode'/' | sed 's/@/ /g' \
      > $library.RNP.$mode.R
    R --no-save --slave < $library.RNP.$mode.R
    rm -f temp2.{fwd,rev}.{cds,utr} temp1
  done
  rm -f temp1.{fwd,rev}.{uniq,total} temp1.{fwd,rev}
done

grep -i Ribo RNPcontamination.txt | \
  awk '{printf "%s %s %s %s %.3f\n",$1,$2,$3,$4,$6/$5}'
#-> This shows the level of RNP contamination (assuming a baseline of zero
# footprint density in UTRs; mostly relevant for virus work where you can
#  compare infected vs mock to get the real baseline in 3'UTRs).

#Note that you will get fewer footprints in the 3'UTRs simply because the UTRs
#  are not all there - i.e. RNASeq 3'UTR(+10,+100) density is less than
#  RNASeq CDS(-10,-100) density. So should scale estimated contamination based
#  on the RiboSeq 3'UTR occupancy by the ratio
#  RNASeq CDS(-10,-100)/3'UTR(+10,+100)
grep -i RNA RNPcontamination.txt | \
  awk '{printf "%s %s %s %s %.3f\n",$1,$2,$3,$4,$6/$5}'
#Need to scale RiboSeq by the reciprocal of this
grep -i RNA RNPcontamination.txt | \
  awk '{printf "%s %s %s %s %.3f\n",$1,$2,$3,$4,$5/$6}' | grep total


#-----------------------------------------------------------------------

# For virus-infected samples, calculate length distribution and framing of
# reads on virus and host RNA for combined plotting

# Note; for refseq mRNAs, some CDSs have "<|>"s in the annotation, which
# could bias frame and start/stop proximal stats. So, here, we'll exclude all
# hits to CDSs with < or > in the CDS annotation ('egrep -v "<|>|join"')

# Print the mapping locations and lengths of vRNA & mRNA mapping reads
# select only forward strand matches to mRNAs ('awk '{if ($3=="+")...')

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt  | grep infected)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    awk '{if ($3=="+") print $4,$5+1,length($6)}' $library.vRNA.bowtie \
       > $library.vRNAhits.total
done

# Get lengths of reads mapping within coding sequences

# Exclude first and last 5 codons (15 nt) of each CDS, as these may contain
# initiating or terminating ribosomes.

# NOTE: the mRNA CDS coordinates have been automatically extracted from the
# RefSeq annotation. However, for the virus, we need to select a specific ORF
# (generally a highly expressed ORF is good) and explicitly specify the
# start and stop coordinates of that ORF, in addition to the
# chromosome accession from the fasta file used to build the virus
# bowtie index (e.g. "NC_001501") -- this must be the accession of the
# sequence containing the ORF of interest.

virus_chr_acc="X" # e.g. "NC_001501"
virus_ORF_start=X #621
virus_ORF_end=X #2237

for library in $(awk '{print $1}' libraries.txt | grep infected)
do
    cat $library.mRNAhits.total | sed 's/_/ /g' | sed 's/\.\./ /' |  \
     sed 's/NM /NM_/' | sed 's/XM /XM_/' | \
    awk '{if (12+$5>=15+$2&&12+$5<=$3-15) print $6}'| sort -n | uniq -c \
      > $library.lenhist.mRNA.CDS
    awk '{if ( ($1=="'$virus_chr_acc'") && (12+$2>=15+'$virus_ORF_start') && (12+$2<='$virus_ORF_end'-15) ) print $3}' \
     $library.vRNAhits.total | sort -n | uniq -c > $library.lenhist.vRNA.CDS
done

# calculate the % of reads at each length

scriptsdir="/home/adinan/Annotated_scripts/"
for library in $(awk '{print $1}' libraries.txt | grep infected)
do
  python $scriptsdir/get_vRNA_plus_mRNA_length_distros.py \
    $library.lenhist.mRNA.CDS $library.lenhist.vRNA.CDS \
    > $library.len.percs.combined
done

# Make a combined virus + host length distribution plot

cat $plotsdir/LengthDistros_top.R > length_distros_combined.R
nsamples=$(awk '/infected/{++c}END{print c}' libraries.txt) #count no. of infected samples
count=0

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt | grep infected)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    count=$((count+1))
    cat $plotsdir/LengthDistros_middle.R | \
     sed 's/lll/'$library'/' | sed 's/ppp/p'$count'/' |
     sed 's/nnn/'$count'/' >> length_distros_combined.R
done

cat $plotsdir/LengthDistros_bottom.R | sed 's/nnn/'$nsamples'/' \
  >> length_distros_combined.R

R --no-save --slave < length_distros_combined.R

#-----------------------------------------------------------------------

# for infected samples, make combined phasing plots for mRNA + vRNA

# get virus phasing
# again, need to specify the virus chromosome and start/stop coords
virus_chr_acc="X" # e.g. "NC_001501"
virus_ORF_start=X #621
virus_ORF_end=X #2237

for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt | grep infected)
do
    library=$(echo $line | awk -F: '{print $1}')
    awk '{if ( ($1=="'$virus_chr_acc'") && (12+$2>=15+'$virus_ORF_start') && (12+$2<='$virus_ORF_end'-15) ) print $3,$2,($2-'$virus_ORF_start')%3}' \
     $library.vRNAhits.total > $library.vRNA.framing
done

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "library" "virus_1" "virus_2" "virus_3" "host_1" "host_2" "host_3" >  framing_summarised.txt
for line in $(awk '{printf "%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5}' libraries.txt | grep infected)
do
    library=$(echo $line | awk -F: '{print $1}')
    python $scriptsdir/summarise_framing.py $library.framing.txt \
     $library.vRNA.framing >> framing_summarised.txt
done

cat $plotsdir/Framing_combined.R > Framing_combined.R
R --no-save --slave < Framing_combined.R
