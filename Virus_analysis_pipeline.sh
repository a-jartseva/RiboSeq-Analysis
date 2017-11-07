#!/bin/bash

# NOTE:
# run this script from within the directory where you performed the
# basic RiboSeq analysis pipeline (i.e. where your vRNA.bowtie files are located)

# Most of the R scripts are given in the Plots directory.
# However, you need to make "head" R files for the virus that contains
# information on the genome (i.e. locations of genes, etc)
# -- follow the structure of the "sample_head.R" and "Sample_head_linear.R" files
# in the Plots directory and use these as templates to make your virus head files


#-------------------------------------------------------------------------------

# Set parameters:

workingdir=$(pwd)

# Directories
databasedir="X" # directory containing bowtie databases
stardbdir="X" # directory containing STAR databases
scriptsdir="X" # directory containing analysis scripts (python scripts, etc.)
plotsdir="X" # directory containing plot scripts (R scripts)

#-------------------------------------------------------------------------------

# create a directory for virus analysis

mkdir virus_analyses && cd virus_analyses

# Make group files - i.e. group together RNAseq and RiboSeq libraries
for group in $(awk '{print $5}' ../libraries.txt | sort | uniq)
do
  awk '{if ($5=="'"$group"'") print $1,$4}' ../libraries.txt > $group.txt
done
awk '{print $5}' ../libraries.txt | sort | uniq > groups.txt

# Histogram files for virus genome coverage

mkdir HistFiles

for library in $(awk '{print $1}' ../libraries.txt)
do
  awk '{if ($3=="+") print $5+1}' ../$library.vRNA.bowtie | \
     sort -n | uniq -c > HistFiles/$library.fwd
  awk '{if ($3=="-") print $5+1}' ../$library.vRNA.bowtie | \
     sort -n | uniq -c > HistFiles/$library.rev
done

# Next, fill in nt positions with no mapped reads as '0'

mkdir HistFiles_wZeros
rm -f temp9
for library in $(awk '{print $1}' ../libraries.txt)
do
  for dir in fwd rev
  do
    cat HistFiles/$library.$dir >> temp9 # concatenate histogram files to a temp file
  done
done

# get min and max positions in the genome with read coverage from temp file
minnt=$(awk '{print $2}' temp9 | sort -n | head -1)
maxnt=$(awk '{print $2}' temp9 | sort -n | tail -1)
rm -f temp9
for library in $(awk '{print $1}' ../libraries.txt)
do
  for dir in fwd rev
  do
    rm -f HistFiles_wZeros/$library.$dir
    count=$minnt
    while [ $count -le $maxnt ]
    do
      awk '{if ($2=="'"$count"'") {print $1} else {print 0}}' \
         HistFiles/$library.$dir | uniq | sort -nr | head -1 | \
          awk '{print $1,"'"$count"'"}' >> HistFiles_wZeros/$library.$dir
      ((count+=1))
    done
  done
done

#------------------------------------------------------------------------------

# calculate genome coverage in sliding windows of various sizes (for plotting)

for hw in 1 7 15 150 #hw=half window size (full window is 2*hw + 1)
do
  sw=$(echo $hw | awk '{print 2*$1+1}') #this is the full slidingwindow (sw) size
  rm -rf SlidingWindow_$sw
  mkdir SlidingWindow_$sw
  for library in $(awk '{print $1}' ../libraries.txt)
  do
    for dir in fwd rev
    do
      awk '{print $1}' HistFiles_wZeros/$library.$dir > temp9 #put the read counts in here
      $scriptsdir/slidingwindow temp9 $hw > temp10
      paste HistFiles_wZeros/$library.$dir temp10 | awk '{print $4,$2}' \
         > SlidingWindow_$sw/$library.$dir
      rm -f temp9 temp10
    done
  done
done

#Separate HistFiles_wZeros files by frame (genomic coord) and remake sw.
#  sw and hw are now measured in 'codons'.
for hw in 1 7 15 150 #hw=half window size (full window is 2*hw + 1)
do
  sw=$(echo $hw | awk '{print 2*$1+1}')
  mkdir SlidingWindow_Frames_$sw
  for library in $(awk '{print $1}' ../libraries.txt)
  do
    for dir in fwd rev
    do
      for frame in 0 1 2
      do
        awk '{if ($2%3=="'"$frame"'") print $0}' HistFiles_wZeros/$library.$dir \
          > temp8 # $2 is the mapping position of the read; hence $2%3=frame
        awk '{print $1}' temp8 > temp9 #put the read counts in here
        $scriptsdir/slidingwindow temp9 $hw > temp10
        paste temp8 temp10 | awk '{print $4,$2}' > SlidingWindow_Frames_$sw/$library.$dir.$frame
        rm -f temp8 temp9 temp10
      done
    done
  done
done

#------------------------------------------------------------------------------

# Log scale genome coverage plots

mkdir Plots_Log

# change all instances of "Sample" below to the name of your virus, e.g. "HIV"

# NB: you must construct the "head.R" file for your virus before running this

for group in $(cat groups.txt)
do
  nlib=$(wc -l $group.txt | awk '{print $1}')
  cat $plotdir/Sample_head.R | sed 's/ggg/'$group'/' | \
    sed 's/nnn/'$nlib'/g' | sed 's/_xxx//'  > temp1
  count=0
  for line in $(awk '{printf "%s:%s\n", $1,$2}' $group.txt)
  do
    ((count+=1))
    library=$(echo $line | awk -F: '{print $1}')
    title=$(echo $line | awk -F: '{print $2}')
    cat $plotdir/log_middle.R | sed 's/ccc/'$count'/' | \
      sed 's/lll/'$library'/' | sed 's/ttt/'$title'/' | \
      sed 's/nnn/'$nlib'/g' | sed 's/ddd/HistFiles_wZeros/' >> temp1
  done
  echo "dev.off()" >> temp1
  R --no-save --slave < temp1
  mv -f Sample_log_$group.jpg Plots_Log/
  mv -f temp1 Plots_Log/Sample_log_$group.R
done


#------------------------------------------------------------------------------

#Make log Sliding Window histogram plots

mkdir Plots_Log_Sw
for hw in 1 7 15 150
do
    sw=$(echo $hw | awk '{print 2*$1+1}')
    for group in $(cat groups.txt)
    do
        nlib=$(wc -l $group.txt | awk '{print $1}')
        cat $plotdir/Sample_head.R | sed 's/ggg/'$group'/' | \
         sed 's/nnn/'$nlib'/g' | sed 's/xxx/sw'$sw'/' > temp1
    count=0
    for line in $(awk '{printf "%s:%s\n", $1,$2}' $group.txt)
    do
      ((count+=1))
      library=$(echo $line | awk -F: '{print $1}')
      title=$(echo $line | awk -F: '{print $2}')
      cat $plotdir/log_middle.R | sed 's/ccc/'$count'/' | \
        sed 's/lll/'$library'/' | sed 's/ttt/'$title'/' | \
        sed 's/nnn/'$nlib'/g' | sed 's/ddd/'SlidingWindow_$sw'/' >> temp1
    done
    echo "dev.off()" >> temp1 #"
    R --no-save --slave < temp1
    mv -f Sample_log_${group}_sw$sw.jpg Plots_Log_Sw/
    mv -f temp1 Plots_Log_Sw/Sample_log_${group}_sw$sw.R
  done
done

#------------------------------------------------------------------------------

#Make linear plots, for separate genome regions

# again, change "Sample" to your virus name

mkdir Plots_Linear_Sw
for hw in 1 7 15 150
do
    sw=$(echo $hw | awk '{print 2*$1+1}')
    for group in $(cat groups.txt)
    do
        nlib=$(wc -l $group.txt | awk '{print $1}')
        cat $plotdir/Sample_head_linear.R  | \
            sed 's/ggg/'$group'/' | sed 's/nnn/'$nlib'/g' | sed 's/xxx/sw'$sw'/' \
            > temp1
        count=0
        lastcount=0
        for line in $(awk '{printf "%s:%s\n", $1,$2}' $group.txt)
        do
            ((count += 1))
            if [ $count == $nlib ]
            then
                lastcount=1
            fi
            library=$(echo $line | awk -F: '{print $1}')
            title=$(echo $line | awk -F: '{print $2}')
            cat $plotdir/linear_middle_sw.R | sed 's/ccc/'$count'/' | \
            sed 's/lll/'$library'/' | \
            sed 's/ttt/'$title'/' | \
            sed 's/ddd/'SlidingWindow_$sw'/' | sed 's/zzz/'$lastcount'/' >> temp1
            #when the last library has been added to the R script, then zzz will be 1. otherwise it will be 0.
        done
        cat $plotdir/linear_end.R >> temp1
        R --no-save --slave < temp1
        mv -f Sample_${region}_${group}_sw$sw.jpg Plots_Linear_Sw/
        mv -f temp1 Plots_Linear_Sw/Sample_${region}_${group}_sw$sw.R
    done
done

#Make linear plots, for separate genome regions, divided by frame.
#  Note that, to make the R file easier to read, I removed the empty-file tests
#  from this script. Also removed the reverse reads. Also )grep -v RNASeq).
mkdir Plots_Linear_Frames_Sw
for hw in 1 7 15 150
do
    sw=$(echo $hw | awk '{print 2*$1+1}')
    for group in $(cat groups.txt | grep -v RNASeq)
    do
        for region in wholegenome
        do
            nlib=$(wc -l $group.txt | awk '{print $1}')
            cat $plotdir/Sample_head_linear.R  | \
            sed 's/ggg/'$group'_frames/' | sed 's/nnn/'$nlib'/g' | \
            sed 's/xxx/sw'$sw'/' > temp1
            count=0
            lastcount=0
            for line in $(awk '{printf "%s:%s\n", $1,$2}' $group.txt)
            do
                ((count += 1))
                if [ $count == $nlib ]
                then
                    lastcount=1
                fi
                library=$(echo $line | awk -F: '{print $1}')
                title=$(echo $line | awk -F: '{print $2}')
                cat $plotdir/linear_middle_frames.R | \
                    sed 's/ccc/'$count'/' | sed 's/lll/'$library'/' | \
                    sed 's/ttt/'$title'/' | \
                    sed 's/ddd/'SlidingWindow_Frames_$sw'/' | sed 's/zzz/'$lastcount'/' \
                    >> temp1
            done
            cat $plotdir/linear_end.R >> temp1
            R --no-save --slave < temp1
            mv -f Sample_${region}_${group}_frames_sw$sw.jpg Plots_Linear_Frames_Sw/
            mv -f temp1 Plots_Linear_Frames_Sw/Sample_${region}_${group}_frames_sw$sw.R
        done
    done
done


#------------------------------------------------------------------------------

#Normalized plots

# the easiest way to normalize is to just rescale ../{SlidingWindow_*/*,
#  SlidingWindow_Frames_*/*,HistFiles/*,HistFiles_wZeros/*} and then redo the
#  plots using the exact same scripts as above.

# create a file "normalizations.txt" which contains library names in column 1
# and total mapped reads (virus + refseq mRNA) in column 2

mkdir Normalized
cd Normalized
ln -s ../groups.txt
ln -s ../../libraries.txt

rm -f temp1
for file in $(cat groups.txt)
do
  ln -s ../$file.txt
  awk '{print $1}' $file.txt >> temp1
done

sort temp1 | uniq > libs.txt
rm -f temp1

# make normalizations file - using only +ve sense mapping reads from virus + host
echo "Library Total_mapped_reads" > normalizations.txt

for line in $(awk '{printf "%s:%s:%s:%s:%s:%s:%s\n", $1,$2,$3,$4,$5,$6,$7}' libraries.txt)
do
    library=$(echo $line | awk -F: '{print $1}')
    virus=$(echo $line | awk -F: '{print $2}')
    hostname=$(echo $line | awk -F: '{print $3}')
    condition=$(echo $line | awk -F: '{print $4}')
    libtype=$(echo $line | awk -F: '{print $5}')
    mrnareads=$(awk '{if ($3=="+") print $0}' \
      ../../$library.mRNA.bowtie | wc -l | awk '{print $1}')
    vrnareads=$(awk '{if ($3=="+") print $0}' \
      ../../$library.vRNA.bowtie | wc -l | awk '{print $1}')
    totalreads=$((mrnareads + vrnareads))
    echo "$library $totalreads" >> normalizations.txt
done

normfile="normalizations.txt"

mkdir HistFiles HistFiles_wZeros
for sw in 3 15 31 301
do
  mkdir SlidingWindow_$sw
  mkdir SlidingWindow_Frames_$sw
done

for line in $(awk '{printf "%s:%s:%s:%s\n", $1,$2,$3,$4}' libraries.txt)
do
  library=$(echo $line | awk -F: '{print $1}')
  scale=$(grep -w $library $normfile | awk '{print 1000000/$2}')
  for dir in fwd rev
  do
    awk '{print $1*"'"$scale"'",$2}' ../HistFiles/$library.$dir \
       > HistFiles/$library.$dir
    awk '{print $1*"'"$scale"'",$2}' ../HistFiles_wZeros/$library.$dir \
       > HistFiles_wZeros/$library.$dir
    for sw in 3 15 31 301
    do
      awk '{print $1*"'"$scale"'",$2}' ../SlidingWindow_$sw/$library.$dir \
        > SlidingWindow_$sw/$library.$dir
      for frame in 0 1 2
      do
        awk '{print $1*"'"$scale"'",$2}' \
          ../SlidingWindow_Frames_$sw/$library.$dir.$frame > \
           SlidingWindow_Frames_$sw/$library.$dir.$frame
      done
    done
  done
done

# now remake the plots as per above commands
