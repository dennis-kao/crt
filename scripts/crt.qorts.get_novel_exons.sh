#!/bin/bash

#parses qorts output looking for novel junctions in muscular genes and novel exons like found in Gonorazky.2016.DMD article
#$1 = sample name for output: muscle2
#$2 = threshold for novel exons discovery, the maximum distance between novel junctions to detect an exon in bp = 100, 200
# the coverage threshold for junction is 5

gunzip -c QC.spliceJunctionCounts.novelSplices.txt.gz | sed 1d | cut -f 1,3,4,5 > $1.novel_junctions_all.bed

#remove novel junctions found in 3 MH controls
#it has been used for muscle samples, however for MH samples themselves we don't have to remove junctions
#bedtools intersect is bad because it removes everything
#cat $1.novel_junctions_all.bed | awk '{print $1"-"$2"-"$3}' | sort > $1.novel_junctions_all.txt
#comm -23 $1.novel_junctions_all.txt ~/work/project_muscular/mh_novel_splices.txt  > $1.novel_junctions.txt

#cat $1.novel_junctions_all.bed | awk '{print $1"-"$2"-"$3"\t"$4}' | sort > $1.novel_junctions_all.txt
#rnaseq.qorts.sublist.pl $1.novel_junctions_all.txt $1.novel_junctions.txt | sed s/"-"/"\t"/g > $1.novel_junctions.bed

cp $1.novel_junctions_all.bed $1.novel_junctions.bed

#5 reads is a coverage threshold for novel junction
bedtools intersect -wa -wb -a $1.novel_junctions.bed -b ~/Desktop/reference_tables/muscular_genes.bed | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4}' | awk '{if ($5>5) print $0}' > $1.results.bed

#get novel exons $2 = 100 or 200
cat $1.results.bed | awk 'BEGIN{prev=0}{print $0"\t"$2-prev;prev=$3}' | awk -v threshold=$2 '{if ($6>0 && $6<threshold) print $0;}' > $1.novel_exons

#report:
echo "Junctions:" `cat $1.results.bed | wc -l`
echo "Exons:" `cat $1.novel_exons | wc -l`

gunzip -c QC.spliceJunctionCounts.knownSplices.txt.gz | sed 1d | cut -f 2,4,5,6 > $1.known_junctions.bed
bedtools intersect -wa -wb -a $1.known_junctions.bed -b ~/Desktop/reference_tables/muscular_genes.bed | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4}' | awk '{if ($5<2) print $0}' > $1.known.low_coverage