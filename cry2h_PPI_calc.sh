#!/bin/bash

#counting raw alignments:
cat $1 $2 | awk '{if($3!~/:/)print $3}' | sort | uniq -c | sort -rn > raw.alignments.totals

############PART 2#############
#quality filter reads##########
###############################

#filter_bowtie'd_reads
echo "aligning Fastqs"
cat $1  | awk '/XM:i:[0-2]/ && /NM:i:[0-2]/{if(((($2==0)&&(substr($6,0,2)~/(1|2)S/)&&($6~/[3-9][0-9]M|1[0-9][0-9])M/))||(($2==0)&&(substr($6,0,4)~/[3-9][0-9]M|1[0-9][0-9]M/)))||(($
cat $2  | awk '/XM:i:[0-2]/ && /NM:i:[0-2]/{if(((($2==0)&&(substr($6,0,2)~/(1|2)S/)&&($6~/[3-9][0-9]M|1[0-9][0-9])M/))||(($2==0)&&(substr($6,0,4)~/[3-9][0-9]M|1[0-9][0-9]M/)))||(($


#join R1 and R2 files
echo "joining R1 and R2"
awk 'NR==FNR {a[$1]=$1" "$2" "$3" "$4" "$6" "$10;next} $1 in a {print a[$1],$2,$3,$4,$6,$10}' R1.nomismatch.sorted R2.nomismatch.sorted> R1R2
wc -l R1R2> R1R2.wc

####################PART 3##########################
#DNA strandedness and fragment size filter##########
####################################################

echo "finding PPIs based on gene sequences performing fragment analysis"

cat $1 | awk -F ":" '{if($1 ~ "@SQ")print ($2" "$3)}' | sort -k1,1n > gene.length

awk 'NR==FNR {a[$1]=$2;next}{if(($2==$7)&&(((a[$3]-$4)+(a[$8]-$9)+116)>=220)&&(((a[$3]-$4)+(a[$8]-$9)+116)<=520))print $0}' gene.lengths R1R2.nomulti > frag.length.filtered.pairs

############PART 4#############
#totaling PPI fragments #######
###############################

echo "totaling reads IDing PPI and AD or DB fusions"

cat frag.length.filtered.pairs| awk '{print $3" "$8}'| awk '{if ($2>$1) print substr($2,0,9)" "substr($1,0,9); else print substr($1,0,9)" "substr($2,0,9)}' |sort -k1,1V -k2,2V | uniq -c | sort -n > frag.length.filtered.pai$

echo "separating PPIs from AD or DB fusions"
cat frag.length.filtered.pair_counts.csv | awk '{if(($0!~/AD/)&&($0!~/DB/))print $0}' > PPI_pairs.csv

echo "finding clones in PPIs"
awk '{if($3>$8)print $1,$3,$4,$8,$9; else print $1,$8,$9,$3,$4}' frag.length.filtered.pairs | sort -u -k2,2V -k4,4V -k3,3n -k5,5n | awk '{print $1}' > R1R2.PPIs.noclones.ids

echo "removing PPI clones"
awk 'NR==FNR {a[$1]=$1;next} $1 in a {print $3" "$8}' R1R2.PPIs.noclones.ids frag.length.filtered.pairs| awk '{if ($2>$1) print substr($2,0,9)" "substr($1,0,9); else print substr($1,0,9)" "substr($2,0,9)}' |sort -k1,1V -k2$

echo "separating non-clonal PPIs from fusions"
cat R1R2.PPIs.noclones | awk '{if(($0!~/AD/)&&($0!~/DB/))print $0}' > PPI_pairs.noclones.csv
