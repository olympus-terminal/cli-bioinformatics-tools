#!/bin/bash

for f in tabbed-GCF_00*;

do cat $f | while read l;
            do echo $l > tmp1;
            cut -d ' ' -f 5 tmp1 > tmp2;
            cut -d '.' -f 1 tmp2 > tmp3;
            grep -f tmp3 EC-Pfam_calculated_associations_Extended.csv > tmp4;
            paste tmp1 tmp4 >> "$f".ECs-out.txt;
            done;
done
