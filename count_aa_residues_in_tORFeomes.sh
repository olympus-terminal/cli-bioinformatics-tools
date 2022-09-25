#!/bin/bash

for f in *.aa.fa; do echo "$f" >1; awk 'NR>1{s+=length()} END{print s}' $f > 2 ;  paste 1 2 >>aa-counts.txt ; done
