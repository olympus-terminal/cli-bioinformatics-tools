#!/bin/bash

mkdir TAGGED
mkdir TABBED
mkdir PFAMs

for f in *.seqtblout
do
 j=$(basename $f .seqtlout)
 sed 's/ \+/\t/g' $f > tabbed-"$f"
 sed "s/^/$j	/" tabbed-"$f">tagged-"$f"
 #rm tabbed-"$f"
 cut -f 1,5,9 tagged-"$f" >"$f".p
 #rm tagged-"$f"
 #tail -n +4 "$f".p | head -n-10 > "$f".pfams
 #rm "$f".p
done

mv tabbed* TABBED/
mv tagged* TAGGED/
mv *.p PFAMs/
