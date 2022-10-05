#!/bin/bash

mkdir TAGGED
mkdir TABBED
mkdir VFAMs

for f in *.seqtblout
do
 j=$(basename $f .seqtlout)
 sed 's/ \+/\t/g' $f > tabbed-"$f"
 sed "s/^/$j	/" tabbed-"$f">tagged-"$f"
 #rm tabbed-"$f"
 cut -f 1,5,9 tagged-"$f" >"$f".v
 #rm tagged-"$f"
 #tail -n +4 "$f".v | head -n-10 > "$f".vfams
 #rm "$f".v
done

mv tabbed* TABBED/
mv tagged* TAGGED/
mv *.v VFAMs/
