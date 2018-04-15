#!/bin/bash

while read l
do
 echo "$l" > tally
 grep -m 1 $l pfam.convert >>tmp1
 paste tally tmp1 >>record
 rm tally
 rm tmp1
done <names
