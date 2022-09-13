#!/bin/bash

cut -d [ -f 2 $1 > "$1".1 ;
cut -d ' ' -f 1 "$1".1 > "$1".2 ;
sort "$1".2 | uniq -c - > "$1"_genus_tallies.txt
rm "$1".1
rm "$1".2
