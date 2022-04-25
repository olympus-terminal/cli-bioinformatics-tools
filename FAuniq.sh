#!/bin/bash

##originally from Pierre Lindenbaum on Biostars.org

sed -e '/^>/s/$/@/' -e 's/^>/#/' $1  |\
tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
sort -u -t '	' -f -k 2,2  |\
sed -e 's/^/>/' -e 's/\t/\n/' > $2
