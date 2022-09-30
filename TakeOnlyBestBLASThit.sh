#!/bin/bash

#!/bin/bash

##where argv(1) are the blastp results in outputfmt 6 (tab delimited)##

export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr "$1" | sort -u -k1,1 --merge > "$1"_TopHitsOnly.txt
