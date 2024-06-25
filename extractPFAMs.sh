#!/bin/bash

## just get PFxxxxx.xx accessions from hmmsearch results or similar data with PFAMs

sed -n 's/.*\(PF[0-9]\{5\}\.[0-9]\{2\}\).*/\1/p' $1 >> "$1"_PFs
