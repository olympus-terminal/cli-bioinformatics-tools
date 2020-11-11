#!/bin/bash

sed -i -E 's/B1_|B2_|B3_/0H_/g' $1
sed -i -E 's/B4_|B5_|B6_/2H_/g' $1
sed -i -E 's/B7_|B8_|B9_/4H_/g' $1
sed -i -E 's/B10_|B11_|B12_/24H_/g' $1
sed -i -E 's/B13_|B14_|B15_/48H_/g' $1
sed -i -E 's/B16_|B17_|B18_/72H_/g' $1

sed -i 's/H_.*TE\:/H/g' $1

sed -i.bak 's/H_.*TE/H/g' $1
