for f in *.csv; do awk '{a[$2,$3]+=$1}END{for(i in a) print i,a[i]}' $f >> "$f"_consolidated.txt; done
