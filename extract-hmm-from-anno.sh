
for f in *-labels.txt; do grep -w -f $f Human.anno > "$f".anno; done
for f in *.anno; do cut -d , -f 9 $f > "$f".proteins.txt; done
for f in *proteins.txt; do sed -i 's/"//g' $f; done
for f in *proteins.txt; do grep -w -f $f hmmsearch-results.domtblout > "$f"-vfams.txt; done
