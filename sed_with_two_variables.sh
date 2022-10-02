#/bin/bash

cat assemblies.txt | while read l; do grep "$l" assembly-div.txt > 1; div=$(cut -f 2 1); assembly=$(cut -f 1 1); sed -iE "s|$assembly|$div|" SpeciesTree_rooted_node_labels.txt; done
