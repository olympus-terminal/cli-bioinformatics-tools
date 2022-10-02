#/bin/bash

cat $1 | while read l; 

do grep "$l" $2 > 1;  ##should be unnecessary if one list (one file) is used for variable references

div=$(cut -f 2 1);
assembly=$(cut -f 1 1);
sed -iE "s|$assembly|$div|" $3; done

##Where argv(1) is the match in the file you want to translate, argv(2) is a two column file having the translations, and argv(3) is file containing the targets to be replaced (e.g., in this case it was a phylogenetic tree)
