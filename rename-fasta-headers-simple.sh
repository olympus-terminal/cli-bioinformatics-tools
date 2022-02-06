#!/bin/bash

for f in *.fa;

do

awk '/^>/{gsub(/^>/,">Seq"i++" ");}1' $f > awked-"$f";

sed 's/ .*//g' awked-"$f" > renamed-"$f";

rm awked-"$f";

mkdir renamed-fa;

mv renamed-"$f" renamed-fa/

done
