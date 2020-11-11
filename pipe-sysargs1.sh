#!/bin/bash

##from target gene list (sysargv1) to extract isoforms and plot with custom plotly boxplots##



cat $1 | while read l ; do grep $l DrugTargetIsoforms.csv > "$l".csv; sed -i '1s;^;Time,Gene,Isoform,TPM\n;' "$l".csv; done
cat $1 | while read l; do python ./boxout-custom.px "$l".csv "$l" "$l".pdf; done
