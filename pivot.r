library(dplyr)
library(reshape2)

pivot = dcast(algallCODE_pfam_counts, X2 ~ X3, value.var = "X1")

write.table(pivot, file='pivot.tsv', quote=FALSE, sep='\t', row.names = FALSE)
