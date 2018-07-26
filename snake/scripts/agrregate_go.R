f=c('~/tmp/go1','~/tmp/go2','~/tmp/go3')
library(tidyverse)
go=data.frame(goa=character())
for (afile in f){
	data = read_tsv(afile, col_names = c('goa',basename(afile)))
	go=full_join(go, data, by='goa')
}
got = as.data.frame(t(go[,-1])) %>%
setNames(go[[1]])
write_tsv(rownames_to_column(got, var='sample'), {snakemake@output[["raw"]]})