options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/determinants/')
library(reshape2)

c_wts_1 = read.table('data_received/C_weights_working1.tsv',sep='\t',header=T)
c_wts = melt(c_wts_1, id.vars = 'X')
c_wts$variable = as.integer(gsub('^X','',c_wts$variable))
colnames(c_wts) = c('animal','dpi','weight')
c_wts = c_wts[!is.na(c_wts$weight),]

c_wts_2 = read.table('data_received/C_weights_working2.tsv',sep='\t',header=T)
c_wts_3 = read.table('data_received/C_weights_working3.tsv',sep='\t',header=T)

c_wts = rbind(c_wts, c_wts_2, c_wts_3)


c_wts = c_wts[with(c_wts, order(dpi, animal)),]

write.table(c_wts, 'data/mri/C_weights.tsv', sep='\t', row.names=F, quote=F, col.names=T)
