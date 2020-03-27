options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/determinants/')
library(sqldf)
library(openxlsx)
library(reshape2)

raw = read.xlsx('data_received/Nat.Hist.Wts.WORKING.xlsx',sheet='Sheet1')

wts1mat = raw[5:28, 2:12]
colnames(wts1mat) = c('animal',raw[4,3:12])
wts1 = melt(wts1mat, id.vars='animal')
colnames(wts1) = c('animal','dpi','wt')

wts2 = raw[5:28, c(2, 15, 13)]
colnames(wts2) = c('animal','dpi','wt')

wts = rbind(wts1, wts2)
wts = wts[!is.na(wts$wt),]
wts = wts[!duplicated(paste0(wts$dpi, '-', wts$animal)),]
wts = wts[with(wts, order(dpi, animal)),]

write.table(wts, 'data/mri/N_weights.tsv', sep='\t', col.names=T, row.names=F, quote=F)