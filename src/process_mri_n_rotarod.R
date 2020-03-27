options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/determinants/')
library(sqldf)
library(openxlsx)
library(reshape2)

learning_runs = c('trial_1','trial_2','trial_3') # per deb's instructions, omit trials 1-3 out of 9
# as the mice are learning the task during these trials

rr = data.frame(dpi=integer(0), animal=integer(0), trial=character(0), latency=numeric(0))
for (dpi in c(-7, 30, 60, 90, 120, 150, 180)) {
  raw = read.xlsx('data_received/N_rotarod_merged.xlsx',sheet=as.character(dpi))
  colnames(raw) = gsub('[^a-z0-9_]','_',tolower(colnames(raw)))
  rr_temp = melt(raw, id.var='mouse')
  rr_temp$dpi = dpi
  rr_temp = rr_temp[,c(4,1,2,3)]
  colnames(rr_temp) = c('dpi','animal','trial','latency')
  rr_temp = rr_temp[!is.na(rr_temp$latency) & !(rr_temp$trial %in% learning_runs),]
  rr = rbind(rr, rr_temp)
}

master = read.table('data/mri/N_survival.tsv',sep='\t',header=T)

rr = rr[rr$animal %in% master$animal,]

rr = rr[with(rr, order(dpi, animal)),]

write.table(rr, 'data/mri/N_rotarod.tsv', sep='\t', col.names=T, row.names=F, quote=F)