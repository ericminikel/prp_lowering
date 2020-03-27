options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/determinants/')
library(openxlsx)
library(reshape2)

bw_raw = read.xlsx('data_received/BLI_WeightData.xlsx',sheet='Sheet1')

bw1 = as.matrix(bw_raw[3:28,2:15])
rownames(bw1) = bw_raw[3:28,1]
colnames(bw1) = gsub('[^0-9]','',bw_raw[2,2:15])
bwdf1 = melt(bw1)
colnames(bwdf1) = c('dpi','animal','weight')
bwdf1$weight = as.numeric(bwdf1$weight)
bwdf1 = bwdf1[!is.na(bwdf1$weight),]

bw2 = as.matrix(bw_raw[32:52,2:10])
rownames(bw2) = bw_raw[32:52,1]
colnames(bw2) = gsub('.* ','',bw_raw[31,2:10])
bwdf2 = melt(bw2)
colnames(bwdf2) = c('dpi','animal','weight')
bwdf2$weight = as.numeric(bwdf2$weight)
bwdf2 = bwdf2[!is.na(bwdf2$weight),] 

bw3 = as.matrix(bw_raw[56:83,2:10])
rownames(bw3) = bw_raw[56:83,1]
colnames(bw3) = gsub('.* ','',bw_raw[55,2:10])
bwdf3 = melt(bw3)
colnames(bwdf3) = c('dpi','animal','weight')
bwdf3$weight = as.numeric(bwdf3$weight)
bwdf3 = bwdf3[!is.na(bwdf3$weight),] 

bw4 = as.matrix(bw_raw[87:105,2:10])
rownames(bw4) = bw_raw[87:105,1]
colnames(bw4) = gsub('.* ','',bw_raw[86,2:10])
bwdf4 = melt(bw4)
colnames(bwdf4) = c('dpi','animal','weight')
bwdf4$weight = as.numeric(bwdf4$weight)
bwdf4 = bwdf4[!is.na(bwdf4$weight),] 

bw5 = as.matrix(bw_raw[109:128,2:10])
rownames(bw5) = bw_raw[109:128,1]
colnames(bw5) = gsub('.* ','',trimws(bw_raw[108,2:10]))
bwdf5 = melt(bw5)
colnames(bwdf5) = c('dpi','animal','weight')
bwdf5$weight = as.numeric(bwdf5$weight)
bwdf5 = bwdf5[!is.na(bwdf5$weight),] 

bw = rbind(bwdf1, bwdf2, bwdf3, bwdf4, bwdf5)

# correct typos in animal IDs
bw$animal[bw$animal == '55771'] = '55771/2'
bw$animal[bw$animal == '5860'] = '55860'

write.table(bw, 'data/mri/B_weights.tsv', sep='\t', col.names=T, row.names=F, quote=F)
