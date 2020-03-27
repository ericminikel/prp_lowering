options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/determinants/')
library(sqldf)
library(openxlsx)
library(reshape2)

template = read.xlsx('data_received/IONIS behavior Natural History.xlsx',sheet='ALL MICE TEMPLATE',startRow=4)
colnames(template) = gsub('[^0-9a-z_]','_',tolower(colnames(template)))
template$category = tolower(c(rep(template$x1[1],21), rep(template$x1[22], 40-22+1), rep(template$x1[41], 65-41+1)))
# check that categories were assigned correctly
# View(template[,c('x1','category','observations','scoring')])
use_rows = !grepl('EXCLUDE',template$scoring) & !is.na(template$observations) & template$observations != 'COMMENTS'
# check that worked correctly:
# template$observations[use_rows]
template$observations = tolower(template$observations)
template$scoring = tolower(template$scoring)
behavs_meta = template[use_rows,c('category','observations','scoring')]

# write out initial draft - went in and manually updated later
# write.table(behavs_meta, 'data/mri/N_behavs_meta.tsv', sep='\t', col.names=T, row.names=F, quote=F)
behavs_meta = read.table('data/mri/N_behavs_meta.tsv', sep='\t', header=T, quote='', comment.char='')

sheetnames = getSheetNames('data_received/IONIS behavior Natural History.xlsx')
behavs_dpi = as.integer(gsub('[^0-9]','',sheetnames)) # note first is NA - leaving it in there so indices match up in loop

behavs = data.frame(dpi=integer(0), animal=integer(0), observation=character(0), score=numeric(0))
for (i in 2:length(sheetnames)) {
  raw = read.xlsx('data_received/IONIS behavior Natural History.xlsx',sheet=sheetnames[i],startRow=4)
  mat = raw[which(use_rows),c(2,4:ncol(raw))]
  mat[,1] = tolower(mat[,1])
  colnames(mat) = suppressWarnings(as.integer(colnames(mat)))
  #rownames(mat) = tolower(raw[which(use_rows),2])
  melted = melt(mat, id.vars=1)
  colnames(melted) = c('observation','animal','score')
  melted$dpi = behavs_dpi[i]
  behavs = rbind(behavs, melted[,colnames(behavs)])
}

surviv = read.table('data/mri/N_survival.tsv',sep='\t',header=T)

behavs_cleaned = sqldf("
select   b.dpi, b.animal, b.observation, b.score
from     behavs b, surviv s, behavs_meta bm
where    b.animal = s.animal
and      b.dpi <= s.dpi
and      b.observation = bm.observations
order by 1, 2, 3
;")

write.table(behavs_cleaned, 'data/mri/N_behavs.tsv', sep='\t', col.names=T, row.names=F, quote=F)