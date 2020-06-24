options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/prp_lowering/')
library(sqldf)
library(reshape2)

#### BEHAVIORALS

# list of phenotypes
signs = data.frame(name=c('scruff_poor_grooming','poor_body_condition','reduced_activity','hunched_posture','irregular_gait_hindlimb_weakness','tremor','blank_stare','difficulty_righting'))

# master data frame for all processed behaviorals
behavs_all = data.frame(animal=character(0), dpi=integer(0), phenotype=character(0), value=integer(0), rater=character(0))

expts = c('D','G','M','Q','R','S','T','V','Z')

for (expt in expts) {

raw = read.table(paste('data/raw/',expt,'_behaviorals_raw.tsv',sep=''),sep='\t',quote='',comment.char='')

raters = data.frame(date=character(0), rater=character(0))

n = 0 # prod(dim(raw)) if you wanted to pre-allocate all the memory instead of rbind()
behavs = data.frame(animal_no=character(n), date=character(n), phenotype=character(n), value=character(n), rater=character(n))

for (i in 4:nrow(raw)) {
  animal = raw[i,1]
  
  for (j in 1:ncol(raw)) {
    if (raw[3,j] %in% signs$name) { # if this is a behavioral 0/1 column
      phenotype = raw[3,j]
      date = substr(raw[2,j],1,10)
      rater = gsub('.* ','',raw[2,j])
      value = raw[i,j]
      cat(paste("\rProcessing animal ",animal,", date ",date,"..."))
      flush.console()
      behavs = rbind(behavs, cbind(animal, date, phenotype, value, rater))
    }
  }
}

colnames(behavs) = c('animal','date','phenotype','value','rater')

# now join to master and remove observations after death date - results from some typos in the date field
master = read.table(paste('data/raw/',expt,'_master.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

behavs_cleaned = sqldf("
select   b.date, b.animal, julianday(b.date) - julianday(m.inoculation_date) dpi, b.phenotype, case when b.value = '1' then 1 else 0 end value, b.rater
from     behavs b, master m
where    b.animal = m.animal
and      b.date is not null
and      b.date != ''
and      b.date < m.death_date
order by 1, 2, 3, 4
;")

behavs_all = rbind(behavs_all, behavs_cleaned[,c('animal','dpi','phenotype','value','rater')])

}

# now integrate the 01, L, and P experiments in their nonstandard data format
raw = read.table('data/raw/01LP_behaviorals_raw.tsv',sep='\t',quote='',comment.char='')
master = read.table('data/raw/01LP_master.tsv',sep='\t',header=T,quote='',comment.char='')

raters = data.frame(date=character(0), rater=character(0))

# propagate dates throughout
last_date = ''
for (i in 3:ncol(raw)) {
  # skip NA cols
  if (is.na(raw[1,i])) {
    next
  }
  # read in new date for a set of cols
  if (raw[1,i] != '') {
    last_date = substr(raw[1,i], 1, 10)
    rater = trimws(substr(raw[1,i],11,20),which='both') # grab the initials of who did the behaviorals
    if (last_date == '') {
      next
    } else {
      if (grepl('/',last_date)) {
        date = as.Date(last_date, format='%m/%d/%y')
      } else {
        date = as.Date(last_date, format='%Y-%m-%d')
      }
      raters = rbind(raters, c(as.character(date),rater))
    }
    next
    # apply that date to the rest of that set
  } else {
    raw[1,i] = last_date
  }
}

colnames(raters) = c('date','rater')

n = 0 # prod(dim(raw)) if you wanted to pre-allocate all the memory instead of rbind()
behavs = data.frame(animal_no=character(n), date=character(n), phenotype=character(n), value=character(n))
# manually melt the table, removing observations without comments as these are often blank meaning not done, rather than blank meaning zero
# this takes several minutes to run
i = 1
for (row in 3:nrow(raw)) {
  animal_no = raw[row,1]
  
  flush.console()
  for (col in 3:1204) { # should be ncol(raw) instead of 1204 but last date is empty
    if (is.na(raw[row,col])) {
      next # skip NA cols
    }
    if (raw[2,col] == 'weight change') {
      next # skip weight change formulae that were manually entered on occasion
    }
    
    raw_date = raw[1,col]
    if (raw_date == '') {
      next
    } else {
      if (grepl('/',raw_date)) {
        date = as.Date(raw_date, format='%m/%d/%y')
      } else {
        date = as.Date(raw_date, format='%Y-%m-%d')
      }
    }
    cat(paste("\rProcessing animal ",animal_no,", date ",date,"..."))
    phenotype = raw[2,col]
    if (raw[row,col] == 1) {
      value = 1
    } else {
      value = 0
    }
    behavs = rbind(behavs, cbind(animal_no, as.character(date), phenotype, value))
  }
}
colnames(behavs) = c('animal','date','phenotype','value')
behavs$rater = raters$rater[match(behavs$date, raters$date)]

master$expt = ''
master$expt[1:84] = '01'
master$expt[85:132] = 'L'
master$expt[133:180] = 'P'
behavs$expt = master$expt[match(behavs$animal, master$animal_no)]

behavs$datetext = behavs$date
behavs$date = as.Date(behavs$datetext)

behavs$inoculation_date = as.Date(master$inoculation_date[match(behavs$animal, master$animal_no)])
behavs$dpi = as.integer(behavs$date - behavs$inoculation_date)

master$death_date = as.Date(master$dod)

behavs_cleaned = sqldf("
select   b.date, b.animal, b.dpi, b.phenotype, case when b.value = '1' then 1 else 0 end value, b.rater
from     behavs b, master m
where    b.animal = m.animal_no
and      b.date is not null
and      b.date != ''
and      b.date < m.death_date
and      (b.expt != '01' or b.dpi not in (139, 160, 166, 170, 181, 190, 255)) -- incomplete monitoring dates for cohort 01 identified by visual inspection
and      (b.expt != 'P' or b.dpi not in (140)) -- incomplete monitoring dates for cohort P identified by visual inspection
and      b.rater not in ('JLB','EVM','PM') -- JLB's behaviorals are systematically different from others; my (EVM) and afternoon (PM) checks are always incomplete
order by 1, 2, 3, 4
;")

behavs_all = rbind(behavs_all, behavs_cleaned[,c('animal','dpi','phenotype','value','rater')])

# handle unspecified raters (blank and dates filled in instead of initials)
behavs_all$rater[nchar(behavs_all$rater) == 0 | nchar(behavs_all$rater) == 10] = 'unknown'
# behavs_all = read.table('data/processed/behavs.tsv',sep='\t',header=T)

write.table(behavs_all,'data/processed/behavs.tsv',sep='\t',row.names=F,col.names=T,quote=F)







#### WEIGHTS

weights_all = data.frame(animal=character(0), dpi=integer(0), weight=numeric(0))

for (expt in expts) {

raw_weights = read.table(paste('data/raw/',expt,'_weights.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
surgery_notes = read.table(paste('data/raw/',expt,'_surgery_notes.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

melted_weights = melt(raw_weights,id.vars=c('animal'))
colnames(melted_weights) = c('animal','datetext','weight')
melted_weights$weight = as.numeric(melted_weights$weight) # this will force text, empty, and NA fields to NA
melted_weights = melted_weights[!is.na(melted_weights$weight),]
melted_weights$date = as.Date(gsub('X','',melted_weights$datetext),format='%Y.%m.%d')
melted_weights$inoculation_date = as.Date(surgery_notes$inoculation_date[match(melted_weights$animal, surgery_notes$animal)],format='%Y-%m-%d')
melted_weights$dpi = as.integer(melted_weights$date - melted_weights$inoculation_date)

melted_weights = melted_weights[!is.na(melted_weights$dpi),]

weights_all = rbind(weights_all, melted_weights[,c('animal','dpi','weight')])

}

# now integrate 01, L, and P expts in their nonstandard format
raw_weights = read.table(paste('data/raw/01LP_weights.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
raw_weights = raw_weights[,-which(colnames(raw_weights) %in% c('group_name','cage','treatment','baseline'))]
master = read.table('data/raw/01LP_master.tsv',sep='\t',header=T,quote='',comment.char='')

melted_weights = melt(raw_weights,id.vars=c('animal_no'))
colnames(melted_weights) = c('animal','datetext','weight')
melted_weights$weight = as.numeric(melted_weights$weight) # this will force text, empty, and NA fields to NA
melted_weights = melted_weights[!is.na(melted_weights$weight),]
melted_weights$date = as.Date(gsub('X','',melted_weights$datetext),format='%Y.%m.%d')
melted_weights = melted_weights[!is.na(melted_weights$date),]
melted_weights$inoculation_date = as.Date(master$inoculation_date[match(melted_weights$animal, master$animal_no)],format='%Y-%m-%d')
melted_weights$dpi = as.integer(melted_weights$date - melted_weights$inoculation_date)
melted_weights = melted_weights[!is.na(melted_weights$dpi) & melted_weights$weight > 0,] # remove NA dates/dpi, and the random zeroes somebody entered at 263 dpi for P cohort

weights_all = rbind(weights_all, melted_weights[,c('animal','dpi','weight')])

write.table(weights_all,paste('data/processed/weights.tsv',sep=''),sep='\t',row.names=F,col.names=T,quote=F)

#### NESTS

nst_all = data.frame(cage=character(0), dpi=integer(0), edry=numeric(0), nslt=numeric(0), comb=numeric(0))

for (expt in expts) {

# for some reason I get "incomplete final line found" warnings for some even though all the data are read in properly - use suppressWarnings to ignore these
nests = suppressWarnings(read.table(paste('data/raw/',expt,'_nests.tsv',sep=''),sep='\t',header=T,quote='',comment.char=''))
master = read.table(paste('data/raw/',expt,'_master.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

nst = melt(nests, id.vars=c('cage'))
# remove entries where the techs noted that either the cage had just been changed or was damp, making it impossible to score:
nst = nst[!grepl('([Cc]hange|[Dd]amp)',nst$variable),]
# parse dates and numeric ratings. suppress warnings - this will throw out invalid dates and comments the techs entered
nst$date = suppressWarnings(as.Date(gsub('X','',nst$variable),format='%Y.%m.%d'))
nst$edry =  suppressWarnings(as.numeric(gsub('\\/.*','',nst$value)))
nst$nslt =  suppressWarnings(as.numeric(gsub('.*\\/','',nst$value)))
nst$comb = (nst$edry + nst$nslt) / 2



# map dpi from master table. note although some cohorts have mixed treatment cages, no cages have different inoculation dates.
nst$inoculation_date = as.Date(master$inoculation_date[match(nst$cage,master$cage)])
nst$dpi = as.integer(nst$date - nst$inoculation_date)

nst = nst[!is.na(nst$comb) & !is.na(nst$dpi), c('cage','dpi','edry','nslt','comb')] # now remove missing data

nst_all = rbind(nst, nst_all)
}

write.table(nst_all, 'data/processed/nests.tsv', sep='\t', col.names=T, row.names=F, quote=F)


### SURVIVAL

# goals:
# - one row per animal
# columns:
#    - endpoint dpi
#    - endpoint type (from defined values = five symptoms, weight loss, fd, moribund, etc.)
#    - acm
#    - tissues

survival = data.frame(animal=character(0), dpi=integer(0), 
                      endpoint=logical(0), endpoint_comments=character(0), acm=logical(0), acm_comments=character(0), 
                      tissues=character(0))

for (expt in expts) {
  master = read.table(paste('data/raw/',expt,'_master.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
  
  master$dpi = as.integer(as.Date(master$death_date) - as.Date(master$inoculation_date))
  master$endpoint = master$prion_endpoint
  master$acm = master$acm1
  
  if (!('tissues' %in% colnames(master))) {
    master$tissues = master$tissues_harvested
  }
  
  survival = rbind(survival, master[,c('animal','dpi','endpoint','endpoint_comments','acm','acm_comments','tissues')])
}

# add in non-standard 01LP data
master = read.table('data/raw/01LP_master.tsv',sep='\t',header=T,quote='',comment.char='')
# note that the definitions of 'include' and 'acm1' are slightly different for the L cohort than for others.
# they are included here for convenience later on but note that they should be handled differently in analysis
master$dpi = as.integer(as.Date(master$dod) - as.Date(master$inoculation_date))
master$endpoint = master$include
master$acm = master$acm1
master$tissues = master$tissues_harvested
master$animal = master$animal_no
master = subset(master, nchar(animal_no) <= 3) # remove the "practice" animals without cage letters assigned. CC numbers are longer.
survival = rbind(survival, master[,c('animal','dpi','endpoint','endpoint_comments','acm','acm_comments','tissues')])

# bug fix in response to reviewer #3 - 2020-06-11
# for analyses restricted to animals reaching pre-specified euthanasia endpoint, "found dead" should be excluded
survival$endpoint[survival$endpoint_comments=='found dead'] = F

write.table(survival, 'data/processed/survival.tsv', sep='\t', col.names=T, row.names=F, quote=F)

# table(survival$endpoint_comments)
# table(survival$acm_comments)


#### SURGERY NOTES

surgery = data.frame(animal=character(0), dpi=integer(0), event=character(0))

for (expt in expts[expts != 'G']) { # exclude G since they had no ICVs
  surgery_notes = read.table(paste('data/raw/',expt,'_surgery_notes.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
  relevant_columns = c('animal',colnames(surgery_notes)[grepl('_date',colnames(surgery_notes))])
  surgery_melted = melt(surgery_notes[,relevant_columns], id.vars='animal')
  surgery_melted$date = as.Date(surgery_melted$value)
  surgery_melted$event = gsub('_date','',surgery_melted$variable)
  
  surgery_dpi = sqldf("
  select   icv.animal, icv.date - inoc.date dpi, icv.event
  from     surgery_melted icv, surgery_melted inoc
  where    icv.animal = inoc.animal
  and      icv.event != 'inoculation'
  and      inoc.event = 'inoculation'
  and      icv.date is not null -- remove surgeries that did not occur, e.g. because mouse had already died by that point
  ;")
  surgery = rbind(surgery, surgery_dpi)
}

# handle 01LP
surgery_notes = read.table('data/raw/01LP_master.tsv',sep='\t',header=T,quote='',comment.char='')
relevant_columns = c("animal_no", "icv_dose_1_date", "inoculation_date", "icv_dose_2_date")
surgery_melted = melt(surgery_notes[,relevant_columns], id.vars='animal_no')
surgery_melted$date = as.Date(surgery_melted$value)
surgery_melted$event = gsub('_dose_','',gsub('_date','',surgery_melted$variable))
surgery_melted = surgery_melted[nchar(surgery_melted$animal_no) <= 3,]

surgery_dpi = sqldf("
  select   icv.animal_no animal, icv.date - inoc.date dpi, icv.event
  from     surgery_melted icv, surgery_melted inoc
  where    icv.animal_no = inoc.animal_no
  and      icv.event != 'inoculation'
  and      inoc.event = 'inoculation'
  and      icv.date is not null -- remove surgeries that did not occur, e.g. because mouse had already died by that point
  ;")
surgery = rbind(surgery, surgery_dpi)

write.table(surgery, 'data/processed/surgery.tsv', sep='\t', row.names=F, col.names=T, quote=F)

#### BLIND & COHORT INFO

blind = read.table('data/raw/blind.tsv', sep='\t', header=T)
blind = blind[blind$expt != 'HK',] 


g_cohort_sex = read.table('data/raw/G_sex.tsv',sep='\t',header=T)
blind$sex = g_cohort_sex$sex[match(blind$animal, g_cohort_sex$animal)]
blind$sex[is.na(blind$sex)] = 'f'

# export cohort list for annotation
# sqldf("select expt, cohort, count(*) n from blind group by 1, 2 order by 1, 2;")

display_params = read.table('data/raw/display_params.tsv', sep='\t', header=T)
cohort_params = read.table('data/raw/cohort_params.tsv', sep='\t', header=T, comment.char='')

blind$aso = display_params$display[match(blind$treatment,display_params$private_name)]
blind$aso[is.na(blind$aso)] = 'no treatment'
blind$treatment = blind$aso

blind$strain = display_params$display[match(blind$inoculum, display_params$private_name)]
blind$strain[is.na(blind$strain)] = blind$inoculum[is.na(blind$strain)]

tpts = sqldf("select animal, min(dpi) timepoint from surgery group by 1 order by 1;")
blind$timepoint = tpts$timepoint[match(blind$animal, tpts$animal)]

blind$display = cohort_params$display[match(blind$cohort, cohort_params$cohort)]
blind$color   =   cohort_params$color[match(blind$cohort, cohort_params$cohort)]
blind$lty     =     cohort_params$lty[match(blind$cohort, cohort_params$cohort)]
blind$x       =       cohort_params$x[match(blind$cohort, cohort_params$cohort)]

blind_output = blind[,c('expt','cohort_alias','animal','sex','cage','strain','treatment','dose','timepoint','display','color','lty','x')]
colnames(blind_output)[2] = 'cohort'

write.table(blind_output, 'data/processed/blind.tsv', sep='\t', col.names=T, row.names=F, quote=F)

cohort_params_redacted = cohort_params[cohort_params$expt != 'HK',c('expt','cohort_alias','n','display','color','lty','x','color2','lty2','display2')]
colnames(cohort_params_redacted)[2] = 'cohort'

write.table(cohort_params_redacted, 'data/processed/params.tsv', sep='\t', col.names=T, row.names=F, quote=F)

