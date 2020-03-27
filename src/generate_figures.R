options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/determinants/')
library(sqldf)
library(survival)
library(mratios)
library(drc)

imgmode = 'png'
# imgmode = 'pdf'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}

expand_range = function(x,by=0.5) {
  xmin = min(x,na.rm=T) - by
  xmax = max(x,na.rm=T) + by
  return (c(xmin, xmax))
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

calculate_ci = function(df) {
  mean_col = colnames(df)[grepl('^mean_',colnames(df))]
  sd_col = colnames(df)[grepl('^sd_',colnames(df))]
  n_col = colnames(df)[grepl('^n_',colnames(df))]
  df$l95 = df[,mean_col] - 1.96*df[,sd_col]/sqrt(df[,n_col])
  df$u95 = df[,mean_col] + 1.96*df[,sd_col]/sqrt(df[,n_col])
  return (df)
}

default_lwd = 3
strain_lwd = 1.5

params = read.table('data/processed/params.tsv',sep='\t',header=T,comment.char='')
blind = read.table('data/processed/blind.tsv',sep='\t',header=T,comment.char='')
surviv = read.table('data/processed/survival.tsv',sep='\t',header=T,comment.char='')
surgery = read.table('data/processed/surgery.tsv',sep='\t',header=T,comment.char='')
behavs = read.table('data/processed/behavs.tsv',sep='\t',header=T,comment.char='')
nests = read.table('data/processed/nests.tsv',sep='\t',header=T,comment.char='')
weights = read.table('data/processed/weights.tsv',sep='\t',header=T,comment.char='')

#### GROUP BY COHORT AND CALCULATE MEAN WEIGHTS AND BEHAVS

rm(v)
rm(b)
v1 = sqldf("
           select   v.animal, v.dpi, sum(value) n_signs
           from     behavs v
           group by 1, 2
           order by 1, 2
           ;")
v2 = sqldf("
           select   b.expt, b.cohort, v1.dpi, avg(n_signs) mean_signs, stdev(n_signs) sd_signs, count(*) n_animals
           from     v1, blind b
           where    b.animal = v1.animal
           group by 1, 2, 3
           having   count(*) > 1
           order by 1, 2, 3
           ;")
v = calculate_ci(v2)
v$color = params$color2[match(v$cohort, params$cohort)]
v$lty = params$lty2[match(v$cohort, params$cohort)]


w = sqldf("
          select   b.expt, b.cohort, w.dpi, avg(w.weight) mean_weight, stdev(w.weight) sd_weight, count(*) n_weight
          from     weights w, blind b
          where    w.animal = b.animal
          group by 1, 2, 3
          having   count(*) > 1
          order by 1, 2, 3
          ;")
w = calculate_ci(w)
w$color = params$color[match(w$cohort, params$cohort)]
w$lty = params$lty2[match(w$cohort, params$cohort)]

# date for first weights
fwdpi = sqldf("
              select   animal, min(dpi) first_dpi
              from     weights
              where    dpi >= 80
              group by 1
              order by 1
              ;")
# first weights
fw = sqldf("
           select   f.animal, f.first_dpi, w.weight
           from     fwdpi f, weights w
           where    f.animal = w.animal
           and      f.first_dpi = w.dpi
           order by 1
           ;")
rw1 = sqldf("
           select   w.animal, w.dpi, w.weight / f.weight relwt
           from     weights w, fw f
           where    w.animal = f.animal
           order by 1, 2
           ;")


rw2 = sqldf("
          select   b.expt, b.cohort, w.dpi, avg(w.relwt) mean_weight, stdev(w.relwt) sd_weight, count(*) n_weight
          from     rw1 w, blind b
          where    w.animal = b.animal
          and      dpi > 80
          group by 1, 2, 3
          having   count(*) > 1
          order by 1, 2, 3
          ;")
rw = calculate_ci(rw2)
rw$color = params$color2[match(rw$cohort, params$cohort)]
rw$lty = params$lty2[match(rw$cohort, params$cohort)]


#### PREPARE NEST DATA
rm(b)
rm(nc)
rm(n)
# remove 'cohorts' where animals simply died or were euthanized early - before nest data were collected
cage_exclude_cohorts = data.frame(cohort = c('', 'V-died', 'path-RML'))

mixed_cages = sqldf("
          select   b.cage, count(distinct b.cohort) n_cohort
          from     blind b
          where    b.cohort not in (select cohort from cage_exclude_cohorts) 
          group by 1
          having   count(distinct b.cohort) > 1
          order by 1
          ;")
# list of cages usable for nest data and which cohort they belong to
nest_cages = sqldf("
select   b.expt, b.cage, b.cohort
from     blind b
where    b.cage not in (select cage from mixed_cages)
and      b.cohort not in (select cohort from cage_exclude_cohorts) 
group by 1, 2, 3
order by 1, 2, 3
;")

n = sqldf("
select   nc.expt, nc.cohort, n.dpi, avg(n.comb) mean_nest, stdev(n.comb) sd_nest, count(*) n_cage
from     nests n, nest_cages nc
where    n.cage = nc.cage
group by 1, 2, 3
order by 1, 2, 3
;")

n = calculate_ci(n)
n$color = params$color[match(n$cohort, params$cohort)]
n$lty = params$lty[match(n$cohort, params$cohort)]










#### FIGURE S1: potency and time-to-effect
imgsave(paste('display_items/figure-s1.',imgmode,sep=''),width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,5,6,7,7),nrow=4,byrow=T)
layout(layout_matrix,heights=c(1.5,1,1,.2))

panel = 1

potency = read.table('data/ionis/potency_screening.tsv',sep='\t',header=T)
psmry = sqldf("select   region, treatment, panel, x, avg(mrna) mean_mrna, stdev(mrna) sd_mrna, count(*) n_mrna
              from     potency
              group by 1, 2, 3, 4
              order by 3, 4;
              ")
psmry = calculate_ci(psmry)
psmry$color = params$color[match(psmry$treatment,params$display)]
psmry$y = 7-psmry$x
potency$y = 7-potency$x
set.seed(1)
potency$yjit = jitter(potency$y, amount=.25)
potency$color = params$color[match(potency$treatment,params$display)]
xlims=c(0,1.5)
ylims=range(psmry$y) + c(-0.5,0.5)

par(mar=c(4,8,4,1))
for (pnl in unique(psmry$panel)) {
  rows = psmry$panel == pnl
  plot(NA,NA,xlim=xlims,ylim=ylims,xaxs='i',yaxs='i',ann=FALSE,axes=FALSE)
  axis(side=1, at=(0:4)/4, labels=percent((0:4)/4), lwd=1, lwd.ticks=1)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, at=psmry$y[rows], text=psmry$treatment[rows], col=psmry$color[rows], las=2, line=0.25)
  abline(v=1, lty=3)
  raw_rows = potency$panel == pnl
  points(x=potency$mrna[raw_rows], y=potency$yjit[raw_rows], pch=19, col=alpha(potency$color[raw_rows], ci_alpha))
  #rect(xleft=rep(0,sum(rows)), xright=psmry$mean_mrna[rows], ytop=psmry$y[rows]+0.33, ybottom=psmry$y[rows]-0.33, col=psmry$color, border=NA)
  segments(x0=psmry$mean_mrna[rows], y0=psmry$y[rows]+0.33, y1=psmry$y[rows]-0.33, lwd=1, col=psmry$color[rows])
  arrows(x0=psmry$l95[rows], x1=psmry$u95[rows], y0=psmry$y[rows], y1=psmry$y[rows], angle=90, length=0.035, code=3, lwd=1, col=psmry$color[rows])
  mtext(side=3, at=.5, line=0.5, text=psmry$region[rows][1])
  mtext(side=1, at=.5, line=2.5, text='mRNA (% saline)')
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.5, line = 0.5)
  panel = panel + 1
}

par(mar=c(4,5,3,1))
tte = read.table('data/ionis/time_to_effect.tsv',sep='\t',header=T)
tte$color = params$color[match(tte$treatment,params$display)]
tte$x = tte$days + seq(-.25, .25, by=1/6)

for (region in unique(tte$region)) {
  
  plot(NA, NA, xlim=c(0,8), ylim=c(0,1.2), axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, labels=NA, lwd=1, lwd.ticks=0)
  axis(side=1, at=c(1,4,7), lwd.ticks=1)
  mtext(side=1, line=2.5, text='days post-ICV')
  axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), las=2)
  mtext(side=2, line=3.0, text='mRNA (% saline)')
  abline(h=1, lwd=0.5, lty=3)
  
  # plot individual data points
  rows = tte$region == region
  points(tte$x[rows], tte$mrna[rows], pch=20, col=tte$color[rows])
  
  # plot means
  smry = sqldf(paste0("select treatment, days, avg(mrna) mean_mrna from tte where region = '",region,"' group by 1, 2 order by 1, 2;"))
  smry$color = params$color[match(smry$treatment,params$display)]
  segments(x0=smry$days-0.33, x1=smry$days+0.33, y0=smry$mean_mrna, lwd=2, col=smry$color)
  
  mtext(side=3, line=0.5, text=region)
  
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  panel = panel + 1
}

leg = sqldf("select treatment, color, avg(mrna) mean_mrna from tte group by 1, 2 order by 3 desc;")

par(mar=rep(.5,4))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), ann=F, axes=F, xaxs='i', yaxs='i')
legend('center', leg$treatment, col=leg$color, text.col=leg$color, pch=20, horiz=T, bty='n', cex=1.5)

dev.off() #### -- END FIGURE S1















#### FIGURE 1: replication of Raymond 2019 data with MOE ASOs

imgsave(paste('display_items/figure-1.',imgmode,sep=''), width=6.5*resx, height=4.8*resx, res=resx)

layout_matrix = matrix(c(1:6,7,7,7),nrow=3,byrow=T)
layout(layout_matrix,heights=c(1,1,.25))

par(mar=c(4,4.5,4,1))

panel = 1

m_full_xlims = c(-14, 350)
m_symp_xlims = c(90, 350)

sf = survfit(Surv(dpi, acm) ~ cohort, data=sqldf("select s.dpi, s.acm, b.cohort from surviv s, blind b where s.animal = b.animal and b.expt='M';"))
sf$cohort = gsub('cohort=','',names(sf$strata))
sf$treatment = params$display[match(sf$cohort, params$cohort)]
sf$color = params$color[match(sf$cohort, params$cohort)]

dosemarky = 1.08

plot(sf, ann=FALSE, axes=FALSE, xlim=m_full_xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$col, lwd=c(3,3,3,3,3))
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0:7)*50, lwd=0, lwd.ticks=1)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
abline(v=0)
par(xpd=TRUE)
text(x=mean(c(-14,76)),y=dosemarky,pos=3,labels='500µg doses')
points(x=c(-14,76),y=c(dosemarky,dosemarky),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.25, text='survival')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='mortality')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
panel = panel + 1

mw = subset(w, expt=='M')
mw$color = params$color[match(mw$cohort, params$cohort)]

plot(NA, NA, xlim=m_symp_xlims, ylim=c(15,35), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90, 330, by=30), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=0:7*5, labels=0:7*5, lwd=1, lwd.ticks=1, las=2)
abline(h=1,lty=3)
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=2, line=2.5, text='body weight (g)')
mtext(side=3, line=1.5, text='body weights')
for (coh in unique(mw$cohort)) {
  s = subset(mw, cohort==coh)
  points(s$dpi, s$mean_weight, col=s$color, pch=20, type='l', lwd=3)
  polygon(x=c(s$dpi[!is.na(s$sd_weight)],rev(s$dpi[!is.na(s$sd_weight)])), y=c(s$u95[!is.na(s$sd_weight)],rev(s$l95[!is.na(s$sd_weight)])), col=alpha(s$color[!is.na(s$sd_weight)],ci_alpha), border=NA)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
panel = panel + 1

mv = subset(v, expt=='M')
mv$color = params$color[match(mv$cohort, params$cohort)]

plot(NA, NA, xlim=m_symp_xlims, ylim=c(0,6), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90, 330, by=30), lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='days post-infection')
axis(side=2, at=0:8, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='symptoms per animal')
mtext(side=3, line=1.5, text='behaviorals')
for (coh in unique(mv$cohort)) {
  b = subset(mv, cohort==coh)
  points(b$dpi, b$mean_signs, type='l', lwd=3, col=b$color)
  polygon(x=c(b$dpi[b$n_animals > 1],rev(b$dpi[b$n_animals > 1])), y=c(b$u95[b$n_animals > 1],rev(b$l95[b$n_animals > 1])), col=alpha(b$color,ci_alpha), border=NA)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
panel = panel + 1

# the M cage data look nice but we were not scoring cages yet when we did expt P so the figure
# has an odd number of panels if we include nests just for M. comment this out for now:
# mn = subset(n, expt=='M' & n_cage > 1)
# par(mar=c(4,4,4,1))
# plot(NA, NA, xlim=m_symp_xlims, ylim=c(0,2.1), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
# axis(side=1, at=seq(90, 330, by=30), lwd=1, lwd.ticks=1, tck=-0.05)
# axis(side=2, at=c(0:4)/2, lwd=1, lwd.ticks=1, las=2)
# mtext(side=1, line=2.5, text='days post-infection')
# mtext(side=2, line=2.5, text='mean nest score')
# mtext(side=3, line=1.5, text='nest-building')
# for (coh in unique(mn$cohort)) {
#   s = subset(mn, cohort==coh)
#   points(s$dpi, s$mean_nest, col=s$color, pch=20, type='l', lwd=4)
# }
# mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
# panel = panel + 1

p_full_xlims = c(0, 300)
p_symp_xlims = c(120, 300)

sf = survfit(Surv(dpi, acm) ~ cohort, data=sqldf("select s.dpi, s.acm, b.cohort from surviv s, blind b where s.animal = b.animal and b.expt='P';"))
sf$cohort = gsub('cohort=','',names(sf$strata))
sf$treatment = params$display[match(sf$cohort, params$cohort)]
sf$color = params$color[match(sf$cohort, params$cohort)]

dosemarky = 1.08
plot(sf, ann=FALSE, axes=FALSE, xlim=m_full_xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$col, lwd=c(3,3,3,3,3))
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0:7)*50, lwd=0, lwd.ticks=1)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
abline(v=0)
par(xpd=TRUE)
text(x=mean(c(119,120)),y=dosemarky,pos=3,labels='500µg dose')
points(x=mean(c(119,120)),y=dosemarky,pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.25, text='survival')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='mortality')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
panel = panel + 1

pw = subset(w, expt=='P' & cohort != 'd120-MOE-sack')
# pw = subset(rw, expt=='P')
# the 'P' cohort was before we were systematic about weights before 120 dpi so they had a baseline at 88 but none again
# until 128 dpi. it is easier to visualize the data re-normalized to the 128 dpi timepoint instead of 88
# 
# prw1 = sqldf("
#              select   w.animal, w.dpi, w.weight / f.weight relwt
#              from     weights w, weights f
#              where    w.animal = f.animal
#              and      f.dpi = 128
#              order by 1, 2
#              ;")
# prw2 = sqldf("
#              select   b.expt, b.cohort, w.dpi, avg(w.relwt) mean_weight, stdev(w.relwt) sd_weight, count(*) n_weight
#              from     prw1 w, blind b
#              where    w.animal = b.animal
#              and      w.dpi > 120
#              and      b.expt = 'P'
#              group by 1, 2, 3
#              having   count(*) > 1
#              order by 1, 2, 3
#              ;")
# prw = calculate_ci(prw2)
# prw$color = params$color[match(prw$cohort, params$cohort)]
# pw = prw


plot(NA, NA, xlim=p_symp_xlims, ylim=c(15,35), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90, 330, by=30), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=0:7*5, labels=0:7*5, lwd=1, lwd.ticks=1, las=2)
abline(h=1,lty=3)
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=2, line=2.5, text='body weight (g)')
mtext(side=3, line=1.5, text='body weights')
for (coh in unique(pw$cohort)) {
  s = subset(pw, cohort==coh)
  points(s$dpi, s$mean_weight, col=s$color, pch=20, type='l', lwd=3)
  polygon(x=c(s$dpi[!is.na(s$sd_weight)],rev(s$dpi[!is.na(s$sd_weight)])), y=c(s$u95[!is.na(s$sd_weight)],rev(s$l95[!is.na(s$sd_weight)])), col=alpha(s$color[!is.na(s$sd_weight)],ci_alpha), border=NA)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
panel = panel + 1

pv = subset(v, expt=='P')
pv$color = params$color[match(pv$cohort, params$cohort)]

plot(NA, NA, xlim=p_symp_xlims, ylim=c(0,6), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90, 330, by=30), lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='days post-infection')
axis(side=2, at=0:8, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='symptoms per animal')
mtext(side=3, line=1.5, text='behaviorals')
for (coh in unique(pv$cohort[pv$cohort %in% params$cohort[!is.na(params$x)]])) {
  b = subset(pv, cohort==coh)
  points(b$dpi, b$mean_signs, type='l', lwd=3, col=b$color)
  polygon(x=c(b$dpi[b$n_animals > 1],rev(b$dpi[b$n_animals > 1])), y=c(b$u95[b$n_animals > 1],rev(b$l95[b$n_animals > 1])), col=alpha(b$color,ci_alpha), border=NA)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
panel = panel + 1

# blank plot for legend
par(mar=c(0,1,0,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', ann=F, axes=F)
par(xpd=F)
legend(x=0,y=1,legend=params$display[params$expt=='M'],horiz=T,col=params$color[params$expt=='M'],text.col=params$color[params$expt=='M'],title.col='#000000',lwd=3,bty='n',cex=1.5)

dev.off()











#### FIGURE 2: dose-response

imgsave(paste('display_items/figure-2.',imgmode,sep=''),width=6.5*resx,height=7*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,4,4,5,6,7),nrow=3,byrow=T)
layout(layout_matrix)

rm(s)
dr = sqldf("
           select   b.cohort, b.dose, avg(s.dpi) mean_dpi, stdev(s.dpi) sd_dpi, count(*) n
           from     surviv s, blind b
           where    s.animal = b.animal
           and      s.endpoint
           and      b.expt = 'D'
           group by 1
           order by 1
           ;")

inctime_table = sqldf("
                      select   d.cohort, d.dose, c.mean_dpi c_mean, c.sd_dpi c_sd, d.mean_dpi d_mean, d.sd_dpi d_sd, d.mean_dpi / c.mean_dpi dpi_ratio
                      from     dr c, dr d
                      where    c.dose = 0
                      order by 2
                      ;")

# we want 95%CI on the ratio of each dose to the 0 dose
# see https://stats.stackexchange.com/a/16354
# using mratios library -> ttestratio function
inctime_table$ttest_p = 1.0
inctime_table$l95 = 1.0
inctime_table$u95 = 1.0

for (i in 1:nrow(inctime_table)) {
  this_dose = inctime_table$dose[i]
  data = sqldf("select s.dpi, b.dose from surviv s, blind b where s.animal = b.animal and b.expt='D' and s.endpoint;")
  tr = ttestratio(data$dpi[data$dose==this_dose], data$dpi[data$dose==0], alternative='two.sided')
  inctime_table$ttest_p[i] = tr$p.value
  inctime_table$l95[i] = tr$conf.int[1]
  inctime_table$u95[i] = tr$conf.int[2]
}
inctime_table$color = params$color[match(inctime_table$cohort,params$cohort)]
inctime_table

dose_qpcr = read.table('data/ionis/doseresp_qpcr.tsv',sep='\t',header=T)
asmry = sqldf("
              select   dq.dose, avg(dq.mrna) mean_mrna, stdev(dq.mrna) sd_mrna, count(*) n_rep
              from     dose_qpcr dq
              where    treatment in ('active ASO 1','saline')
              and      region = 'cortex'
              group by 1
              order by 1
              ;")
asmry = calculate_ci(asmry)
asmry$cohort = paste('DR-',formatC(as.integer(asmry$dose),width=3,flag='0'),sep='')
asmry$color = params$color[match(asmry$cohort, params$cohort)]
asmry$color[asmry$dose==700] = '#000000'

rows500 = dose_qpcr$treatment=='active ASO 1' & dose_qpcr$region=='cortex' & dose_qpcr$dose==500
rows700 = dose_qpcr$treatment=='active ASO 1' & dose_qpcr$region=='cortex' & dose_qpcr$dose==700
t.test(dose_qpcr$mrna[rows500], dose_qpcr$mrna[rows700], alternative='two.sided')

inctime_table$qpcr_mean = asmry$mean_mrna[match(inctime_table$dose, asmry$dose)]
inctime_table$qpcr_l95 = asmry$l95[match(inctime_table$dose, asmry$dose)]
inctime_table$qpcr_u95 = asmry$u95[match(inctime_table$dose, asmry$dose)]

sf = survfit(Surv(dpi, acm) ~ dose, data=sqldf("select s.dpi, s.acm, b.cohort, b.dose from surviv s, blind b where s.animal = b.animal and b.expt='D';"))
sf$dose = gsub('dose=','',names(sf$strata))
sf$cohort = paste('DR-',formatC(as.integer(sf$dose),width=3,flag='0'),sep='')
sf$color = params$color[match(sf$cohort, params$cohort)]

# check log-rank test significance of 0 vs. 30 dose
just0_30 = sqldf("select s.dpi, s.acm, b.cohort, b.dose from surviv s, blind b where s.animal = b.animal and b.expt='D' and b.dose in (0,30);")
survdiff(Surv(dpi, acm) ~ dose, data=just0_30) 
survfit(Surv(dpi, acm) ~ dose, data=just0_30) 

# pseudopanel to hold master legend
par(mar=c(4,4,2,0))
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ann=FALSE,xaxs='i',yaxs='i')
par(xpd=T)
legend(x=0, y=1, legend=asmry$dose, lwd=3, col=asmry$color, bty='n', text.col=asmry$color, text.font=2, title='active ASO 1\ndose (µg)', title.col='black', cex=1.4)
par(xpd=F)

# A - dose vs. knockdown
xats = c(1:9, 10*(1:9), 100*(1:9), 1000)
xbigs = c(1, 10, 100, 1000)

par(mar=c(4,4,4,2))
plot(NA, NA, xlim=c(0,750),  ylim=c(0,1.1), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=seq(0,1000,by=100), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0,500), lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='dose (µg)')
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='PrP mRNA (% PBS)')
abline(h=c(.5,1),lwd=1,lty=3)
m = drm(mean_mrna ~ dose, data=asmry, fct=LL.4())
x = 0:2000
y = predict(m, newdata=data.frame(dose=x))
points(x, y, type='l', lwd=1)
par(xpd=T)
points(asmry$dose, asmry$mean_mrna, pch=19, col=asmry$color, cex=1.25)
arrows(x0=asmry$dose, y0=asmry$l95, y1=asmry$u95, lwd=2, col=asmry$color, angle=90, code=3, length=0.05)
par(xpd=F)
mtext(side=3, line=1.0, text='dose vs.\nknockdown')
mtext('A', side=3, cex=2, adj = -0.3, line = 1.5)

par(mar=c(4,4,4,2))
plot(NA, NA, xlim=c(0,1.2), ylim=c(.9,1.6), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
abline(v=1, lwd=1, lty=3)
abline(h=1, lwd=1, lty=3)
axis(side=1, at=(0:3)*.5, labels=percent((0:3)*.5), lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='PrP mRNA (% saline)')
axis(side=2, at=seq(.9,2,.1), labels=percent(seq(.9,2,.1)-1), lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='Δ time to endpoint')
points(inctime_table$qpcr_mean, inctime_table$dpi_ratio, pch=19, col=inctime_table$color)
segments(x0=inctime_table$qpcr_mean, x1=inctime_table$qpcr_mean, y0=inctime_table$l95, y1=inctime_table$u95, lwd=3, col=inctime_table$color)
segments(x0=inctime_table$qpcr_l95, x1=inctime_table$qpcr_u95, y0=inctime_table$dpi_ratio, y1=inctime_table$dpi_ratio, lwd=3, col=inctime_table$color)
par(xpd=T)
text(x=inctime_table$qpcr_mean*1.08, y=inctime_table$dpi_ratio*1.01, labels=paste(inctime_table$dose, 'µg'), col=inctime_table$col, font=2, adj=c(0,0), srt=20, cex=.9)
par(xpd=F)
mtext(side=3, line=1.0, text='knockdown vs.\nincubation time')
mtext('B', side=3, cex=2, adj = -0.3, line = 1.5)

dosemarky = 1.08
par(mar=c(4,6,4,4))
plot(sf, ann=FALSE, axes=FALSE, xlim=c(-14,300), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$col, lwd=c(3,3,3,3,3))
axis(side=1, at=c(0:7)*50, lwd=0, lwd.ticks=1)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
abline(h=0)
abline(v=0)
par(xpd=TRUE)
text(x=c(-14,76),y=c(dosemarky,dosemarky),pos=3,labels=c('dose 1','dose 2'))
points(x=c(-14,76),y=c(dosemarky,dosemarky),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.5, text='survival')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='mortality')
mtext('C', side=3, cex=2, adj = -0.1, line = 1.5)

#dw = subset(rw, expt=='D')
dw = subset(w, expt=='D')
dw$color = params$color[match(dw$cohort, params$cohort)]
par(mar=c(4,5,4,1))
plot(NA, NA, xlim=c(120,250), ylim=c(15,35), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90, 330, by=30), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=0:7*5, labels=0:7*5, lwd=1, lwd.ticks=1, las=2)
abline(h=1,lty=3)
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=2, line=3.0, text='body weight (g)')
mtext(side=3, line=1.5, text='body weights')
for (coh in unique(dw$cohort)) {
  s = subset(dw, cohort==coh)
  points(s$dpi, s$mean_weight, col=s$color, pch=20, type='l', lwd=3)
  polygon(x=c(s$dpi[!is.na(s$sd_weight)],rev(s$dpi[!is.na(s$sd_weight)])), y=c(s$u95[!is.na(s$sd_weight)],rev(s$l95[!is.na(s$sd_weight)])), col=alpha(s$color[!is.na(s$sd_weight)],ci_alpha), border=NA)
}
mtext('D', side=3, cex=2, adj = -0.1, line = 1.5)

dv = subset(v, expt=='D')
dv$color = params$color[match(dv$cohort, params$cohort)]

plot(NA, NA, xlim=c(120,250), ylim=c(0,6), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=seq(90, 330, by=30), lwd=1, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='days post-infection')
abline(v=120,lwd=2)
axis(side=2, at=0:8, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='symptoms per animal')
mtext(side=3, line=1.5, text='behaviorals')
for (coh in unique(dv$cohort)) {
  b = subset(dv, cohort==coh)
  points(b$dpi, b$mean_signs, type='l', lwd=3, col=b$color)
  polygon(x=c(b$dpi[b$n_animals > 1],rev(b$dpi[b$n_animals > 1])), y=c(b$u95[b$n_animals > 1],rev(b$l95[b$n_animals > 1])), col=alpha(b$color,ci_alpha), border=NA)
}
mtext('E', side=3, cex=2, adj = -0.1, line = 1.5)

dn = subset(n, expt=='D' & n_cage > 1)

par(mar=c(4,4,4,1))
plot(NA, NA, xlim=c(120,250), ylim=c(0,2.1), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=seq(90, 330, by=30), lwd=1, lwd.ticks=1, tck=-0.05)
axis(side=2, at=c(0:4)/2, lwd=1, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=2, line=2.5, text='mean nest score')
mtext(side=3, line=1.5, text='nest-building')
for (coh in unique(dn$cohort)) {
  s = subset(dn, cohort==coh)
  points(s$dpi, s$mean_nest, col=s$color, pch=20, type='l', lwd=4)
}
mtext('F', side=3, cex=2, adj = -0.1, line = 1.5)

dev.off()





























#### TABLE 1: effects of ASO treatment and genetic knockout across 5 strains

table1_left = sqldf("
select   p.x,
         b.strain, 
         avg(case when b.treatment = 'saline' then s.dpi else null end) saline_mean_dpi,
         stdev(case when b.treatment = 'saline' then s.dpi else null end) saline_sd_dpi,
         sum(case when b.treatment = 'saline' then 1 else 0 end) saline_n,
         avg(case when b.treatment = 'active ASO 1' then s.dpi else null end) active_mean_dpi,
         stdev(case when b.treatment = 'active ASO 1' then s.dpi else null end) active_sd_dpi,
         sum(case when b.treatment = 'active ASO 1' then 1 else 0 end) active_n,
         avg(case when b.treatment = 'active ASO 1' then s.dpi else null end) /  avg(case when b.treatment = 'saline' then s.dpi else null end) delta
from     surviv s, blind b, params p
where    s.animal = b.animal
and      b.cohort = p.cohort
and      b.expt in ('S','Q','R') and (not (b.expt = 'S' and b.strain='22L'))
and      s.endpoint
group by 1, 2
order by 1
;")
table1_left$saline_itime = paste(formatC(table1_left$saline_mean_dpi,format='f',digits=0),
                                 '±',formatC(table1_left$saline_sd_dpi,format='f',digits=0),sep='')
table1_left$active_itime = paste(formatC(table1_left$active_mean_dpi,format='f',digits=0),
                                 '±',formatC(table1_left$active_sd_dpi,format='f',digits=0),sep='')
table1_left$a_delta_pct = paste('+',percent(table1_left$delta - 1),sep='')

table1_right = sqldf("
select   p.x,
         b.strain, 
         avg(  case when b.cohort like '%WT%'  then s.dpi else null end) wt_mean_dpi,
         stdev(case when b.cohort like '%WT%'  then s.dpi else null end) wt_sd_dpi,
         sum(  case when b.cohort like '%WT%'  then 1 else 0 end) wt_n,
         avg(  case when b.cohort like '%ZH3%' then s.dpi else null end) zh3_mean_dpi,
         stdev(case when b.cohort like '%ZH3%' then s.dpi else null end) zh3_sd_dpi,
         sum(  case when b.cohort like '%ZH3%' then 1 else 0 end) zh3_n,
         avg(  case when b.cohort like '%ZH3%' then s.dpi else null end) /  avg(case when b.cohort like '%WT%'  then s.dpi else null end) delta
from     surviv s, blind b, params p
where    s.animal = b.animal
and      b.cohort = p.cohort
and      b.expt = 'G'
and      s.endpoint
group by 1, 2
order by 1
;")


table1_right$wt_itime = paste(formatC(table1_right$wt_mean_dpi,format='f',digits=0),
                                 '±',formatC(table1_right$wt_sd_dpi,format='f',digits=0),sep='')
table1_right$zh3_itime = paste(formatC(table1_right$zh3_mean_dpi,format='f',digits=0),
                                 '±',formatC(table1_right$zh3_sd_dpi,format='f',digits=0),sep='')
table1_right$g_delta_pct = paste('+',percent(table1_right$delta - 1),sep='')
table1_right = rbind(table1_right, rep(NA,ncol(table1_right))) # add a blank row for RML[ASO]

table1 = cbind(table1_left[,c('strain','saline_itime','saline_n','active_itime','active_n','a_delta_pct')],
               table1_right[,c('wt_itime','wt_n','zh3_itime','zh3_n','g_delta_pct')])

write.table(table1, 'display_items/table-1.tsv', sep='\t', col.names=T, row.names=F, quote=F)


#### TABLE 2

z_params = sqldf("
                 select   p.x, b.strain,
                 avg(surg.dpi) dpi_treatment,
                 avg(  case when b.treatment='saline'       and s.acm then s.dpi else null end) saline_mean_dpi,
                 stdev(case when b.treatment='saline'       and s.acm then s.dpi else null end) saline_sd_dpi,
                 sum(  case when b.treatment='saline'       and s.acm then 1 else 0 end) saline_n,
                 avg(  case when b.treatment='active ASO 6' and s.acm then s.dpi else null end) active_mean_dpi,
                 stdev(case when b.treatment='active ASO 6' and s.acm then s.dpi else null end) active_sd_dpi,
                 sum(  case when b.treatment='active ASO 6' and s.acm then 1 else 0 end) active_n
                 from     surviv s, blind b, surgery surg, params p
                 where    s.animal = b.animal
                 and      s.animal = surg.animal
                 and      surg.event = 'icv1'
                 and      b.expt = 'Z'
                 and      b.cohort = p.cohort
                 group by 1, 2
                 order by 1
                 ;")

z_params$timepoint = percent(z_params$dpi_treatment / z_params$saline_mean_dpi)
z_params$delta = paste('+',percent(z_params$active_mean_dpi / z_params$saline_mean_dpi - 1),sep='')

z_params$saline_itime = paste(formatC(z_params$saline_mean_dpi,format='f',digits=0),
                              '±',formatC(z_params$saline_sd_dpi,format='f',digits=0),sep='')
z_params$active_itime = paste(formatC(z_params$active_mean_dpi,format='f',digits=0),
                              '±',formatC(z_params$active_sd_dpi,format='f',digits=0),sep='')


z_params2 = sqldf("
                  select   p.x, b.strain,
                  sum(case when s.acm and s.dpi > z.saline_mean_dpi * 1.1 then 1 else 0 end) n_outlive,
                  sum(case when s.acm then 1 else 0 end) n_total,
                  avg(case when s.acm and s.dpi > z.saline_mean_dpi * 1.1 then s.dpi else null end) mean_dpi_outlive,
                  stdev(case when s.acm and s.dpi > z.saline_mean_dpi * 1.1 then s.dpi else null end) sd_dpi_outlive,
                  avg(case when s.acm and s.dpi > z.saline_mean_dpi * 1.1 then s.dpi else null end)/z.saline_mean_dpi-1 mean_delta_outlive
                  from     surviv s, blind b, params p, z_params z
                  where    s.animal = b.animal
                  and      b.expt = 'Z'
                  and      b.cohort = p.cohort and b.strain = z.strain
                  and      b.treatment = 'active ASO 6'
                  group by 1, 2
                  order by 1
                  ;")
z_params2$proportion_outlive = paste(z_params2$n_outlive,z_params2$n_total,sep='/')
z_params2$dpi_outlive = paste(round(z_params2$mean_dpi_outlive,digits=0),'±',round(z_params2$sd_dpi_outlive,digits=0),sep='')
z_params2$delta_outlive = paste0('+',percent(z_params2$mean_delta_outlive))

table2 = z_params[,c('strain','dpi_treatment','timepoint','saline_itime','saline_n','active_itime','active_n','delta')]
table2$proportion_outlive = z_params2$proportion_outlive[match(table2$strain, z_params2$strain)]
table2$dpi_outlive = z_params2$dpi_outlive[match(table2$strain, z_params2$strain)]
table2$delta_outlive = z_params2$delta_outlive[match(table2$strain, z_params2$strain)]

write.table(table2, 'display_items/table-2.tsv', sep='\t', col.names=T, row.names=F, quote=F)

# do the proportions differ?
fisher.test(z_params2[,c('n_outlive','n_total')], alternative='two.sided')










### FIGURE S2: survival curves from table 1 ASO treatments

imgsave(paste('display_items/figure-s2.',imgmode,sep=''),width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,5,6,7,8), nrow=4, byrow=T)
layout(layout_matrix)

ss = sqldf("
select   b.expt, b.strain, b.cohort, b.treatment, s.dpi, s.acm, p.color, p.lty
from     surviv s, blind b, params p
where    b.animal = s.animal
and      b.cohort = p.cohort
and      b.expt in ('S','Q','R')
order by 1, 2, 3, 4, 5
;")

ss$title = ss$strain
ss$title[ss$strain=='22L' & ss$expt=='S'] = '22L (original)'
ss$title[ss$strain=='22L' & ss$expt=='Q'] = '22L (repeat)'

par(mar=c(4,5,4,1))

s_full_xlims = c(-25,325)
dosemarky = 1.08
panel = 1
for (title in c('RML','22L (original)','22L (repeat)','Fukuoka-1','ME7','OSU','RML[ASO]')) {
  subs = ss[ss$title==title,]
  sf = survfit(Surv(dpi, acm) ~ treatment, data=subs)
  sf$cohort = gsub('treatment=','',names(sf$strata))
  sf$color = params$color[match(sf$cohort, params$display)]
  plot(sf, ann=FALSE, axes=FALSE, xlim=s_full_xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=3)
  axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=c(0:7)*50, lwd=0, lwd.ticks=1)
  axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
  abline(v=0)
  par(xpd=TRUE)
  text(x=mean(c(-14,76)),y=dosemarky,pos=3,labels='500µg doses')
  points(x=c(-14,76),y=c(dosemarky,dosemarky),pch=25,col='black',bg='black')
  par(xpd=FALSE)
  mtext(side=2, line=3.25, text='survival')
  mtext(side=1, line=2.5, text='days post-infection')
  mtext(side=3, line=1.5, text=subs$title[1])
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
  panel = panel + 1
}

plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', ann=F, axes=F)
legend('center',sf$cohort,col=sf$color,lwd=3,text.col=sf$color,bty='n',title.col='black',title='treatment')

dev.off()




#### FIGURE S3: do weights, symptoms, and nests support the differences between strains in G cohort being real?

imgsave(paste('display_items/figure-s3.',imgmode,sep=''),width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,1,2,2,3,4,4,5,5,6), nrow=2, byrow=T)
layout(layout_matrix)

panel = 1

rm(s)
gs = sqldf("
select   b.expt, b.cohort, s.dpi, s.acm, p.color2, p.lty2
from     surviv s, blind b, params p
where    b.animal = s.animal
and      b.cohort = p.cohort
and      b.expt = 'G'
order by 1, 2, 3
;")

gsurv = survfit(Surv(dpi, acm) ~ cohort, data=gs)
gsurv$cohort = gsub('cohort=','',names(gsurv$strata))
gsurv$color = params$color2[match(gsurv$cohort, params$cohort)]
gsurv$lty = params$lty2[match(gsurv$cohort, params$cohort)]

par(mar=c(4,5,3,1))
plot(gsurv, ann=FALSE, axes=FALSE, xlim=c(0,450), ylim=c(0,1.05), xaxs='i', yaxs='i', lwd=strain_lwd, col=gsurv$color, lty=gsurv$lty)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=(0:9)*50, labels=NA, lwd=1, lwd.ticks=1, tck=-0.05)
axis(side=1, at=c(0:5)*100, lwd=0, lwd.ticks=0, line=-0.5)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='proportion surviving')
mtext(side=3, line=0, text='mortality')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


grw = subset(rw, expt == 'G')

par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(90,450), ylim=c(0.5,1.25), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90,450,by=30), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(90,450,by=30), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=(2:6)/4, labels=percent((2:6)/4), las=2)
mtext(side=2, line=3, text='relative weight')
abline(h=1, lty=3)
for (coh in unique(grw$cohort)) {
  s = subset(grw, cohort==coh)
  points(s$dpi, s$mean_weight, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='body weight')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# blank panel for legend
par(mar=c(4,1,3,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
par(xpd=T)
leg1 = subset(params, expt=='G' & grepl('22L',cohort))
legend(x=0,y=1,legend=leg1$display2,col='black',lty=leg1$lty2,text.col='black',lwd=strain_lwd,text.font=2,bty='n',cex=1)
leg2 = subset(params, expt=='G' & grepl('WT',cohort))
leg2 = leg2[with(leg2, order(x)),]
legend(x=0,y=0.7,legend=leg2$display,col=leg2$color2,lty=leg2$lty2,text.col=leg2$color2,lwd=strain_lwd,text.font=2,bty='n',cex=1)
par(xpd=F)

gv = subset(v, expt == 'G')

par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(90,450), ylim=c(0,6), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90,450,by=30), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(90,450,by=30), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=0:6, labels=0:6, las=2)
mtext(side=2, line=3, text='symptoms per animal')
for (coh in unique(gv$cohort)) {
  s = subset(gv, cohort==coh)
  points(s$dpi, s$mean_signs, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='symptoms')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


gn = subset(n, expt=='G')
gn$color = params$color2[match(gn$cohort, params$cohort)]
gn$lty = params$lty2[match(gn$cohort, params$cohort)]


par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(90,450), ylim=c(0,2.1), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90,450,by=30), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(90,450,by=30), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=0:2, labels=0:2, las=2)
mtext(side=2, line=3, text='mean nest score')
for (coh in unique(gn$cohort)) {
  s = subset(gn, cohort==coh)
  points(s$dpi, s$mean_nest, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='nest-building')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# blank panel to match legend above
par(mar=c(4,1,3,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')

dev.off()




















#### FIGURE S4: Z cohort

imgsave(paste('display_items/figure-s4.',imgmode,sep=''),width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,1,2,2,3,4,4,5,5,6), nrow=2, byrow=T)
layout(layout_matrix)

panel = 1

rm(s)
zs = sqldf("
           select   b.expt, b.cohort, s.dpi, s.acm, p.color2, p.lty2
           from     surviv s, blind b, params p
           where    b.animal = s.animal
           and      b.cohort = p.cohort
           and      b.expt = 'Z'
           order by 1, 2, 3
           ;")

zsurv = survfit(Surv(dpi, acm) ~ cohort, data=zs)
zsurv$cohort = gsub('cohort=','',names(zsurv$strata))
zsurv$color = params$color2[match(zsurv$cohort, params$cohort)]
zsurv$lty = params$lty2[match(zsurv$cohort, params$cohort)]

par(mar=c(4,5,3,1))
plot(zsurv, ann=FALSE, axes=FALSE, xlim=c(0,300), ylim=c(0,1.05), xaxs='i', yaxs='i', lwd=strain_lwd, col=zsurv$color, lty=zsurv$lty)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(0,450,by=50), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(0,500,by=100), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='proportion surviving')
mtext(side=3, line=0, text='mortality')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


zw = subset(rw, expt == 'Z')
zw$lty = params$lty2[match(zw$cohort,params$cohort)]
zw$color = params$color2[match(zw$cohort,params$cohort)]


par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(90,270), ylim=c(.5,1.25), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90,450,by=30), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(90,450,by=30), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=(2:6)/4, labels=percent((2:6)/4), las=2)
mtext(side=2, line=3, text='relative weight')
abline(h=1, lty=3)
for (coh in unique(zw$cohort)) {
  s = subset(zw, cohort==coh)
  points(s$dpi, s$mean_weight, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='body weight')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# blank panel for legend
par(mar=c(4,1,3,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')
par(xpd=T)
leg1 = subset(params, expt=='Z' & grepl('22L',cohort))
legend(x=0,y=1,legend=leg1$display2,col='black',lty=leg1$lty2,text.col='black',lwd=strain_lwd,text.font=2,bty='n',cex=1)
leg2 = subset(params, expt=='Z' & grepl('ACTIVE6',cohort))
leg2 = leg2[with(leg2, order(x)),]
legend(x=0,y=0.7,legend=leg2$display,col=leg2$color2,lty=leg2$lty2,text.col=leg2$color2,lwd=strain_lwd,text.font=2,bty='n',cex=1)
par(xpd=F)

zv = subset(v, expt == 'Z')
zv$lty = params$lty2[match(zv$cohort,params$cohort)]
zv$color = params$color2[match(zv$cohort,params$cohort)]

par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(90,270), ylim=c(0,6), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90,450,by=30), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(90,450,by=30), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=0:6, labels=0:6, las=2)
mtext(side=2, line=3, text='symptoms per animal')
for (coh in unique(zv$cohort)) {
  s = subset(zv, cohort==coh)
  points(s$dpi, s$mean_signs, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='symptoms')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


zn = subset(n, expt=='Z')
zn$color = params$color2[match(zn$cohort, params$cohort)]
zn$lty = params$lty2[match(zn$cohort, params$cohort)]


par(mar=c(4,5,3,1))
plot(NA, NA, xlim=c(90,270), ylim=c(0,2.1), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(90,450,by=30), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(90,450,by=30), line=-0.5, lwd=0, lwd.ticks=0)
mtext(side=1, line=2, text='days post-infection')
axis(side=2, at=0:2, labels=0:2, las=2)
mtext(side=2, line=3, text='mean nest score')
for (coh in unique(zn$cohort)) {
  s = subset(zn, cohort==coh)
  points(s$dpi, s$mean_nest, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='nest-building')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# blank panel to match legend above
par(mar=c(4,1,3,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F, xaxs='i', yaxs='i')

dev.off()
























#### FIGURE 4: MRI bioluminescence study + control ASO NfL study

imgsave(paste('display_items/figure-4.',imgmode,sep=''), width=6.5*resx, height=3.25*resx, res=resx)

#layout_matrix = matrix(c(1,4,2,5,3,6),nrow=3,byrow=T)
layout_matrix = matrix(c(1,3,2,4),nrow=2,byrow=T)
layout(layout_matrix,heights=c(1,1))

panel = 1

##### 'C' cohort

expt = 'C'
master = read.table(paste('data/mri/',expt,'_survival.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
cohorts = read.table(paste('data/mri/',expt,'_cohorts.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
cnfl = read.table(paste('data/mri/',expt,'_nfl.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')


cnfl_smry = sqldf("
                  select   m.treatment cohort, n.dpi, avg(n.nfl_pgml) mean_nfl, stdev(n.nfl_pgml) sd_nfl, count(*) n_nfl
                  from     cnfl n, master m
                  where    n.animal = m.animal
                  group by 1, 2 
                  order by 1, 2
                  ;")
cnfl_smry$color = cohorts$color[match(cnfl_smry$cohort, cohorts$cohort)]
cnfl_smry = calculate_ci(cnfl_smry)
aso_dpi = 120
xlims = c(0,280)
ylims = c(0,2000)
marky = 1750

par(mar=c(1,4,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(min(xlims),max(xlims),50), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(min(xlims),max(xlims),50), lwd=0, line=-0.75)
axis(side=2, at=(0:(max(ylims)/100))*100, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=(0:ceiling(max(ylims/500)))*500, labels=NA, lwd=0, lwd.ticks=1, tck=-0.050)
axis(side=2, at=(0:ceiling(max(ylims/500)))*500, labels=formatC((0:ceiling(max(ylims/500)))*500, big.mark=','), lwd=0, lwd.ticks=0, las=2, line=-0.5)
mtext(side=2, line=3.0, text='plasma NfL (pg/mL)')
mtext(side=3, line=1.5, text='neurofilament')
# axis(side=2, at=(0:(max(ylims)/100))*100, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
# axis(side=2, at=(0:ceiling(max(ylims/500)))*500, labels=formatC((0:ceiling(max(ylims/500)))/2, format='f', digits=1), lwd=0, lwd.ticks=1, tck=-0.050, las=2)
# mtext(side=2, line=3.5, text='plasma NfL (ng/mL)')
par(xpd=TRUE)
text(x=aso_dpi,y=marky,pos=3,labels=c('300 µg dose'),cex=0.8)
points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
segments(x0=aso_dpi,x1=aso_dpi,y0=0,y1=marky,lwd=1)
par(xpd=FALSE)

# for (ms in unique(cnfl$animal)) {
#   subs = subset(cnfl, animal==ms)
#   points(subs$timepoint, subs$nfl_pgml, type='l', col=subs$color)
# }
for (coh in unique(cnfl_smry$cohort)) { # do them in reverse order so pre-treatment is on top at 81 dpi
  subs = subset(cnfl_smry, cohort==coh)
  color = subs$color[1]
  points(subs$dpi, subs$mean_nfl, type='l', lwd=3, col=color)
  points(subs$dpi, subs$mean_nfl, pch=20, col=color)
  polysubs = subset(subs, !is.na(sd_nfl)) 
  polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
}

par(xpd=T)
xleg = 175
yleg = 2200
legend(x=xleg,y=yleg,legend=cohorts$cohort,col=cohorts$color,text.col=cohorts$color,text.font=2,lwd=3,lty=cohorts$lty,bty='n',title='cohorts',title.col='black',cex=0.65)
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# 
# cw = read.table('data/mri/C_weights.tsv',sep='\t',header=T)
# cw$cohort = master$treatment[match(cw$animal, master$animal)]
# crw = sqldf("
#             select   b.cohort, b.animal, b.dpi, b.weight, b.weight/baseline.weight relwt
#             from     cw b, cw baseline
#             where    b.animal = baseline.animal
#             and      baseline.dpi in (51,53)
#             order by 1, 2, 3
#             ;")
# cw_smry = sqldf("
#                 select   cohort, dpi, avg(weight) mean_wt, stdev(weight) sd_wt, count(*) n_wt
#                 from     cw
#                 group by 1, 2
#                 having   count(*) > 1
#                 order by 1, 2
#                 ;")
# crw_smry = sqldf("
#                  select   cohort, dpi, avg(relwt) mean_wt, stdev(relwt) sd_wt, count(*) n_wt
#                  from     crw
#                  group by 1, 2
#                  having   count(*) > 1
#                  ;")
# crw_smry$color = cohorts$color[match(crw_smry$cohort, cohorts$cohort)]
# crw_smry$lty = cohorts$lty[match(crw_smry$cohort, cohorts$cohort)]
# crw_smry = calculate_ci(crw_smry)
# cw_smry = calculate_ci(cw_smry)
# 
# par(mar=c(0,5,3,2))
# plot(NA, NA, xlim=xlims, ylim=c(15,35), ann=F, axes=F, xaxs='i', yaxs='i')
# axis(side=1, at=seq(min(xlims),max(xlims),10), labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
# axis(side=1, at=seq(min(xlims),max(xlims),50), labels=NA, tck=-0.05)
# axis(side=2, at=0:7*5, las=2)
# mtext(side=2, line=3.5, text='weight (g)')
# abline(h=1, lty=3)
# cohorts_for_loop = unique(cw_smry$cohort)
# for (coh in cohorts_for_loop) {
#   s = subset(cw_smry, cohort == coh)
#   color = cohorts$color[cohorts$cohort==coh]
#   lty = cohorts$lty[cohorts$cohort==coh]
#   points(s$dpi, s$mean_wt, type='l', lwd=default_lwd, col=color, lty=lty)
#   polysubs = subset(s, !is.na(sd_wt))
#   polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
# }
# marky=30
# par(xpd=T)
# text(x=aso_dpi,y=marky,pos=3,labels=c('300 µg dose'))
# points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
# segments(x0=aso_dpi,x1=aso_dpi,y0=0,y1=marky,lwd=1)
# par(xpd=F)
# mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
# panel = panel + 1

sf = survfit(Surv(dpi) ~ treatment, data=master)
sf$treatment = gsub('treatment=','',names(sf$strata))
sf$color = cohorts$color[match(sf$treatment, cohorts$cohort)]

marky = 1.05
par(mar=c(2,4,3,1))
plot(sf, ann=FALSE, axes=FALSE, xlim=xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(min(xlims),max(xlims),50), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(min(xlims),max(xlims),50), lwd=0, line=-.75)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=3, line=1.0, text='mortality')
par(xpd=TRUE)
text(x=aso_dpi,y=marky,pos=3,labels=c('300 µg dose'),cex=0.8)
points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.0, text='survival')
mtext(side=1, line=1.0, text='days post-infection')

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

aso_dpi = 84
marky = 8

expt = 'B'
master = read.table(paste('data/mri/',expt,'_survival.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
bli = read.table(paste('data/mri/',expt,'_bli.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
cohorts = read.table(paste('data/mri/',expt,'_cohorts.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

bli$cohort = master$treatment[match(bli$animal, master$animal)]
bli$cohort[bli$cohort!='uninfected' & bli$img_dpi < aso_dpi] = 'pre-treatment'

bsmry = sqldf("
              select   img_dpi, cohort, avg(flux) mean_bli, stdev(flux) sd_bli, sum(case when flux is not null then 1 else 0 end) n_bli
              from     bli b
              group by 1, 2
              order by 1, 2
              ;")
bsmry = calculate_ci(bsmry)

par(mar=c(1,4,3,1))
plot(NA,NA,xlim=c(0,260),ylim=c(0,12.5),xaxs='i',yaxs='i',ann=FALSE,axes=FALSE)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(min(xlims),max(xlims),50), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(min(xlims),max(xlims),50), lwd=0, line=-0.75)
axis(side=2, at=(0:3)*5, las=2, tck=-0.05)
axis(side=2, at=0:15, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
mtext(side=3, line=1.5, text='bioluminescence')
par(xpd=TRUE)
text(x=aso_dpi,y=marky,pos=3,labels=c('500 µg dose'),cex=0.8)
points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
segments(x0=aso_dpi,x1=aso_dpi,y0=0,y1=marky,lwd=1)
par(xpd=FALSE)
mtext(side=2, line=2.5, text='flux (10^7 photons)')
cohorts_for_loop = c("active ASO 1", "saline", "uninfected", "pre-treatment")
for (coh in cohorts_for_loop) { # do them in reverse order so pre-treatment is on top at 81 dpi
  if (coh %in% c('uninfected','pre-treatment')) {
    subs = subset(bsmry, cohort == coh)
  } else {
    subs = subset(bsmry, cohort == coh | cohort == 'pre-treatment' & img_dpi == 81)
  }
  color = cohorts$color[cohorts$cohort==coh]
  lty = cohorts$lty[cohorts$cohort==coh]
  points(subs$img_dpi, subs$mean, type='l', lwd=3, lty=lty, col=color)
  points(subs$img_dpi, subs$mean, pch=20, col=color)
  polysubs = subset(subs, !is.na(sd_bli))
  polygon(x=c(polysubs$img_dpi,rev(polysubs$img_dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
}

par(xpd=T)
xleg = 175
yleg = 17.5
legend(x=xleg,y=yleg,legend=cohorts$cohort,col=cohorts$color,text.col=cohorts$color,text.font=2,lwd=3,lty=cohorts$lty,bty='n',title='cohorts',title.col='black',cex=0.65)
rml_bounds = c(6.75,14)
arrlen = 4
segments(x0=xleg, y0=min(rml_bounds), y1=max(rml_bounds), lwd=1)
segments(x0=xleg, x1=xleg+arrlen, y0=min(rml_bounds), lwd=1)
segments(x0=xleg, x1=xleg+arrlen, y0=max(rml_bounds), lwd=1)
text(x=xleg-5, y=mean(rml_bounds), srt=90, labels='prion-infected',cex=0.65)
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# 
# 
# # weights
# bw = read.table('data/mri/B_weights.tsv',sep='\t',header=T)
# bw$cohort = master$treatment[match(bw$animal, master$animal)]
# bw$cohort[bw$cohort!='uninfected' & bw$dpi < aso_dpi] = 'pre-treatment'
# # check that all had baseline at same dpi:
# # sqldf("select animal, min(dpi) from bw where weight is not null group by 1 order by 2 desc;")
# brw = sqldf("
#             select   b.cohort, b.animal, b.dpi, b.weight, b.weight/baseline.weight relwt
#             from     bw b, bw baseline
#             where    b.animal = baseline.animal
#             and      baseline.dpi = 11
#             order by 1, 2, 3
#             ;")
# bw_smry = sqldf("
#                 select   cohort, dpi, avg(weight) mean_wt, stdev(weight) sd_wt, count(*) n_wt
#                 from     bw
#                 group by 1, 2
#                 having   count(*) > 1
#                 order by 1, 2
#                 ;")
# bw_smry = calculate_ci(bw_smry)
# brw_smry = sqldf("
#                  select   cohort, dpi, avg(relwt) mean_wt, stdev(relwt) sd_wt, count(*) n_wt
#                  from     brw
#                  group by 1, 2
#                  having   count(*) > 1
#                  union all -- allow pre-treatment plot red line to extend up to where indiv cohort lines start
#                  select   'pre-treatment' as cohort, dpi, avg(relwt) mean_wt, stdev(relwt) sd_wt, count(*) n_wt
#                  from     brw
#                  where    cohort in ('saline','active ASO 1')
#                  and      dpi = 89
#                  group by 1, 2
#                  having count(*) > 1
#                  order by 1, 2
#                  ;")
# brw_smry$color = cohorts$color[match(brw_smry$cohort, cohorts$cohort)]
# brw_smry$lty = cohorts$lty[match(brw_smry$cohort, cohorts$cohort)]
# brw_smry = calculate_ci(brw_smry)
# # 
# # plot(NA, NA, xlim=c(0,260), ylim=c(0.75,1.5), ann=F, axes=F, xaxs='i', yaxs='i')
# # axis(side=1, at=c(0:22)*25, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
# # axis(side=1, at=c(0:10)*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
# # axis(side=2, at=(2:6)/4, labels=percent((2:6)/4), las=2)
# # mtext(side=2, line=3.5, text='relative weight')
# # abline(h=1, lty=3)
# # cohorts_for_loop = c("active ASO 1", "saline", "uninfected", "pre-treatment")
# # for (coh in cohorts_for_loop) {
# #   s = subset(brw_smry, cohort == coh)
# #   color = cohorts$color[cohorts$cohort==coh]
# #   lty = cohorts$lty[cohorts$cohort==coh]
# #   points(s$dpi, s$mean_wt, type='l', lwd=default_lwd, col=color, lty=lty)
# #   polysubs = subset(s, !is.na(sd_wt))
# #   polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
# # }
# # marky=1.3
# # text(x=aso_dpi,y=marky,pos=3,labels=c('500 µg dose'))
# # points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
# # segments(x0=aso_dpi,x1=aso_dpi,y0=0,y1=marky,lwd=1)
# # mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
# # panel = panel + 1
# 
# par(mar=c(0,5,3,2))
# plot(NA, NA, xlim=c(0,260), ylim=c(15,35), ann=F, axes=F, xaxs='i', yaxs='i')
# axis(side=1, at=c(0:30)*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
# axis(side=1, at=c(0:10)*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
# axis(side=2, at=0:7*5, las=2)
# mtext(side=2, line=3.5, text='weight (g)')
# abline(h=1, lty=3)
# cohorts_for_loop = c("active ASO 1", "saline", "uninfected", "pre-treatment")
# for (coh in cohorts_for_loop) {
#   s = subset(bw_smry, cohort == coh)
#   color = cohorts$color[cohorts$cohort==coh]
#   lty = cohorts$lty[cohorts$cohort==coh]
#   points(s$dpi, s$mean_wt, type='l', lwd=default_lwd, col=color, lty=lty)
#   polysubs = subset(s, !is.na(sd_wt))
#   polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
# }
# marky=35
# par(xpd=T)
# text(x=aso_dpi,y=marky,pos=3,labels=c('500 µg dose'))
# points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
# segments(x0=aso_dpi,x1=aso_dpi,y0=0,y1=marky,lwd=1)
# par(xpd=F)
# mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
# panel = panel + 1

sf = survfit(Surv(dpi, acm) ~ treatment, data=subset(master, treatment != 'uninfected'))
sf$treatment = gsub('treatment=','',names(sf$strata))
sf$color = cohorts$color[match(sf$treatment, cohorts$cohort)]

marky = 1.05
par(mar=c(2,4,3,1))
plot(sf, ann=FALSE, axes=FALSE, xlim=c(0,260), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
points(x=c(0,aso_dpi), y=c(1,1), type='l', lwd=3, col=cohorts$color[cohorts$cohort=='pre-treatment'])
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=seq(min(xlims),max(xlims),50), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=seq(min(xlims),max(xlims),50), lwd=0, line=-0.75)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=3, line=1.0, text='mortality')
par(xpd=TRUE)
text(x=aso_dpi,y=marky,pos=3,labels=c('500 µg dose'),cex=0.8)
points(x=aso_dpi,y=marky,pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.0, text='survival')
mtext(side=1, line=1.0, text='days post-infection')

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

dev.off() #### -- end Figure 4





























#### TABLE "T" (just for presentations, in-line text)

t_params_ep = sqldf("
                 select   p.x, b.timepoint,
                          avg(surg.dpi) actual_dpi_treatment,
                          avg(  case when b.treatment='saline'       and s.endpoint then s.dpi else null end) saline_mean_dpi,
                          stdev(case when b.treatment='saline'       and s.endpoint then s.dpi else null end) saline_sd_dpi,
                          sum(  case when b.treatment='saline'       and s.endpoint then 1 else 0 end) saline_n,
                          avg(  case when b.treatment='active ASO 6' and s.endpoint then s.dpi else null end) active_mean_dpi,
                          stdev(case when b.treatment='active ASO 6' and s.endpoint then s.dpi else null end) active_sd_dpi,
                          sum(  case when b.treatment='active ASO 6' and s.endpoint then 1 else 0 end) active_n
                 from     surviv s, blind b, surgery surg, params p
                 where    s.animal = b.animal
                 and      s.animal = surg.animal
                 and      surg.event = 'icv1'
                 and      b.expt = 'T'
                 and      b.cohort = p.cohort
                 and      b.timepoint is not null
                 group by 1, 2
                 order by 1
                 ;")

# inspect the data
# View(sqldf("select * from blind b, surviv s where b.animal = s.animal and b.timepoint in (105,120) and b.expt ='T' order by b.animal;"))
# View(sqldf("select * from blind b, surviv s where b.animal = s.animal and b.timepoint in (-7,1,28,56,78) and b.expt ='T' order by b.animal;"))

t_params_acm = sqldf("
                    select   p.x, b.timepoint,
                    avg(surg.dpi) actual_dpi_treatment,
                    avg(  case when b.treatment='saline'       and s.acm then s.dpi else null end) saline_mean_dpi,
                    stdev(case when b.treatment='saline'       and s.acm then s.dpi else null end) saline_sd_dpi,
                    sum(  case when b.treatment='saline'       and s.acm then 1 else 0 end) saline_n,
                    avg(  case when b.treatment='active ASO 6' and s.acm then s.dpi else null end) active_mean_dpi,
                    stdev(case when b.treatment='active ASO 6' and s.acm then s.dpi else null end) active_sd_dpi,
                    sum(  case when b.treatment='active ASO 6' and s.acm then 1 else 0 end) active_n
                    from     surviv s, blind b, surgery surg, params p
                    where    s.animal = b.animal
                    and      s.animal = surg.animal
                    and      surg.event = 'icv1'
                    and      b.expt = 'T'
                    and      b.cohort = p.cohort
                    and      b.timepoint is not null
                    group by 1, 2
                    order by 1
                    ;")

t_params_ep$delta = paste('+',percent(t_params_ep$active_mean_dpi / t_params_ep$saline_mean_dpi - 1),sep='')
t_params_ep$saline_itime = paste(formatC(t_params_ep$saline_mean_dpi,format='f',digits=0),
                                 '±',formatC(t_params_ep$saline_sd_dpi,format='f',digits=0),sep='')
t_params_ep$active_itime = paste(formatC(t_params_ep$active_mean_dpi,format='f',digits=0),
                                 '±',formatC(t_params_ep$active_sd_dpi,format='f',digits=0),sep='')

t_params_acm$delta = paste('+',percent(t_params_acm$active_mean_dpi / t_params_acm$saline_mean_dpi - 1),sep='')
t_params_acm$saline_itime = paste(formatC(t_params_acm$saline_mean_dpi,format='f',digits=0),
                                 '±',formatC(t_params_acm$saline_sd_dpi,format='f',digits=0),sep='')
t_params_acm$active_itime = paste(formatC(t_params_acm$active_mean_dpi,format='f',digits=0),
                                 '±',formatC(t_params_acm$active_sd_dpi,format='f',digits=0),sep='')

table_t_ep = t_params_ep[,c('timepoint','saline_itime','saline_n','active_itime','active_n','delta')]
table_t_acm = t_params_acm[,c('timepoint','saline_itime','saline_n','active_itime','active_n','delta')]

write.table(table_s1_acm, 'display_items/table-s3.tsv', sep='\t', col.names=T, row.names=F, quote=F)


















#### FIGURE 5: timepoint dependence
imgsave(paste('display_items/figure-5.',imgmode,sep=''),width=6.5*resx,height=7.5*resx,res=resx)

#layout_matrix = matrix(c(1,1,1,2,2,2,3,4,5,6,7,8), nrow=4, byrow=T) # from when this and ST were all one figure
#layout(layout_matrix)
#par(mfrow=c(2,1))
# layout_matrix = matrix(c(1:7,rep(8,7)),nrow=2,byrow=T)
# layout(layout_matrix, widths=c(1.25,rep(1,6)))

layout_matrix = matrix(1:8,nrow=8,byrow=T)
layout(layout_matrix, heights=c(1.5,rep(1,5),1.5,4))

panel = 1

t_surv = sqldf("
               select   b.timepoint, b.treatment, s.animal, s.dpi, s.acm, b.cohort, p.color, p.lty
               from     surviv s, blind b, params p
               where    s.animal = b.animal
               and      b.cohort = p.cohort
               and      b.expt = 'T' and b.timepoint is not null
               order by 1, 2, 4
               ;")

infect_col = '#FF2017'
treat_col = '#000000'

xlims = c(-25,500)


for (tpt in unique(t_surv$timepoint)) {
  
  if (panel == 1) {
    par(mar=c(0.5,6,3,2))
  } else if (panel == 7) {
    par(mar=c(3,6,1.25,2))
  } else {
    par(mar=c(0.5,6,1.25,2))
  }
  
  tsurv = survfit(Surv(dpi, acm) ~ cohort, data=subset(t_surv, timepoint==tpt))
  tsurv$cohort = gsub('cohort=','',names(tsurv$strata))
  tsurv$color = params$color[match(tsurv$cohort, params$cohort)]
  
  anmls = blind$animal[blind$cohort %in% tsurv$cohort]
  surg_dpi = sort(unique(surgery$dpi[grepl('icv',surgery$event) & surgery$animal %in% anmls]))
  
  plot(tsurv, ann=FALSE, axes=FALSE, xlim=xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', lwd=default_lwd, col=tsurv$color, lty=1)

  axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  if (panel == 7) {
    axis(side=1, at=0:10*50, lwd=0, line=-0.5)
  }
  
    axis(side=2, at=c(0:4)/4, labels=NA, lwd=1, lwd.ticks=1)
    axis(side=2, at=c(0:2)/2, labels=percent(0:2/2), lwd=0, las=2, line=-0.5)
    mtext(side=2, line=2.5, text='survival')
  
  segments(x0=0, y0=0, y1=1, col=infect_col)
  par(xpd=T)
  points(x=0,y=1.0,pch=25,col=infect_col,bg=infect_col)
  for (dpi in surg_dpi) {
    points(x=dpi,y=1.05,pch=25,col=treat_col,bg=treat_col)
  }
  mtext(side=3, line=0.1, at=tpt, adj=0, text=paste0(tpt,' dpi \u2192'), cex=0.9)
  par(xpd=F)
  # if (panel == 1) {
  #   mtext('A-G', side=3, cex=2, adj=-0.08, line = 0.5)
  # }
  mtext(LETTERS[panel], side=2, las=2, cex=2, line=4)
  panel = panel + 1
  
}


# 
# 
# # go with this for now
# t_params = t_params_acm
# 
# t_data = sqldf("
# select   s.animal, b.cohort, b.timepoint, b.treatment, p.color, s.dpi
# from     surviv s, blind b, params p
# where    s.animal = b.animal
# and      b.cohort = p.cohort
# and      b.expt = 'T'
# and      s.acm
# ;")
# 
# t_cohs = sqldf("
# select   cohort, timepoint, color, avg(dpi) mean_dpi, stdev(dpi) sd_dpi, count(*) n_dpi
# from     t_data
# group by 1, 2, 3
# order by 2
# ;")
# t_cohs = calculate_ci(t_cohs)
# 
# saline_mean = round(sqldf("
# select   avg(dpi) mean
# from     surviv s, blind b
# where    s.animal = b.animal
# and      b.expt = 'T'
# and      s.acm
# and      b.treatment = 'saline'
# ;")$mean)
# saline_col = params$color[params$cohort == 'T-PRE7-PBS']
# 
# t_intervention_lims = c(-14, 170)
# t_surv_lims = c(0, 510)
# ylims = c(0,550)
# xlabs = c(unique(blind$timepoint[blind$expt=='T' & !is.na(blind$timepoint)]), round(saline_mean))
# 
# par(mar=c(3,4,3,2))
# plot(NA, NA, xlim=t_intervention_lims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
# mtext(side=3, line=1.5, text='efficacy by treatment timepoint')
# axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
# axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
# axis(side=1, at=0:10*50, lwd=0, lwd.ticks=0, line=-0.5)
# # axis(side=1, at=xlabs, labels=NA, lwd=0, lwd.ticks=1, tck=-0.0125, line=0.5)
# # for (xlab in xlabs[1:7]) { # loop to make sure they all get plotted no matter what
# #   axis(side=1, at=xlab, labels=xlab, lwd=0, lwd.ticks=0, cex.axis=1)
# # }
# # axis(side=1, at=saline_mean, labels=saline_mean, lwd=0, lwd.ticks=0, cex.axis=1, col.axis=saline_col)
# mtext(side=1, line=1.5, text='timepoint of chronic treatment initiation (dpi)')
# abline(v=0, lwd=1, lty=3)
# mtext(side=3, at=0, line=0.25, text='inoculation', cex=0.8)
# par(xpd=T)
# points(x=0, y=max(t_surv_lims)*1.05, pch=25, bg='black', col='black')
# segments(x0=saline_mean, y0=-25, y1=max(t_surv_lims), lwd=1.5, lty=3, col=alpha(saline_col,.5))
# segments(x0=0, x1=saline_mean+3, y0=saline_mean, lwd=1.5, lty=3, col=alpha(saline_col,.5))
# text(x=saline_mean, y=max(t_surv_lims), pos=3, labels='saline mean', cex=0.8, col=saline_col)
# text(x=150, y=150, pos=4, labels='saline mean', cex=0.8, col=saline_col)
# par(xpd=F)
# axis(side=2, at=(0:50)*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
# axis(side=2, at=(0:5)*100, lwd=0, lwd.ticks=1, las=2, tck=-0.05)
# mtext(side=2, line=2.5, text='survival (dpi)')
# # raw data points
# points(x=t_data$timepoint, y=t_data$dpi, col=alpha(t_data$color, 0.5), pch=19, lwd=0)
# # means & 95%CIs
# arrows(x0=t_cohs$timepoint, y0=t_cohs$l95, y1=t_cohs$u95, col=t_cohs$color, angle=90, code=3, length=0.1, lwd=1.5)
# segments(x0=t_cohs$timepoint-7, x1=t_cohs$timepoint+7, y0=t_cohs$mean_dpi, col=t_cohs$color, lwd=1.5)
# leg = subset(params, expt=='T' & grepl('T-PRE7',cohort))
# par(xpd=T)
# legend(x=30,y=750,horiz=T,legend=leg$display,col=leg$color,lwd=1.5,text.col=leg$color,text.font=2,bty='n',cex=1)
# par(xpd=F)
# mtext(LETTERS[panel], side=3, cex=2, adj = -0.05, line = 0.5)
# panel = panel + 1

# is there a significant difference among the 5 early or among the 2 late cohorts?
t_grouped_surv = sqldf("
select   s.animal, b.cohort, b.timepoint, b.treatment, p.color, s.dpi, s.acm,
         case when b.timepoint in (105, 120) and b.treatment = 'active ASO 6' then 'late ASO'
              when b.timepoint in (-7,1,28,54,78) and b.treatment = 'active ASO 6' then 'early ASO'
              when b.treatment = 'saline' then 'saline' end meta_cohort
from     surviv s, blind b, params p
where    s.animal = b.animal
and      b.cohort = p.cohort
and      b.expt = 'T'
;")

# special plot parameters for T

t_pparams = data.frame(meta_cohort=c('early ASO','late ASO','saline'),
                       alt_disp=c('early ASO (≤78 dpi)', 'late ASO (≥105 dpi)', 'saline'),
                      color=c(rep(params$color[params$cohort=='T-PRE7-ACTIVE6'],2),params$color[params$cohort=='T-PRE7-PBS']),
                      lty=c(5,3,1))
t_cohorts = sqldf("
           select   ts.cohort, ts.meta_cohort, tp.color, tp.lty
           from     t_grouped_surv ts, t_pparams tp
           where    ts.meta_cohort = tp.meta_cohort
           group by 1, 2, 3, 4
           order by 1, 2, 3, 4
           ;")

ks.test(t_grouped_surv$dpi[t_grouped_surv$timepoint == 105 & t_grouped_surv$treatment=='active ASO 6'], t_grouped_surv$dpi[t_grouped_surv$timepoint == 120 & t_grouped_surv$treatment=='active ASO 6'], alternative='two.sided')
m = lm(dpi ~ timepoint, data=subset(t_grouped_surv, meta_cohort=='early ASO'))
summary(m)
anova = aov(dpi ~ timepoint, data=subset(t_grouped_surv, meta_cohort=='early ASO'))
summary(anova)
# nope, no differences
# ok, prepare simplified survival curve with just 3 curves
sf = survfit(Surv(dpi, acm) ~ meta_cohort, data=t_grouped_surv)
sf$meta_cohort = gsub('meta_cohort=','',names(sf$strata))
sf$col = t_pparams$color[match(sf$meta_cohort, t_pparams$meta_cohort)]
sf$lty = t_pparams$lty[match(sf$meta_cohort, t_pparams$meta_cohort)]

# add in the G RML cohort as a control for comparison
g_rml_surv = sqldf("
select   s.animal, b.cohort, p.color, s.dpi, s.acm
from     surviv s, blind b, params p
where    s.animal = b.animal
and      b.cohort = p.cohort
and      b.expt = 'G' and b.strain = 'RML'
;")
gsf = survfit(Surv(dpi, acm) ~ cohort, data=g_rml_surv)
gsf$cohort = gsub('cohort=','',names(gsf$strata))
#gsf$col = params$color2[match(gsf$cohort, params$cohort)]
#gsf$lty = params$lty2[match(gsf$cohort, params$cohort)]
#gsf$lty[gsf$lty==3] = 2
gsf$lty = c(1,2)
gc_col = '#E6550D'
gsf$col = gc_col

par(mar=c(3,6,3,2))
plot(NA, NA, ann=FALSE, axes=FALSE, xlim=xlims, ylim=c(0,1.05), xaxs='i', yaxs='i')
par(new=T)
plot(gsf, ann=FALSE, axes=FALSE, xlim=xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=gsf$col, lty=gsf$lty, lwd=default_lwd)
par(new=T)
plot(sf, ann=FALSE, axes=FALSE, xlim=xlims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$col, lty=sf$lty, lwd=default_lwd)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:5*100, lwd=0, lwd.ticks=0, line=-0.5)

axis(side=2, at=c(0:4)/4, labels=NA, lwd=1, lwd.ticks=1)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, las=2, line=-0.5)
mtext(side=2, line=2.5, text='survival')

mtext(side=1, line=1.5, text='days post-infection')
mtext(side=3, line=0.5, text='combined survival')
segments(x0=0, y0=0, y1=1, col=infect_col)
par(xpd=T)
points(x=0,y=1.0,pch=25,col=infect_col,bg=infect_col)
par(xpd=F)
par(xpd=T)
legend(x=435,y=1.15,legend=rev(t_pparams$meta_cohort), col=rev(t_pparams$color), lty=rev(t_pparams$lty), lwd=default_lwd/2, text.col=rev(t_pparams$color), text.font=2, bty='n', cex=1, title.col='black', title='treatment')
legend(x=435,y=.78,legend=c('wild-type','Prnp+/-'), col=gsf$col, lty=gsf$lty, lwd=default_lwd/2, text.col=gsf$col, text.font=2, bty='n', cex=1, title.col='black', title='genetic control')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj=-0.05, line = 0.5)
panel = panel + 1

dev.off()










#### FIGURE S7 - THE SUPP HALF OF FIGURE 5
imgsave(paste('display_items/figure-s7.',imgmode,sep=''),width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1:8),byrow=T,nrow=2)
layout(layout_matrix, widths=c(.75,1,1,1))
panel = 1

# meta cohort legend panel
meta_cohorts = sqldf("select meta_cohort, color, lty from t_cohorts group by 1 order by 1 desc;")
par(mar=c(1,1,1,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', axes=F, ann=F)
legend(x=0,y=.75,legend=meta_cohorts$meta_cohort,col=meta_cohorts$color,text.col=meta_cohorts$color,text.font=2,lwd=3,lty=meta_cohorts$lty, bty='n',
       title = 'intervention group', title.col='black')


par(mar=c(3,3,3,1))
# redo all the T grouping based on *meta* cohorts rather than individual cohorts, 
# and require >4 mice (these are large meta-cohorts; gets noisy once you are down to 1-2 mice)
blind$meta_cohort = t_cohorts$meta_cohort[match(blind$cohort, t_cohorts$cohort)]
tw = sqldf("
            select   b.expt, b.meta_cohort, w.dpi, avg(w.weight) mean_weight, stdev(w.weight) sd_weight, count(*) n_weight
            from     weights w, blind b
            where    w.animal = b.animal
            and      dpi > 80
            and      b.expt = 'T'
            group by 1, 2, 3
            having   count(*) > 4
            order by 1, 2, 3
            ;")
tw = calculate_ci(tw)
tw$color = t_pparams$color[match(tw$meta_cohort, t_pparams$meta_cohort)]
tw$lty = t_pparams$lty[match(tw$meta_cohort, t_pparams$meta_cohort)]


plot(NA, NA, xlim=c(90,510), ylim=c(15,35), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:10*50, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days post-infection')
axis(side=2, at=0:10*5, las=2)
mtext(side=2, line=2, text='weight (g)')
abline(h=1, lty=3)
for (coh in unique(tw$meta_cohort)) {
  s = subset(tw, meta_cohort==coh)
  points(s$dpi, s$mean_weight, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
  polygon(x=c(s$dpi[!is.na(s$sd_weight)],rev(s$dpi[!is.na(s$sd_weight)])), y=c(s$u95[!is.na(s$sd_weight)],rev(s$l95[!is.na(s$sd_weight)])), col=alpha(s$color[!is.na(s$sd_weight)],ci_alpha), border=NA)
}
mtext(side=3, line=0.5, text='weights')
par(xpd=T)
legend(x=510,y=1.5,legend=rev(t_pparams$meta_cohort), col=rev(t_pparams$color), lty=rev(t_pparams$lty), lwd=1.5, text.col=rev(t_pparams$color), text.font=2, bty='n', cex=1.2)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


# "meta cohort" breakdown for behaviorals
tv = sqldf("
           select   b.expt, b.meta_cohort, v1.dpi, avg(n_signs) mean_signs, stdev(n_signs) sd_signs, count(*) n_animals
           from     v1, blind b
           where    b.animal = v1.animal
           group by 1, 2, 3
           having   count(*) > 4
           order by 1, 2, 3
           ;")
tv = calculate_ci(tv)
tv$color = t_pparams$color[match(tv$meta_cohort, t_pparams$meta_cohort)]
tv$lty = t_pparams$lty[match(tv$meta_cohort, t_pparams$meta_cohort)]

plot(NA, NA, xlim=c(90,510), ylim=c(0,6), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:10*50, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days post-infection')
axis(side=2, at=0:8, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2, text='symptoms per animal')
mtext(side=3, line=0.5, text='behaviorals')
for (coh in unique(tv$meta_cohort)) {
  b = subset(tv, meta_cohort==coh)
  points(b$dpi, b$mean_signs, type='l', lwd=strain_lwd, col=b$color, lty=b$lty)
  polygon(x=c(b$dpi[b$n_animals > 1],rev(b$dpi[b$n_animals > 1])), y=c(b$u95[b$n_animals > 1],rev(b$l95[b$n_animals > 1])), col=alpha(b$color,ci_alpha), border=NA)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# "meta cohort" breakdown for nests
t_nest_cages = sqldf("
                     select   b.expt, b.cage, b.meta_cohort
                     from     blind b
                     where    b.cage not in (select cage from mixed_cages)
                     and      b.cohort not in (select cohort from cage_exclude_cohorts) 
                     and      b.expt = 'T'
                     group by 1, 2, 3
                     order by 1, 2, 3
                     ;")

tn = sqldf("
           select   nc.expt, nc.meta_cohort, n.dpi, avg(n.comb) mean_nest, stdev(n.comb) sd_nest, count(*) n_cage
           from     nests n, t_nest_cages nc
           where    n.cage = nc.cage
           group by 1, 2, 3
           order by 1, 2, 3
           ;")

tn = calculate_ci(tn)
tn$color = t_pparams$color[match(tn$meta_cohort, t_pparams$meta_cohort)]
tn$lty = t_pparams$lty[match(tn$meta_cohort, t_pparams$meta_cohort)]

plot(NA, NA, xlim=c(90,510), ylim=c(0,2.1), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:10*50, lwd=0, line=-0.5)
axis(side=2, at=c(0:4)/2, lwd=1, lwd.ticks=1, las=2)
mtext(side=1, line=1.5, text='days post-infection')
mtext(side=2, line=2, text='mean nest score')
mtext(side=3, line=0.5, text='nests')
for (coh in unique(tn$meta_cohort)) {
  s = subset(tn, meta_cohort==coh)
  points(s$dpi, s$mean_nest, col=s$color, pch=20, type='l', lwd=strain_lwd, lty=s$lty)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


# compare early cohorts
early_breakdown = data.frame(cohort=c('T-PRE7-ACTIVE6','T-1DPI-ACTIVE6','T-28DPI-ACTIVE6','T-54DPI-ACTIVE6','T-78DPI-ACTIVE6'),
                             display=c('-7','1','28','54','78'),
                             color=c("#008837", "#A6DBA0", "#575757", "#C2A5CF", "#7B3294"),
                             lty=c(1,1,1,1,1))


# breakdown legend panel
par(mar=c(1,1,1,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', axes=F, ann=F)
legend(x=0,y=.75,legend=paste(early_breakdown$display,'dpi'),col=early_breakdown$color,text.col=early_breakdown$color,text.font=2,lwd=3,lty=early_breakdown$lty, bty='n',
       title = 'early intervention\ntimepoint', title.col='black')

par(mar=c(3,3,3,1))

# code to use individual weights:
# tw = sqldf("
#            select   w.*, e.cohort, e.color, e.lty
#            from     weights w, blind b, early_breakdown e
#            where    w.animal = b.animal
#            and      b.cohort = e.cohort
#            and      b.expt = 'T'
#            ;")
# grouped weights are easier to see:
tw = sqldf("
           select   e.cohort, w.dpi, e.color, e.lty, avg(w.weight) mean_weight, stdev(w.weight) sd_weight, count(*) n_weight
           from     weights w, blind b, early_breakdown e
           where    w.animal = b.animal
           and      b.cohort = e.cohort
           and      b.expt = 'T'
group by 1, 2, 3, 4 order by 1, 2
           ;")

plot(NA, NA, xlim=c(90,510), ylim=c(15,35), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:10*50, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days post-infection')
axis(side=2, at=0:10*5, labels=0:10*5, las=2)
mtext(side=2, line=2.0, text='weight (g)')
for (coh in unique(tw$cohort)) {
  s = subset(tw, coh==cohort)
  points(s$dpi, s$mean_weight, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='weights')
par(xpd=T)
legend(x=510,y=1.5,legend=paste(early_breakdown$display,'dpi'), col=early_breakdown$color, lwd=1.5, text.col=early_breakdown$color, text.font=2, bty='n', cex=1.2)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


m = lm(weight ~ dpi + cohort, data=tw)
summary(m)
# "early breakdown" by timepoint for behaviorals
tve = subset(v, cohort %in% early_breakdown$cohort)
tve$color = early_breakdown$color[match(tve$cohort, early_breakdown$cohort)]
tve$lty = early_breakdown$lty[match(tve$cohort, early_breakdown$cohort)]

plot(NA, NA, xlim=c(90,510), ylim=c(0,6), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:10*50, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days post-infection')
axis(side=2, at=0:6, labels=0:6, las=2)
mtext(side=2, line=2, text='symptoms per animal')
for (coh in unique(tve$cohort)) {
  s = subset(tve, cohort==coh)
  points(s$dpi, s$mean_signs, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='behaviorals')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

# "early breakdown" by timepoint for nests
tne = subset(n, cohort %in% early_breakdown$cohort)
tne$color = early_breakdown$color[match(tne$cohort, early_breakdown$cohort)]
tne$lty = early_breakdown$lty[match(tne$cohort, early_breakdown$cohort)]

plot(NA, NA, xlim=c(90,510), ylim=c(0,2.1), ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=0:10*50, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days post-infection')
axis(side=2, at=0:2, labels=0:2, las=2)
mtext(side=2, line=2, text='mean nest score')
for (coh in unique(tne$cohort)) {
  s = subset(tne, cohort==coh)
  points(s$dpi, s$mean_nest, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
}
mtext(side=3, line=0, text='nests')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

dev.off()


















# 
# #### taller aspect ratio versions for powerpoints
# imgsave(paste('display_items/figure-p1.',imgmode,sep=''),width=7.5*resx,height=3.5*resx,res=resx)
# par(mar=c(4,5,4,10))
# plot(NA, NA, xlim=t_intervention_lims, ylim=t_surv_lims, xaxs='i', yaxs='i', axes=F, ann=F)
# axis(side=1, at=t_intervention_lims, labels=NA, lwd=1, lwd.ticks=0)
# xlabs = c(unique(blind$timepoint[blind$expt=='T' & !is.na(blind$timepoint)]), round(saline_mean))
# axis(side=1, at=xlabs, labels=NA, lwd=0, lwd.ticks=1, cex.axis=1)
# for (xlab in xlabs) {
#   axis(side=1, at=xlab, labels=xlab, lwd=0, lwd.ticks=0, cex.axis=1)
# }
# mtext(side=1, line=2.5, text='timepoint of chronic treatment initiation (dpi)')
# abline(v=0, lwd=1, lty=3)
# mtext(side=3, at=0, line=0, text='inoculation')
# abline(v=saline_mean, lwd=1, lty=3)
# mtext(side=3, at=saline_mean, line=0.5, text='saline mean\nsurvival')
# axis(side=2, at=(0:5)*100, lwd=0, lwd.ticks=1, las=2)
# mtext(side=2, line=3.5, text='survival (dpi)')
# abline(h=saline_mean, lwd=1.5, lty=3, col=alpha(saline_col,.5))
# mtext(side=4, at=saline_mean, line=0.25, las=2, adj=0, col=saline_col, font=1, cex=1, text='saline mean')
# # mtext(side=3, line=0, text='mortality')
# # raw data points
# points(x=t_data$timepoint, y=t_data$dpi, col=alpha(t_data$color, 0.5), pch=19, lwd=0)
# # means & 95%CIs
# arrows(x0=t_cohs$timepoint, y0=t_cohs$l95, y1=t_cohs$u95, col=t_cohs$color, angle=90, code=3, length=0.1, lwd=1.5)
# segments(x0=t_cohs$timepoint-7, x1=t_cohs$timepoint+7, y0=t_cohs$mean_dpi, col=t_cohs$color, lwd=1.5)
# leg = subset(params, expt=='T' & grepl('T-PRE7',cohort))
# par(xpd=T)
# legend(x=170,y=600,legend=leg$display,col=leg$color,lwd=1.5,text.col=leg$color,text.font=2,bty='n',cex=1.2)
# par(xpd=F)
# dev.off()
# 
# 
# #### taller aspect ratio versions for powerpoints
# imgsave(paste('display_items/figure-p2.',imgmode,sep=''),width=7.5*resx,height=3.5*resx,res=resx)
# 
# sf = survfit(Surv(dpi, acm) ~ meta_cohort, data=t_surv)
# sf$meta_cohort = gsub('meta_cohort=','',names(sf$strata))
# sf$col = t_pparams$color[match(sf$meta_cohort, t_pparams$meta_cohort)]
# sf$lty = t_pparams$lty[match(sf$meta_cohort, t_pparams$meta_cohort)]
# par(mar=c(4,5,4,10))
# plot(sf, ann=FALSE, axes=FALSE, xlim=t_surv_lims, ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$col, lty=sf$lty, lwd=1.5)
# axis(side=1, at=c(0:5)*100, lwd=0, lwd.ticks=1)
# axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
# abline(h=0)
# abline(v=0)
# mtext(side=2, line=3.5, text='proportion surviving')
# mtext(side=1, line=2.5, text='days post-infection')
# mtext(side=3, line=0, text='mortality')
# par(xpd=T)
# legend(x=450,y=1.25,legend=rev(t_pparams$alt_disp), col=rev(t_pparams$color), lty=rev(t_pparams$lty), lwd=1.5, text.col=rev(t_pparams$color), text.font=2, bty='n', cex=1.1)
# par(xpd=F)
# dev.off()















# V cohort - table S2 and figure 6C+
v_params_acm = sqldf("
                     select   p.x, b.timepoint,
                     avg(surg.dpi) actual_dpi_treatment,
                     avg(  case when b.treatment='saline'       and s.acm then s.dpi else null end) saline_mean_dpi,
                     stdev(case when b.treatment='saline'       and s.acm then s.dpi else null end) saline_sd_dpi,
                     sum(  case when b.treatment='saline'       and s.acm then 1 else 0 end) saline_n,
                     avg(  case when b.treatment='active ASO 6' and s.acm then s.dpi else null end) active_mean_dpi,
                     stdev(case when b.treatment='active ASO 6' and s.acm then s.dpi else null end) active_sd_dpi,
                     sum(  case when b.treatment='active ASO 6' and s.acm then 1 else 0 end) active_n
                     from     surviv s, blind b, surgery surg, params p
                     where    s.animal = b.animal
                     and      s.animal = surg.animal
                     and      surg.event = 'icv1'
                     and      b.expt = 'V'
                     and      b.cohort = p.cohort
                     and      b.timepoint is not null
                     group by 1, 2
                     order by 1
                     ;")

v_params_acm$delta = paste('+',percent(v_params_acm$active_mean_dpi / v_params_acm$saline_mean_dpi - 1),sep='')
v_params_acm$saline_itime = paste(formatC(v_params_acm$saline_mean_dpi,format='f',digits=0),
                                  '±',formatC(v_params_acm$saline_sd_dpi,format='f',digits=0),sep='')
v_params_acm$active_itime = paste(formatC(v_params_acm$active_mean_dpi,format='f',digits=0),
                                  '±',formatC(v_params_acm$active_sd_dpi,format='f',digits=0),sep='')

table_s4 = v_params_acm[,c('timepoint','saline_itime','saline_n','active_itime','active_n','delta')]
write.table(table_s4, 'display_items/table-s4.tsv', sep='\t', col.names=T, row.names=F, quote=F)

## examine data
# View(sqldf("select * from blind b, surviv s where b.animal = s.animal and b.expt ='V' order by b.animal;"))
# View(sqldf("select * from blind b, surviv s where b.animal = s.animal and b.timepoint in (-7,1,28,56,78) and b.expt ='T' order by b.animal;"))

v_surv = sqldf("
               select   b.timepoint, b.treatment, s.animal, s.dpi, s.acm, b.cohort, p.color, p.lty
               from     surviv s, blind b, params p
               where    s.animal = b.animal
               and      b.cohort = p.cohort
               and      b.expt = 'V' and b.timepoint is not null
               order by 1, 2, 4
               ;")

# in animals that lived longer, how much longer did they live?
saline_mean = mean(v_surv$dpi[v_surv$timepoint %in% c(132, 143) & v_surv$treatment=='saline' & v_surv$acm])
extended_mean = mean(v_surv$dpi[v_surv$timepoint %in% c(132, 143) & v_surv$treatment=='active ASO 6' & v_surv$acm & v_surv$dpi > 1.1*saline_mean])
extended_mean - saline_mean
extended_proportion = mean(v_surv$dpi[v_surv$timepoint %in% c(132, 143) & v_surv$treatment=='active ASO 6' & v_surv$acm] > 1.1*saline_mean)
sum(v_surv$dpi[v_surv$timepoint %in% c(132, 143) & v_surv$treatment=='active ASO 6' & v_surv$acm] > 1.1*saline_mean)
sum(v_surv$dpi[v_surv$timepoint %in% c(132, 143) & v_surv$treatment=='active ASO 6' & v_surv$acm] <= 1.1*saline_mean)

#### FIGURE 6: V cohort - very late intervention
imgsave(paste('display_items/figure-6.',imgmode,sep=''),width=6.5*resx,height=8*resx,res=resx)

#layout_matrix = matrix(1:8, nrow=4, byrow=T)
#layout_matrix = matrix(c(1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,8,8,8,9,9,9,10,11,11,11,12,12,12),nrow=4,byrow=T)

layout_matrix = matrix(c(1:20,rep(21,5)), nrow=5, byrow=T)
layout(layout_matrix, widths=c(.5,1,1,1,1), heights=c(1,1,1,1,.15))

v_surv = sqldf("
               select   b.timepoint, b.treatment, s.animal, s.dpi, s.acm, b.cohort, p.color, p.lty
               from     surviv s, blind b, params p
               where    s.animal = b.animal
               and      b.cohort = p.cohort
               and      b.expt = 'V' and b.timepoint is not null
               order by 1, 2, 4
               ;")

# use raw, individual weights for cohort v - because we want to show the mice are symptomatic before treatment
vw = sqldf("
           select   b.timepoint, b.treatment, w.animal, w.dpi, w.weight, b.cohort, p.color, p.lty
           from     weights w, blind b, params p
           where    w.animal = b.animal
           and      b.cohort = p.cohort
           and      b.expt = 'V'
           ;")


panel = 1

for (tpt in unique(blind$timepoint[blind$expt=='V' & !is.na(blind$timepoint)])) {

  # blank panel for timepoint
  par(mar=c(1,1,1,1))
  plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ann=FALSE,xaxs='i',yaxs='i')
  par(xpd=T)
  text(x=.5, y=.5, labels=paste(tpt,'\ndpi'), cex=1.5)
  
  par(xpd=F)
  
  # set margins for rest of panels
  if (tpt == 120) {
    par(mar=c(1,3,4,1))
  } else if (tpt == 156) {
    par(mar=c(4,3,1,1))
  } else {
    par(mar=c(2.5,3,2.5,1))
  }

  if (tpt < 156) {

    vsurv = survfit(Surv(dpi, acm) ~ cohort, data=subset(v_surv, timepoint==tpt))
    vsurv$cohort = gsub('cohort=','',names(vsurv$strata))
    vsurv$color = params$color[match(vsurv$cohort, params$cohort)]
    vsurv$lty = params$lty[match(vsurv$cohort, params$cohort)]

    plot(vsurv, ann=FALSE, axes=FALSE, xlim=c(0,300), ylim=c(0,1.05), xaxs='i', yaxs='i', lwd=default_lwd, col=vsurv$color, lty=vsurv$lty)
    axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
    axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
    axis(side=1, at=0:10*50, lwd=0, line=-0.5)
    axis(side=2, at=c(0:4)/4, labels=NA, lwd=1, lwd.ticks=1)
    axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, las=2, line=-0.5)
    mtext(side=2, line=2.25, text='survival')
    # abline(v=tpt)
    par(xpd=T)
    points(x=tpt,y=1.05,pch=25,col='black',bg='black')
    par(xpd=F)
  }
  if (tpt == 156) {
    # make a special 156 dpi V plot that shows how many died before treatment group assignment
    vsurv = sqldf("
                  select   s.animal, s.dpi, s.acm, b.cohort, b.expt, p.color, p.lty
                  from     surviv s, blind b, params p
                  where    b.expt = 'V'
                  and      s.animal = b.animal
                  and      b.cohort = p. cohort
                  and      b.cage in (select cage from blind where cohort in ('V-156-PBS','V-156-ACTIVE6'))
                  and      s.dpi > 0
                  order by 4, 1
                  ;")
    vsurv$status1 = vsurv$dpi <= tpt
    vs1 = survfit(Surv(dpi, status1) ~ 1, data=vsurv)
    vs2 = survfit(Surv(dpi, acm) ~ cohort, data=subset(vsurv, dpi > 156))
    vs2$cohort = gsub('cohort=','',names(vs2$strata))
    vs2$color = params$color[match(vs2$cohort, params$cohort)]
    vs2$lty = params$lty[match(vs2$cohort, params$cohort)]
    
    plot(NA, NA, ann=FALSE, axes=FALSE, xlim=c(0,300), ylim=c(0,1.05), xaxs='i', yaxs='i')
    axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
    axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
    axis(side=1, at=0:10*50, lwd=0, line=-0.5)
    mtext(side=1, line=1.5, text='dpi')
    axis(side=2, at=c(0:4)/4, labels=NA, lwd=1, lwd.ticks=1)
    axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, las=2, line=-0.5)
    mtext(side=2, line=2.25, text='survival')
    abline(v=tpt)
    par(xpd=T)
    points(x=tpt,y=1.05,pch=25,col='black',bg='black')
    par(xpd=F)
    lines(vs1, lwd=default_lwd, col='#FF2017', lty=1, conf.int=F, xmax=156)
    prop = sum(vs1$n.censor)/sum(vs1$n) # proportion that survived to participate in the study
    x_aso = c(156,vs2$time[1:vs2$strata['cohort=V-156-ACTIVE6']])
    x_aso = rep(x_aso, each=2)
    y_aso = prop*c(1,vs2$surv[1:vs2$strata['cohort=V-156-ACTIVE6']])
    y_aso = c(prop, rep(y_aso, each=2))[1:8]
    lines(x_aso, y_aso, lwd=default_lwd, col=params$color[params$cohort=='V-156-ACTIVE6'], lty=params$lty[params$cohort=='V-156-ACTIVE6'])
    x_pbs = c(156,vs2$time[(vs2$strata['cohort=V-156-ACTIVE6']+1):(vs2$strata['cohort=V-156-ACTIVE6']+vs2$strata['cohort=V-156-PBS'])])
    x_pbs = rep(x_pbs, each=2)
    y_pbs = prop*c(1,vs2$surv[(vs2$strata['cohort=V-156-ACTIVE6']+1):(vs2$strata['cohort=V-156-ACTIVE6']+vs2$strata['cohort=V-156-PBS'])])
    y_pbs = c(prop, rep(y_pbs, each=2))[1:12]
    lines(x_pbs, y_pbs, lwd=default_lwd, col=params$color[params$cohort=='V-156-PBS'], lty=params$lty[params$cohort=='V-156-PBS'])
  }
  if (tpt==120) {
    mtext(side=3, line=1, text='mortality')
  }
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.5, line = 0.5)
  panel = panel + 1
  
  hvw = subset(vw, timepoint==tpt)
  plot(NA, NA, xlim=c(90,300), ylim=c(15,35), ann=F, axes=F, xaxs='i', yaxs='i')
  axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=1, at=0:10*50, lwd=0, line=-0.5)
  if (tpt == 156) {
    mtext(side=1, line=1.5, text='dpi')
  }
  axis(side=2, at=0:7*5, labels=0:7*5, las=2)
  mtext(side=2, line=2, text='weight (g)')
  abline(h=1, lty=3)
  abline(v=tpt)
  par(xpd=T)
  par(xpd=T)
  points(x=tpt,y=35,pch=25,col='black',bg='black')
  par(xpd=F)
  par(xpd=F)
  for (anml in unique(hvw$animal)) {
    s = subset(hvw, animal==anml)
    points(s$dpi, s$weight, type='l', lwd=0.75, col=s$color, lty=s$lty)
  }
  if (tpt==120) {
    mtext(side=3, line=1, text='\nindividual\nweights')
  }
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.5, line = 0.5)
  panel = panel + 1
  
  # symptom count
  vv = subset(v, expt == 'V' & cohort %in% blind$cohort[blind$timepoint==tpt])
  vv$lty = params$lty[match(vv$cohort,params$cohort)]
  vv$color = params$color[match(vv$cohort,params$cohort)]
  plot(NA, NA, xlim=c(90,300), ylim=c(0,6), ann=F, axes=F, xaxs='i', yaxs='i')
  axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=1, at=0:10*50, lwd=0, line=-0.5)
  if (tpt == 156) {
    mtext(side=1, line=1.5, text='dpi')
  }
  axis(side=2, at=0:6, labels=0:6, las=2)
  mtext(side=2, line=2, text='mean symptoms')
  for (coh in unique(vv$cohort)) {
    s = subset(vv, cohort==coh)
    points(s$dpi, s$mean_signs, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
  }
  abline(v=tpt)
  par(xpd=T)
  points(x=tpt,y=6,pch=25,col='black',bg='black')
  par(xpd=F)
  if (tpt==120) {
    mtext(side=3, line=1, text='behaviorals')
  }
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.5, line = 0.5)
  panel = panel + 1
  
  # nests
  vn = subset(n, expt=='V'  & cohort %in% blind$cohort[blind$timepoint==tpt])
  vn$color = params$color[match(vn$cohort, params$cohort)]
  vn$lty = params$lty[match(vn$cohort, params$cohort)]
  plot(NA, NA, xlim=c(90,300), ylim=c(0,2.1), ann=F, axes=F, xaxs='i', yaxs='i')
  axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=0:10*50, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=1, at=0:10*50, lwd=0, line=-0.5)
  if (tpt == 156) {
    mtext(side=1, line=1.5, text='dpi')
  }
  axis(side=2, at=0:2, labels=0:2, las=2)
  mtext(side=2, line=2, text='mean score')
  for (coh in unique(vn$cohort)) {
    s = subset(vn, cohort==coh)
    points(s$dpi, s$mean_nest, type='l', lwd=strain_lwd, col=s$color, lty=s$lty)
    points(s$dpi, s$mean_nest, pch=20, cex=0.6, col=s$color)
  }
  abline(v=tpt)
  par(xpd=T)
  points(x=tpt,y=2.1,pch=25,col='black',bg='black')
  par(xpd=F)
  if (tpt==120) {
    mtext(side=3, line=1, text='nests')
  }
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.5, line = 0.5)
  panel = panel + 1
  
}

par(mar=c(0,0,0,0))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', ann=F, axes=F)
leg = subset(params, expt=='V' & grepl('V-120',cohort))
legend(x=0.25,y=1.5,horiz=T,legend=leg$display,col=leg$color,text.col=leg$color,lwd=3,bty='n',cex=1.75)

dev.off()


# how much was survival increased among those who did survive longer at 132 and 143?
v_saline_mean = mean(v_surv$dpi[v_surv$treatment=='saline' & v_surv$timepoint %in% c(132,143)])
v_aso_outlived_mean = mean(v_surv$dpi[v_surv$treatment=='active ASO 6' & v_surv$timepoint %in% c(132,143) & v_surv$dpi > v_saline_mean*1.1])
v_aso_outlived_mean - v_saline_mean

# how many individual animals had weight loss prior to 132 or 143
sqldf("
select   tmt.timepoint,
         sum(case when tmt.weight < prev.max_wt then 1 else 0 end) n_that_had_lost_weight,
         avg(case when tmt.weight < prev.max_wt then 1 else 0 end) proportion_that_had_lost_weight
from     vw tmt, (select animal, max(weight) max_wt from vw group by 1) prev
where    tmt.animal = prev.animal
and      tmt.dpi = tmt.timepoint
and      tmt.timepoint in (132, 143)
group by 1
order by 1
;")











imgsave(paste('display_items/figure-3.',imgmode,sep=''),width=6.5*resx,height=4.9*resx,res=resx)

layout_matrix = matrix(1:6,nrow=3,byrow=F)
layout(layout_matrix)

expt = 'N'

master = read.table(paste('data/mri/',expt,'_survival.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
cohorts = read.table(paste('data/mri/',expt,'_cohorts.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
rr = read.table(paste('data/mri/',expt,'_rotarod.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
nfl = read.table(paste('data/mri/',expt,'_nfl.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
nests = read.table(paste('data/mri/',expt,'_nests.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
weights = read.table(paste('data/mri/',expt,'_weights.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
behavs = read.table(paste('data/mri/',expt,'_behavs.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
behavs_meta = read.table(paste('data/mri/',expt,'_behavs_meta.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

write.table(behavs_meta, 'display_items/table-s2.tsv', sep='\t', row.names=F, col.names=T, quote=F)

linecol = '#AAAAAA'


rr$cohort = master$cohort[match(rr$animal,master$animal)]

rr_baseline = sqldf("
                    select   cohort, animal, avg(latency) mean_baseline_latency
                    from     rr
                    where    dpi = -7
                    group by 1, 2
                    order by 1, 2
                    ")

rr_smry1 = sqldf("
                 select   r.cohort, r.animal, r.dpi, avg(r.latency)/b.mean_baseline_latency change_latency
                 from     rr r, rr_baseline b
                 where    r.cohort = b.cohort and r.animal = b.animal
                 group by 1, 2, 3
                 order by 1, 2, 3
                 ;")

rr_smry2 = sqldf("
                 select   cohort, dpi, avg(change_latency) mean_change_lat, stdev(change_latency) sd_change_lat, count(*) n_change_lat
                 from     rr_smry1
                 group by 1, 2
                 having   count(*) > 1
                 order by 1, 2
                 ;")
rr_smry2 = calculate_ci(rr_smry2)
rr_smry2$color = cohorts$color[match(rr_smry2$cohort, cohorts$cohort)]
rr_smry2$lty = cohorts$lty[match(rr_smry2$cohort, cohorts$cohort)]

# add P vals
rr_smry2$ks_p = as.numeric(NA)
for (i in 1:nrow(rr_smry2[rr_smry2$cohort=='RML',])) {
  dpi = rr_smry2$dpi[i]
  rml_rows  = rr_smry1$dpi == dpi & rr_smry1$cohort == 'RML'
  unin_rows = rr_smry1$dpi == dpi & rr_smry1$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(rr_smry1$change_latency[rml_rows], rr_smry1$change_latency[unin_rows], alternative='two.sided'))
  rr_smry2$ks_p[i] = ks_obj$p.value
}
rr_smry2$p_symb = ''
rr_smry2$p_symb[!is.na(rr_smry2$ks_p) & rr_smry2$ks_p < 0.05] = '*'
rr_smry2$p_symb[!is.na(rr_smry2$ks_p) & rr_smry2$ks_p < 0.01] = '**'
rr_smry2$p_symb[!is.na(rr_smry2$ks_p) & rr_smry2$ks_p < 0.001] = '***'

nfl$color = cohorts$color[match(nfl$cohort, cohorts$cohort)]
nfl_smry = sqldf("
                 select   cohort, timepoint, avg(nfl_pgml) mean_nfl, stdev(nfl_pgml) sd_nfl, count(*) n_nfl
                 from     nfl
                 group by 1, 2
                 order by 1, 2
                 ;")
nfl_smry = calculate_ci(nfl_smry)
nfl_smry$color = cohorts$color[match(nfl_smry$cohort, cohorts$cohort)]
nfl_smry$lty =     cohorts$lty[match(nfl_smry$cohort, cohorts$cohort)]

# add P vals
nfl_smry$ks_p = as.numeric(NA)
for (i in 1:nrow(nfl_smry[nfl_smry$cohort=='RML',])) {
  dpi = nfl_smry$timepoint[i]
  rml_rows = nfl$timepoint == dpi & nfl$cohort == 'RML'
  unin_rows = nfl$timepoint == dpi & nfl$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(nfl$nfl_pgml[rml_rows], nfl$nfl_pgml[unin_rows], alternative='two.sided'))
  nfl_smry$ks_p[i] = ks_obj$p.value
}
nfl_smry$p_symb = ''
nfl_smry$p_symb[!is.na(nfl_smry$ks_p) & nfl_smry$ks_p < 0.05] = '*'
nfl_smry$p_symb[!is.na(nfl_smry$ks_p) & nfl_smry$ks_p < 0.01] = '**'
nfl_smry$p_symb[!is.na(nfl_smry$ks_p) & nfl_smry$ks_p < 0.001] = '***'

nst_smry = sqldf("
                 select   cohort, dpi, avg(comb) mean_comb, stdev(comb) sd_comb, count(*) n_comb
                 from     nests
                 where    comb is not null
                 group by 1, 2
                 having   count(*) > 1
                 order by 1, 2
                 ;")
nst_smry = calculate_ci(nst_smry)
nst_smry$color = cohorts$color[match(nst_smry$cohort, cohorts$cohort)]
nst_smry$lty =     cohorts$lty[match(nst_smry$cohort, cohorts$cohort)]

# add P vals
nst_smry$ks_p = as.numeric(NA)
for (i in 1:nrow(nst_smry[nst_smry$cohort=='RML',])) {
  dpi = nst_smry$dpi[i]
  rml_rows  = nests$dpi == dpi & nests$cohort == 'RML'
  unin_rows = nests$dpi == dpi & nests$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(nests$comb[rml_rows], nests$comb[unin_rows], alternative='two.sided'))
  nst_smry$ks_p[i] = ks_obj$p.value
}
nst_smry$p_symb = ''
nst_smry$p_symb[!is.na(nst_smry$ks_p) & nst_smry$ks_p < 0.05] = '*'
nst_smry$p_symb[!is.na(nst_smry$ks_p) & nst_smry$ks_p < 0.01] = '**'
nst_smry$p_symb[!is.na(nst_smry$ks_p) & nst_smry$ks_p < 0.001] = '***'



weights$cohort = master$cohort[match(weights$animal, master$animal)]
weighall_timepoints = data.frame(dpi=unique(weights$dpi[weights$cohort=='uninoculated']))
weights$plot_dpi = weights$dpi
needmatch = !(weights$plot_dpi %in% weighall_timepoints$dpi)
weights$plot_dpi[needmatch] = weighall_timepoints$dpi[findInterval(x=weights$dpi[needmatch], vec=weighall_timepoints$dpi, left.open=T, rightmost.closed=T)+1]
# check that rounding up worked as expected:
# weights[weights$dpi != weights$plot_dpi,] # yes - looks good
weights_baseline = sqldf("
                         select   cohort, animal, wt baseline_wt
                         from     weights
                         where    dpi = 78
                         order by 1, 2
                         ;")
weight_change = sqldf("
                      select   w.cohort, w.animal, w.plot_dpi, w.wt/b.baseline_wt weight_change
                      from     weights w, weights_baseline b
                      where    w.animal = b.animal
                      ;")
weights_smry = sqldf("
                     select   w.cohort, w.plot_dpi, avg(w.wt/b.baseline_wt) mean_weight_change, stdev(w.wt/b.baseline_wt) sd_weight_change, count(*) n_wt_change
                     from     weights w, weights_baseline b
                     where    w.animal = b.animal
                     group by 1, 2
                     having   count(*) > 1 -- only include dpi with >1 mouse measured
                     order by 1, 2
                     ;")
weights_smry = calculate_ci(weights_smry)
weights_smry$color = cohorts$color[match(weights_smry$cohort, cohorts$cohort)]
weights_smry$lty =     cohorts$lty[match(weights_smry$cohort, cohorts$cohort)]

# add P vals
weights_smry$ks_p = as.numeric(NA)
for (i in 1:nrow(weights_smry[weights_smry$cohort=='RML',])) {
  dpi = weights_smry$plot_dpi[i]
  rml_rows  = weight_change$plot_dpi == dpi & weight_change$cohort == 'RML'
  unin_rows = weight_change$plot_dpi == dpi & weight_change$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(weight_change$weight_change[rml_rows], weight_change$weight_change[unin_rows], alternative='two.sided'))
  weights_smry$ks_p[i] = ks_obj$p.value
}
weights_smry$p_symb = ''
weights_smry$p_symb[!is.na(weights_smry$ks_p) & weights_smry$ks_p < 0.05] = '*'
weights_smry$p_symb[!is.na(weights_smry$ks_p) & weights_smry$ks_p < 0.01] = '**'
weights_smry$p_symb[!is.na(weights_smry$ks_p) & weights_smry$ks_p < 0.001] = '***'


# summarize behaviorals
smry_by_animal = sqldf("
                       select   b.dpi,
                       b.animal,
                       m.cohort,
                       sum(case when score > 0 then 1 else 0 end) n_signs
                       from     behavs b, master m
                       where    b.animal = m.animal
                       group by 1, 2, 3
                       order by 2, 1
                       ;")

smry_by_cohort = sqldf("
                       select   m.cohort, s.dpi, avg(s.n_signs) mean_signs, stdev(s.n_signs) sd_signs, count(*) n_animals
                       from     smry_by_animal s, master m
                       where    s.animal = m.animal
                       group by 1, 2
                       having   count(*) > 1
                       order by 1, 2
                       ;")
smry_by_cohort = calculate_ci(smry_by_cohort)
behavs_smry = smry_by_cohort
behavs_smry$color = cohorts$color[match(behavs_smry$cohort, cohorts$cohort)]
behavs_smry$lty = cohorts$lty[match(behavs_smry$cohort, cohorts$cohort)]
behavs_smry$ks_p = as.numeric(NA)
for (i in 1:nrow(behavs_smry[behavs_smry$cohort=='RML',])) {
  dpi = behavs_smry$dpi[i]
  rml_rows  = smry_by_animal$dpi == dpi & smry_by_animal$cohort == 'RML'
  unin_rows = smry_by_animal$dpi == dpi & smry_by_animal$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(smry_by_animal$n_signs[rml_rows], smry_by_animal$n_signs[unin_rows], alternative='two.sided'))
  behavs_smry$ks_p[i] = ks_obj$p.value
}
behavs_smry$p_symb = ''
behavs_smry$p_symb[!is.na(behavs_smry$ks_p) & behavs_smry$ks_p < 0.05] = '*'
behavs_smry$p_symb[!is.na(behavs_smry$ks_p) & behavs_smry$ks_p < 0.01] = '**'
behavs_smry$p_symb[!is.na(behavs_smry$ks_p) & behavs_smry$ks_p < 0.001] = '***'

panel = 1


xlims = c(85,190)
ylims = c(0,16)
par(mar=c(1,5,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5)
axis(side=2, at=0:3*5, las=2)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
mtext(side=2.5, line=2.5, text='signs per animal')
mtext(side=3, line=0.5, text='behaviorals')
for (coh in unique(behavs_smry$cohort)) { # do them in reverse order so pre-treatment is on top at 81 dpi
  subs = subset(behavs_smry, cohort==coh)
  color = subs$color[1]
  points(subs$dpi, subs$mean_signs, type='l', lwd=3, col=color, lty=subs$lty)
  points(subs$dpi, subs$mean_signs, pch=19, col=color)
  polysubs = subset(subs, !is.na(sd_signs)) 
  polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
}
text(x=behavs_smry$dpi[behavs_smry$cohort=='RML'], y=rep(15, sum(behavs_smry$cohort=='RML')), labels=behavs_smry$p_symb[behavs_smry$cohort=='RML'], cex=1.0)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

ylims = c(0,2.25)
par(mar=c(1,5,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
abline(v=0)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=0:5/2, lwd=1, lwd.ticks=1, tck=-0.025, las=2)
mtext(side=2, line=2.5, text='combined score')
mtext(side=3, line=0.5, text='nests')

for (coh in unique(nst_smry$cohort)) { # do them in reverse order so pre-treatment is on top at 81 dpi
  subs = subset(nst_smry, cohort==coh)
  color = subs$color[1]
  points(subs$dpi, subs$mean_comb, type='l', lwd=3, col=color, lty=subs$lty)
  points(subs$dpi, subs$mean_comb, pch=19, col=color)
  polysubs = subset(subs, !is.na(sd_comb)) 
  polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
}
text(x=nst_smry$dpi[nst_smry$cohort=='RML'], y=rep(max(ylims), sum(nst_smry$cohort=='RML')), labels=nst_smry$p_symb[nst_smry$cohort=='RML'])
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


ylims = c(.5,1.15)
par(mar=c(4,5,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
mtext(side=1, line=2.0, text='days post-infection')
abline(v=0)
abline(h=1, lty=3)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=0:20/10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=0:8/4, labels=NA, lwd=1, lwd.ticks=1, tck=-0.050)
axis(side=2, at=0:8/4, labels=percent(0:8/4 - 1), lwd=0, las=2, line=-0.5)
mtext(side=2, line=3.0, text='Δ body weight (%)')
mtext(side=3, line=0.5, text='weights')
for (coh in unique(weights_smry$cohort)) {
  s = subset(weights_smry, cohort==coh)
  points(s$plot_dpi, s$mean_weight_change, col=s$color, type='l', lwd=3, lty=s$lty)
  points(s$plot_dpi, s$mean_weight_change, col=s$color, pch=20)
  polygon(x=c(s$plot_dpi[!is.na(s$sd_weight_change)],rev(s$plot_dpi[!is.na(s$sd_weight_change)])), y=c(s$u95[!is.na(s$sd_weight_change)],rev(s$l95[!is.na(s$sd_weight_change)])), col=alpha(s$color[!is.na(s$sd_weight_change)],ci_alpha), border=NA)
}
par(xpd=T)
text(x=weights_smry$plot_dpi[weights_smry$cohort=='RML'], y=rep(max(ylims), sum(weights_smry$cohort=='RML')), labels=weights_smry$p_symb[weights_smry$cohort=='RML'], cex=0.75)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


xlims = c(-15,200)
ylims = c(.1,1.1)

par(mar=c(1,5,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
abline(v=0)
abline(h=1, lty=3)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=0:20/10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=0:8/4, labels=NA, lwd=1, lwd.ticks=1, tck=-0.050)
axis(side=2, at=0:8/4, labels=percent(0:8/4 - 1), lwd=0, las=2, line=-0.5)
mtext(side=2, line=3.5, text='Δ drop latency (%)')
mtext(side=3, line=0.5, text='rotarod')

for (coh in unique(rr_smry2$cohort)) { # do them in reverse order so pre-treatment is on top at 81 dpi
  subs = subset(rr_smry2, cohort==coh)
  color = subs$color[1]
  points(subs$dpi, subs$mean_change_lat, type='l', lwd=3, col=color, lty=subs$lty)
  points(subs$dpi, subs$mean_change_lat, pch=19, col=color)
  polysubs = subset(subs, !is.na(sd_change_lat)) 
  polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
}
par(xpd=T)
text(x=rr_smry2$dpi[rr_smry2$cohort=='RML'], y=rep(max(ylims), sum(rr_smry2$cohort=='RML')), labels=rr_smry2$p_symb[rr_smry2$cohort=='RML'], cex=1.5)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



ylims = c(0,2750)

par(mar=c(1,5,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
abline(v=0)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=(0:(max(ylims)/100))*100, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=(0:ceiling(max(ylims/500)))*500, labels=formatC((0:ceiling(max(ylims/500)))*500, big.mark=','), lwd=0, lwd.ticks=1, tck=-0.050, las=2)
mtext(side=2, line=3.5, text='plasma NfL (pg/mL)')
mtext(side=3, line=0.5, text='neurofilament')
# axis(side=2, at=(0:(max(ylims)/100))*100, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
# axis(side=2, at=(0:ceiling(max(ylims/500)))*500, labels=formatC((0:ceiling(max(ylims/500)))/2, format='f', digits=1), lwd=0, lwd.ticks=1, tck=-0.050, las=2)
# mtext(side=2, line=3.5, text='plasma NfL (ng/mL)')

# for (ms in unique(nfl$animal)) {
#   subs = subset(nfl, animal==ms)
#   points(subs$timepoint, subs$nfl_pgml, type='l', col=subs$color, lwd=0.25)
# }
for (coh in unique(nfl_smry$cohort)) { # do them in reverse order so pre-treatment is on top at 81 dpi
  subs = subset(nfl_smry, cohort==coh)
  color = subs$color[1]
  points(subs$timepoint, subs$mean_nfl, type='l', lwd=3, col=color, lty=subs$lty)
  points(subs$timepoint, subs$mean_nfl, pch=19, col=color)
  polysubs = subset(subs, !is.na(sd_nfl)) 
  polygon(x=c(polysubs$timepoint,rev(polysubs$timepoint)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
}
par(xpd=T)
text(x=nfl_smry$timepoint[nfl_smry$cohort=='RML'], y=rep(max(ylims), sum(nfl_smry$cohort=='RML')), labels=nfl_smry$p_symb[nfl_smry$cohort=='RML'], cex=1.5)
xleg = 0
yleg = 2250
legend(x=xleg,y=yleg,legend=cohorts$cohort,col=cohorts$color,text.col=cohorts$color,text.font=2,lwd=3,lty=cohorts$lty,bty='n',title='cohorts',title.col='black')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


sf = survfit(Surv(dpi, status) ~ cohort, data=master)
sf$cohort = gsub('cohort=','',names(sf$strata))
sf$color = cohorts$color[match(sf$cohort, cohorts$cohort)]
sf$lty = cohorts$lty[match(sf$cohort, cohorts$cohort)]

ylims = c(0, 1.05)
par(mar=c(4,5,3,1))
plot(sf, ann=FALSE, axes=FALSE, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
abline(v=0)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3.5, text='survival')
mtext(side=1, line=2.0, text='days post-infection')
mtext(side=3, line=0.5, text='mortality')

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

dev.off() #### END FIGURE 3















  
#### FIGURE S6: some different ways of slicing and dicing figure 3
imgsave(paste('display_items/figure-s6.',imgmode,sep=''),width=6.5*resx,height=2.3*resx,res=resx)

layout_matrix = matrix(1:3,nrow=1,byrow=T)
layout(layout_matrix)

expt = 'N'

master = read.table(paste('data/mri/',expt,'_survival.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
cohorts = read.table(paste('data/mri/',expt,'_cohorts.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
rr = read.table(paste('data/mri/',expt,'_rotarod.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
nfl = read.table(paste('data/mri/',expt,'_nfl.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
nests = read.table(paste('data/mri/',expt,'_nests.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
weights = read.table(paste('data/mri/',expt,'_weights.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

linecol = '#AAAAAA'

nfl$color = cohorts$color[match(nfl$cohort, cohorts$cohort)]
nfl$lty = cohorts$lty[match(nfl$cohort, cohorts$cohort)]

rr$cohort = master$cohort[match(rr$animal,master$animal)]
rr$color = cohorts$color[match(rr$cohort, cohorts$cohort)]

rr_smry1 = sqldf("
                 select   cohort, color, animal, dpi, avg(latency) mean_latency
                 from     rr
                 group by 1, 2, 3, 4
                 order by 1, 2, 3, 4
                 ;")

rr_smry2 = sqldf("
                 select   cohort, dpi, avg(mean_latency) mean_lat, stdev(mean_latency) sd_lat, count(*) n_lat
                 from     rr_smry1
                 group by 1, 2
                 having   count(*) > 1
                 order by 1, 2
                 ;")
rr_smry2 = calculate_ci(rr_smry2)
rr_smry2$color = cohorts$color[match(rr_smry2$cohort, cohorts$cohort)]
rr_smry2$lty = cohorts$lty[match(rr_smry2$cohort, cohorts$cohort)]

# add P vals
rr_smry2$ks_p = as.numeric(NA)
for (i in 1:nrow(rr_smry2[rr_smry2$cohort=='RML',])) {
  dpi = rr_smry2$dpi[i]
  rml_rows  = rr_smry1$dpi == dpi & rr_smry1$cohort == 'RML'
  unin_rows = rr_smry1$dpi == dpi & rr_smry1$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(rr_smry1$mean_latency[rml_rows], rr_smry1$mean_latency[unin_rows], alternative='two.sided'))
  rr_smry2$ks_p[i] = ks_obj$p.value
}
rr_smry2$p_symb = ''
rr_smry2$p_symb[!is.na(rr_smry2$ks_p) & rr_smry2$ks_p < 0.05] = '*'
rr_smry2$p_symb[!is.na(rr_smry2$ks_p) & rr_smry2$ks_p < 0.01] = '**'
rr_smry2$p_symb[!is.na(rr_smry2$ks_p) & rr_smry2$ks_p < 0.001] = '***'

weights$cohort = master$cohort[match(weights$animal, master$animal)]
weights$color = cohorts$color[match(weights$cohort, cohorts$cohort)]
weighall_timepoints = data.frame(dpi=unique(weights$dpi[weights$cohort=='uninoculated']))
weights$plot_dpi = weights$dpi
needmatch = !(weights$plot_dpi %in% weighall_timepoints$dpi)
weights$plot_dpi[needmatch] = weighall_timepoints$dpi[findInterval(x=weights$dpi[needmatch], vec=weighall_timepoints$dpi, left.open=T, rightmost.closed=T)+1]
# check that rounding up worked as expected:
# weights[weights$dpi != weights$plot_dpi,] # yes - looks good
weights_smry = sqldf("
                     select   cohort, plot_dpi, avg(wt) mean_weight, stdev(wt) sd_weight, count(*) n_wt
                     from     weights
                     group by 1, 2
                     having   count(*) > 1 -- only include dpi with >1 mouse measured
                     order by 1, 2
                     ;")
weights_smry = calculate_ci(weights_smry)
weights_smry$color = cohorts$color[match(weights_smry$cohort, cohorts$cohort)]
weights_smry$lty =     cohorts$lty[match(weights_smry$cohort, cohorts$cohort)]

# add P vals
weights_smry$ks_p = as.numeric(NA)
for (i in 1:nrow(weights_smry[weights_smry$cohort=='RML',])) {
  dpi = weights_smry$plot_dpi[i]
  rml_rows  = weights$plot_dpi == dpi & weights$cohort == 'RML'
  unin_rows = weights$plot_dpi == dpi & weights$cohort == 'uninoculated'
  ks_obj = suppressWarnings(ks.test(weights$wt[rml_rows], weights$wt[unin_rows], alternative='two.sided'))
  weights_smry$ks_p[i] = ks_obj$p.value
}
weights_smry$p_symb = ''
weights_smry$p_symb[!is.na(weights_smry$ks_p) & weights_smry$ks_p < 0.05] = '*'
weights_smry$p_symb[!is.na(weights_smry$ks_p) & weights_smry$ks_p < 0.01] = '**'
weights_smry$p_symb[!is.na(weights_smry$ks_p) & weights_smry$ks_p < 0.001] = '***'

panel = 1


xlims = c(75,185)
ylims = c(10,40)
par(mar=c(4,4,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=-3:60*10, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5, cex.axis=0.8)
mtext(side=1, line=2.0, text='days post-infection')
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
abline(v=0)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=0:8*5, labels=0:8*5, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='body weight (g)')
mtext(side=3, line=0.5, text='weights')
for (anml in unique(master$animal)) {
  s = subset(weights, animal==anml)
  points(s$dpi, s$wt, col=s$color, type='l', lwd=0.5, lty=s$lty)
  points(s$dpi, s$wt, col=s$color, pch=20, cex=0.75)
}
# for (coh in unique(weights_smry$cohort)) {
#   s = subset(weights_smry, cohort==coh)
#   points(s$plot_dpi, s$mean_weight, col=s$color, type='l', lwd=3, lty=s$lty)
#   points(s$plot_dpi, s$mean_weight, col=s$color, pch=19)
#   polygon(x=c(s$plot_dpi[!is.na(s$sd_weight)],rev(s$plot_dpi[!is.na(s$sd_weight)])), y=c(s$u95[!is.na(s$sd_weight)],rev(s$l95[!is.na(s$sd_weight)])), col=alpha(s$color[!is.na(s$sd_weight)],ci_alpha), border=NA)
# }
par(xpd=T)
text(x=weights_smry$plot_dpi[weights_smry$cohort=='RML'], y=rep(max(ylims), sum(weights_smry$cohort=='RML')), labels=weights_smry$p_symb[weights_smry$cohort=='RML'], cex=0.9)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

xlims = c(-15,190)

ylims = c(0,700)

par(mar=c(4,4,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5, cex.axis=0.8)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
mtext(side=1, line=2.0, text='days post-infection')
abline(v=0)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=(0:(max(ylims)/100))*100, lwd=1, lwd.ticks=1, tck=-0.025, las=2)
mtext(side=2, line=3.5, text='drop latency (s)')
mtext(side=3, line=0.5, text='rotarod')
for (anml in unique(master$animal)) {
  s = subset(rr_smry1, animal==anml)
  points(s$dpi, s$mean_latency, col=s$color, type='l', lwd=0.5, lty=s$lty)
  points(s$dpi, s$mean_latency, col=s$color, pch=20, cex=0.75)
}
# for (coh in unique(rr_smry2$cohort)) { # do them in reverse order so pre-treatment is on top at 81 dpi
#   subs = subset(rr_smry2, cohort==coh)
#   color = subs$color[1]
#   points(subs$dpi, subs$mean_lat, type='l', lwd=3, col=color, lty=subs$lty)
#   points(subs$dpi, subs$mean_lat, pch=19, col=color)
#   polysubs = subset(subs, !is.na(sd_lat)) 
#   polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
# }
par(xpd=T)
text(x=rr_smry2$dpi[rr_smry2$cohort=='RML'], y=rep(max(ylims), sum(rr_smry2$cohort=='RML')), labels=rr_smry2$p_symb[rr_smry2$cohort=='RML'], cex=1.5)
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

#ylims = c(0,2750)
ylims=c(10,10000)
yats = rep(1:9, 4) * rep(10^(1:4), each=9)

par(mar=c(4,4,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F, log='y')
axis(side=1, at=seq(-20,200,10), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(-7,seq(0,180,30)), labels=NA, lwd=0, lwd.ticks=1,  tck=-0.05)
axis(side=1, at=c(-7,seq(0,180,30)), lwd=0, lwd.ticks=0, line=-0.5, cex.axis=0.8)
abline(v=c(-7,seq(0,180,30)), lwd=0.125, col=linecol)
mtext(side=1, line=2.0, text='days post-infection')
abline(v=0)
points(x=0,y=max(ylims),pch=25,col='black',bg='black')
axis(side=2, at=10^(1:4), labels=NA, lwd=1, lwd.ticks=1, tck=-0.05)
axis(side=2, at=10^(1:4), labels=c('10','100','1K','10K'), lwd=0, las=2, line=-0.5)
axis(side=2, at=yats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
#axis(side=2, at=(0:(max(ylims)/100))*100, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
#axis(side=2, at=(0:ceiling(max(ylims/500)))*500, labels=formatC((0:ceiling(max(ylims/500)))*500, big.mark=','), lwd=0, lwd.ticks=1, tck=-0.050, las=2)
mtext(side=2, line=2.5, text='plasma NfL (pg/mL)')
mtext(side=3, line=0.5, text='neurofilament')
for (animal in unique(nfl$animal)) {
  rows = nfl$animal == animal
  points(nfl$timepoint[rows], nfl$nfl_pgml[rows], type='l', lwd=0.5, col=nfl$color[rows], lty=nfl$lty[rows])
  points(nfl$timepoint[rows], nfl$nfl_pgml[rows], pch=20, cex=0.75, col=nfl$color[rows])
}

par(xpd=T)
xleg = 0
yleg = 10000
legend(x=xleg,y=yleg,legend=cohorts$cohort,col=cohorts$color,text.col=cohorts$color,text.font=2,lwd=3,lty=cohorts$lty,bty='n',title='cohorts',title.col='black',cex=0.75)
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

dev.off()




















# more detail on N behaviorals

expt = 'N'
behavs = read.table(paste('data/mri/',expt,'_behavs.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')
behavs_meta = read.table(paste('data/mri/',expt,'_behavs_meta.tsv',sep=''),sep='\t',header=T,quote='',comment.char='')

smry_by_obs = sqldf("
select   m.cohort, b.dpi, b.observation, avg(b.score) mean_score, stdev(b.score) sd_score, count(*) n_animals
from     behavs b, master m
where    b.animal = m.animal
group by 1, 2, 3
order by 1, 2, 3
;")
smry_by_obs = calculate_ci(smry_by_obs)
smry_by_obs$col = cohorts$color[match(smry_by_obs$cohort,cohorts$cohort)]
smry_by_obs$lty = cohorts$lty[match(smry_by_obs$cohort,cohorts$cohort)]

imgsave(paste('display_items/figure-s5.',imgmode,sep=''),width=6.5*resx,height=8.5*resx,res=resx)
par(mfrow=c(8,5))
# create <= 33 character display titles where needed
behavs_meta$observations_display = behavs_meta$observations
to_fix = nchar(behavs_meta$observations) > 33
behavs_meta$observations_display[to_fix] = gsub('\\(.*\\)','',behavs_meta$observations_display[to_fix])
# set params
xlims = c(85,190)
xats = -3:60*10
xbigs = c(90, 120, 150, 180)
par(mar=c(2, 2, 1, 0.5))
for (i in 1:nrow(behavs_meta)) {
  smry_rows = smry_by_obs$observation == behavs_meta$observations[i]
  ylims = c(0, max(max(smry_by_obs$u95[smry_rows], na.rm=T),1)*1.05)
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
  axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=xbigs, labels=NA, lwd=0, lwd.ticks=1, tck=-0.050)
  axis(side=1, at=xbigs, lwd=0, line=-1, cex.axis=0.5)
  abline(v=xbigs, lwd=0.125, col=linecol)
  axis(side=2, at=seq(0, max(ylims)*2, by=1), lwd=1, lwd.ticks=1, las=2, cex.axis=0.5)
  for (coh in cohorts$cohort) {
    plot_rows = smry_rows & smry_by_obs$cohort == coh
    color = smry_by_obs$col[plot_rows][1]
    points(x=smry_by_obs$dpi[plot_rows], y=smry_by_obs$mean_score[plot_rows], type='l', lwd=1, col=smry_by_obs$col[plot_rows], lty=smry_by_obs$lty[plot_rows])
    polysubs = subset(smry_by_obs, plot_rows & !is.na(smry_by_obs$sd_score))
    polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
  }
  mtext(side=3, line=0, text=behavs_meta$observations_display[i], cex=0.5)
}
dev.off()

