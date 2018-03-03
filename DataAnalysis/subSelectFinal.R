#########
###             Set dir and data file names
#########

# library(ggplot2)
# library(plyr)
# library(psych)
# library(readxl)
# library(zoo)
library(data.table)
library(sas7bdat)
library(doBy)
library(lme4)

nd = "C:\\Users\\deadkins\\Dropbox\\trauma_methyl\\data_anly\\"
sd = "C:\\Users\\deadkins\\Dropbox\\trauma_methyl\\data_anly\\trauma_data\\"
setwd(nd)


#########
###          Read and manage data 
#########

## My first revision of picks
fin = read.table(paste0(nd,'gsmsNewandOldSpots.csv'), sep=',', header=T)


#########
###        Manage Bill's 1st data file
#########

evt = read.sas7bdat(paste0(sd,'methylgrant.sas7bdat'))
names(evt) = tolower(gsub(" ", "", names(evt)))
ev1 = subset(evt, select =c(driskl73,driskl74,driskl77,drisk70,driskl80,driskl81,driskl82,driskl84,driskl85,drisk89,drisk15,schoolbully1,maltreat,j4nda,jtotanx,events,malby16:group,sample_id,gsmsid,age,cohage,wave,jsex,race3))

table(ev1$driskl73,useNA='ifany')

is.nan.data.frame <- function(ev1)
  do.call(cbind, lapply(ev1, is.nan))
ev1[is.nan(ev1)] <- 0

setDT(ev1)
setkey(ev1,'gsmsid','cohage')
itm = names(ev1)[1:13]

for (i in itm) {
  eval(parse(text = paste0('ev1[, ',i,'r := as.numeric(ifelse(cumsum(',i,')>0, 1, 0)), 
                           by=gsmsid ]')))
}

### Bill's event measure is off due to including drisk70 in the measure twice
ev1$e2 = rowSums(subset(ev1, select=driskl73r:schoolbully1r)) + ev1$drisk70r
ev1$e1 = rowSums(subset(ev1, select=driskl73r:schoolbully1r))
table(ev1$events,ev1$e2, useNA = 'ifany')


#########
###         Manage Bill's 2nd data, merge his 1st and 2nd datafiles
#########

evn2= read.sas7bdat(paste0(sd,'daniel2_1.sas7bdat'))
names(evn2) = tolower(gsub(" ", "", names(evn2)))
evn2 = evn2[order(evn2$gsmsid, evn2$wave),]

is.nan.data.frame <- function(evn2)
  do.call(cbind, lapply(evn2, is.nan))
evn2[is.nan(evn2)] <- 0

evn2 = subset(evn2, select=-c(drisk53,drisk75,drisk90,driskl90,j4ncd3))
 
itm = names(evn2)[3:(length(evn2))]
setDT(evn2)
setkey(evn2,'gsmsid')
for (i in itm) {
  eval(parse(text = paste0('evn2[, ',i,'r := as.numeric(ifelse(cumsum(',i,')>0, 1, 0)), 
                           by=gsmsid ]')))
}
evn2$evn2 = 1


### Merge Bill 1st and Bill 2nd
ev12 = merge(ev1, evn2,  by=c('gsmsid','wave'), all=T)
table(ev12$evn2, useNA = 'ifany')


#########
###         Manage Bill's 3nd data, merge his 1_2 and 3rd datafiles
#########

ev3= read.sas7bdat(paste0(sd,'daniel3.sas7bdat'))
names(ev3) = tolower(gsub(" ", "", names(ev3)))
ev3 = ev3[order(ev3$gsmsid, ev3$wave),]

is.nan.data.frame <- function(ev3)
  do.call(cbind, lapply(ev3, is.nan))
ev3[is.nan(ev3)] <- 0

ev3 = subset(ev3,select=c(gsmsid,wave,j3los4,j3los5))

itm = names(ev3)[3:(length(ev3))]
setDT(ev3)
setkey(ev3,'gsmsid')
for (i in itm) {
  eval(parse(text = paste0('ev3[, ',i,'r := as.numeric(ifelse(cumsum(',i,')>0, 1, 0)), 
                           by=gsmsid ]')))
}

ev3$ev3 = 1
ev123 = merge(ev12, ev3,  by=c('gsmsid','wave'), all=T)
table(ev123$ev3, useNA = 'ifany')

ev123[, dup := .N, by = c('gsmsid','wave','cohage','sample_id')]
table(ev123$dup,useNA='ifany')  ## 3 dups, the 'unique' fx below drops them
ev123 = as.data.frame(ev123)
ev123 = unique(ev123)
setDT(ev123)

ev123 = ev123[, dup := .N, by = c('sample_id')]
ev123 = ev123[, flag := ifelse(sample_id!='' & dup>1,1,0)]
table(ev123$flag, useNA='ifany')   ## 32 obs with dup sample_ids (=16 mismatched ids)
ev123$b4p = 1

ev123[, dup := .N, by = c('gsmsid','wave','cohage')]
table(ev123$dup)   ## 23 obs dup on these 3 vars
ev123[, dup := .N, by = c('gsmsid','wave','age')]
table(ev123$dup)   ## 12 obs dup on these 3 vars
##   Visual inspection suggests some subs gave more than 1 blood samp per wave, sometimes same year (the 2nd table) & sometimes different years (1st minus 2nd table) 
##        (altho prob also data error in files we received)



###
### Merge in Bill's adult event data
###

ev4= read.sas7bdat(paste0(sd,'daniel6.sas7bdat'))
names(ev4) = tolower(gsub(" ", "", names(ev4)))
ev4 = ev4[order(ev4$gsmsid, ev4$wave),]

is.nan.data.frame <- function(ev4)
  do.call(cbind, lapply(ev4, is.nan))
ev4[is.nan(ev4)] <- 0

ev4=subset(ev4,select=c(gsmsid,wave,cohage,jnvio1,jnvio2,jnvio3,jnvio4,jnl7,jnl8,jnl11,jnl13,jnl14,jnl15,jnl16,jnl17,jnl18,jnl19,jnl20,jnl21,jnl22,jnl23,jnl24))

itm = names(ev4)[4:(length(ev4))]
setDT(ev4)
setkey(ev4,'gsmsid')
for (i in itm) {
  
  eval(parse(text = paste0('ev4[, ',i,'y := as.numeric(ifelse(cohage>=18, 0, ',i,')), by=gsmsid ]')))
  eval(parse(text = paste0('ev4[, ',i,'y := as.numeric(ifelse(cumsum(',i,'y)>0, 1, 0)), by=gsmsid ]')))
  eval(parse(text = paste0('ev4[, ',i,'o := as.numeric(ifelse(cohage<18, 0, ',i,')), by=gsmsid ]')))
  eval(parse(text = paste0('ev4[, ',i,'o := as.numeric(ifelse(cumsum(',i,'o)>0, 1, 0)), by=gsmsid ]')))
  eval(parse(text = paste0('ev4[, ',i,'r := as.numeric(ifelse(',i,'o - ',i,'y>0, 1, 0)), by=gsmsid ]')))
  print(i)
  print(eval(parse(text = paste0('table(ev4$',i,'y,ev4$',i,'o,useNA=\'ifany\')'))))
  print(eval(parse(text = paste0('table(ev4$',i,'y,ev4$',i,'r,useNA=\'ifany\')'))))
  print(eval(parse(text = paste0('table(ev4$',i,'o,ev4$',i,'r,useNA=\'ifany\')'))))
}


ev4s=subset(ev4,select=c(gsmsid,wave,jnvio1r,jnvio2r,jnl11r,jnl13r,jnl14r,jnl15r,jnl16r,jnl17r,jnl18r,jnl19r,jnl20r,jnl21r,jnl22r,jnl23r,jnl24r))
#ev4s=subset(ev4,select=c(gsmsid,wave,jnvio1r,jnvio2r,jnl11r,jnl13r,jnl14r,jnl15r,jnl16r,jnl17r,jnl18r,jnl19r,jnl20r,jnl21r,jnl22r,jnl23r,jnl24r,jnvio1,jnvio2,jnvio3,jnvio4,jnl7,jnl8,jnl11,jnl13,jnl14,jnl15,jnl16,jnl17,jnl18,jnl19,jnl20,jnl21,jnl22,jnl23,jnl24))
ev4s$ev4 = 1
ev123 = merge(ev123, ev4s,  by=c('gsmsid','wave'), all=T)
ev123 = subset(ev123, !is.na(cohage))      #### do not have bloodspots for any of these 145 observations, so nothing is lost
#ev123 = subset(ev123, is.na(cohage))     

###
###   Clean and inspect merged bill/paula. Flag potentially problematic samples. recalc 'eventr' vars to include new
###     obs from paula
###

pau = read.table(paste0(nd,'gsms_paula_cuts_clean.csv'), sep=',', header=T, stringsAsFactors = F)
names(pau)[3]='cohage_pau'
setDT(pau,key=c('gsmsid','wave','cohage_pau','sample_id_pau'))
pau = unique(pau)

ev123p = merge(ev123, pau,  by=c('gsmsid','wave'), all=T)


table(ev123p$b4p, useNA = 'ifany')
x = ev123p[is.na(ev123p$b4p),]        ## pau brings 10 brand new obs from final wave of data collection, can use just recalc cum evts
ev123p = subset(ev123p, !(gsmsid=='A01212' & sample_id!=sample_id_pau))   ## drop 2 dups gen by pau merge
x = subset(ev123p, !is.na(sample_id_pau) & sample_id!=sample_id_pau)
table(x$cohage,x$cohage_pau,useNA = 'ifany')    ## 114 obs of which 113 match cohage and add sample_id's to obs which previously had no record of bloodspot


ev123p[, cohage := ifelse(!is.na(sample_id_pau) & sample_id!=sample_id_pau, cohage_pau,cohage)]   ## corrects one obviously miscoded cohage (based on wave, age and cohage_pau)
ev123p[, sample_id := ifelse(!is.na(sample_id_pau) & sample_id!=sample_id_pau, as.character(sample_id_pau),as.character(sample_id))]

ev123p[, cohage := ifelse(is.na(ev123p$b4p), cohage_pau,cohage)]
ev123p[, sample_id := ifelse(is.na(ev123p$b4p), as.character(sample_id_pau),as.character(sample_id))]
ev123p[, cohage := ifelse(cohage==0, age,cohage)]
ev123p[, cohage := ifelse(sample_id=='T00291', 30,cohage)]
table(ev123p$cohage,ev123p$cohage_pau,useNA = 'ifany')

##
## Append 6 bloodspots from Wanda
##

wnd = read.table(paste0(nd,'gsms_wanda_6bloodspots.csv'), sep=',', header=T, stringsAsFactors = F)
wnd$wave = 19
ev123p = merge(ev123p, wnd,  by=c('gsmsid','wave','cohage','sample_id'), all=T)
names(ev123p)


itm = c("driskl73","driskl74","driskl77","drisk70","driskl80","driskl81","driskl82","driskl84","driskl85","drisk89","drisk15","schoolbully1","drisk57","drisk72","drisk76","j3los4","j3los5") 
for (i in itm) {
  eval(parse(text = paste0('ev123p[, ',i,' := ifelse(cohage>=18 & is.na(',i,'),0,',i,')]')))          ## coding may be problematic for controls
  eval(parse(text = paste0('ev123p[, ',i,'r := as.numeric(ifelse(cumsum(',i,')>0, 1, 0)), by=gsmsid ]')))
}

itm = c('jnvio1r','jnvio2r','jnl11r','jnl13r','jnl14r','jnl15r','jnl16r','jnl17r','jnl18r','jnl19r','jnl20r','jnl21r','jnl22r','jnl23r','jnl24r') 
for (i in itm) {
  eval(parse(text = paste0('ev123p[, ',i,' := ifelse(cohage>=18 & is.na(',i,'),0,',i,')]')))          ## coding may be problematic for controls
  eval(parse(text = paste0('ev123p[, ',i,' := as.numeric(ifelse(cumsum(',i,')>0, 1, 0)), by=gsmsid ]')))
}

x = unname(unlist((unique(subset(ev123p,is.na(driskl73r),select=gsmsid)))))
x = subset(ev123p,gsmsid %in% x,select=c(gsmsid,cohage,wave,j3los4,j3los4r,sample_id))  ## confirms loop above carries fwd 1's and, when missing occurs <age 18, carries fwd NA
table(ev123p$jnvio2r,ev123p$cohage,useNA='ifany')

### Based on convo with Edwin and Bill on 10/15/15, final trauma criteria:
gude = c("driskl73r","driskl74r","driskl77r","drisk70r","driskl80r","driskl81r","driskl82r","driskl84r","driskl85r","drisk89r","drisk15r","schoolbully1r","drisk57r","drisk72r","drisk76r","j3los4r","j3los5r") 

### Missign data file for Bill to fill-in
# ddd = c('jnvio1','jnvio2','jnvio3','jnvio4','jnl7','jnl8','jnl11','jnl13','jnl14','jnl15','jnl16','jnl17','jnl18','jnl19','jnl20','jnl21','jnl22','jnl23','jnl24')
# itm = c("driskl73","driskl74","driskl77","drisk70","driskl80","driskl81","driskl82","driskl84","driskl85","drisk89","drisk15","schoolbully1","drisk57","drisk72","drisk76","j3los4","j3los5") 
# ev123p$efin = rowSums(subset(ev123p, select=c(ddd,itm)))
# xxx = subset(ev123p,is.na(efin), select =c('gsmsid','cohage','wave','sample_id','efin',ddd,itm))
# write.table(xxx, file = 'gsms_trauma_missing_data.csv', sep=',', col.names = T, row.names = F)


#########
###         New trauma selections
#########


ev123p$efin = rowSums(subset(ev123p, select=gude))
ev123p$eorg = rowSums(subset(ev123p, select=driskl73r:schoolbully1r)) + ev123p$drisk70r
ev123p[, evmal := as.integer(ifelse(max(maltreat)==1,1L,0L)), by=gsmsid]
ev123p[, maltreatr := as.numeric(ifelse(cumsum(maltreat)>0, 1, 0)), by=gsmsid ]
ev123p[, m16 := as.integer(ifelse(maltreat==1 &  cohage<=16,1L,0L))]
ev123p[, m16id := as.integer(ifelse(max(m16)==1,1L,0L)), by=gsmsid]
ev123p[, group := mean(group,na.rm=T), by=gsmsid]
ev123p = ev123p[with(ev123p, order(gsmsid, cohage)), ]

### Examine breakdown of trauma counts by definition (bill org, my fin)
summary(subset(ev123p, select=c("efin","eorg")))
table(ev123p$efin,ev123p$eorg, useNA = 'ifany')
table(ev123p$efin, useNA = 'ifany')
table(ev123p$eorg, useNA = 'ifany')
## 53 missing efin, but present in eorg


x = subset(ev123p,is.na(maltreat))
x = subset(ev123p,is.na(m16))
table(ev123p$maltreat,ev123p$group,useNA='ifany')
table(ev123p$wave,ev123p$cohage,useNA='ifany')

#########
###         Test Bill's trauma selections numbers
#########

#Group - "This is my <ie,  Bill's> variable for group status for this study. 3 = maltreatment reported by age 16 (336 individual subjects); 2 = no maltreatment but traumatic events between initial assessment and last childhood assessment (270 individual subjects); 1=no report of traumatic events or maltreatment by age 16 (353 individual subjects); 0=all other subjects (461 individual subjects)."
table(ev123p$group, useNA = 'ifany')
length(unique(ev123p[ev123p$group==2 & sample_id!='',gsmsid]))
length(unique(ev123p[ev123p$group==2 ,gsmsid]))
#264. if i don't require bloodspot present = 270 (= number bill gives in 9/3/15 email)

x = ev123p[ev123p$group==2,]
x[, gt1 := as.integer(ifelse(sample_id!='' & events==0 &  cohage<=16,1L,0L))]
x[, gud1 := as.integer(ifelse(max(gt1)==1,1L,0L)), by=gsmsid]
x = subset(x, gud1==1)
x[, gt2 := as.integer(ifelse(sample_id!='' & events>0 & cohage<=16,1L,0L))]
x[, gud2 := as.integer(ifelse(max(gt2)==1,1L,0L)), by=gsmsid]
x = subset(x, gud2==1)
x[, gt3 := as.integer(ifelse(sample_id!='' & events>0 & cohage>=18,1L,0L))]
x[, gud3 := as.integer(ifelse(max(gt3)==1,1L,0L)), by=gsmsid]
x = subset(x, gud3==1)
length(unique(x$gsmsid))
## starting with bill's group==2, i get 180 trauma subs


table(ev123p$group,ev123p$m16id, useNA='ifany')
x = subset(ev123p, m16id==0)   ## confirmed that m16id = malby16
x[, gt1 := as.integer(ifelse(sample_id!='' & events==0 &  cohage<=16,1L,0L))]
x[, gud1 := as.integer(ifelse(max(gt1)==1,1L,0L)), by=gsmsid]
x = subset(x, gud1==1)
x[, gt2 := as.integer(ifelse(sample_id!='' & events>0 & cohage<=16,1L,0L))]
x[, gud2 := as.integer(ifelse(max(gt2)==1,1L,0L)), by=gsmsid]
x = subset(x, gud2==1)
x[, gt3 := as.integer(ifelse(sample_id!='' & events>0 & cohage>=18,1L,0L))]
x[, gud3 := as.integer(ifelse(max(gt3)==1,1L,0L)), by=gsmsid]
x = subset(x, gud3==1)
length(unique(x$gsmsid))
## starting with bill's malby16==1, i get 181 trauma subs 

table(ev123p$events,ev123p$eorg,useNA='ifany')  
x = subset(ev123p, events!=eorg)  ## the 30 values where events==0 and eorg>0 all had cohage recoded from 0 to age above in Paula-merge section


#########
###        Exclude maltreat if occurs by 16 ('m16id')
#########

x = subset(ev123p, m16id==0)
x[, gt1 := as.integer(ifelse(sample_id!='' & efin==0 &  cohage<=16,1L,0L))]
x[, gud1 := as.integer(ifelse(max(gt1)==1,1L,0L)), by=gsmsid]
x = subset(x, gud1==1)
x[, gt2 := as.integer(ifelse(sample_id!='' & efin>0 & cohage<=16,1L,0L))]
x[, gud2 := as.integer(ifelse(max(gt2)==1,1L,0L)), by=gsmsid]
x = subset(x, gud2==1)
x[, gt3 := as.integer(ifelse(sample_id!='' & efin>0 & cohage>=18,1L,0L))]
x[, gud3 := as.integer(ifelse(max(gt3)==1,1L,0L)), by=gsmsid]
x = subset(x, gud3==1)
length(unique(x$gsmsid))
## N=192 using final event definition 
O00289 = subset(x,sample_id=='O00289', select=c(gsmsid,sample_id,cohage,wave,efin,m16,gt1,gt2,gt3,flag,jsex,race3))

y = subset(ev123p, m16id==0)
y[, gt1 := as.integer(ifelse(sample_id!='' & efin==0 &  cohage<=16,1L,0L))]
y[, gud1 := as.integer(ifelse(max(gt1)==1,1L,0L)), by=gsmsid]
y = subset(y, gud1==1)
y[, gt2 := as.integer(ifelse(sample_id!='' & efin>0 & cohage<=16,1L,0L))]
y[, gud2 := as.integer(ifelse(max(gt2)==1,1L,0L)), by=gsmsid]
y = subset(y, gud2==1)
length(unique(y$gsmsid))
## N=207 if exclude t3 criteria

xid = (unique(as.character(x$gsmsid)))
yid = (unique(as.character(y$gsmsid)))
length(xid[xid %in% yid])   # All 192 subs from 3 timepoint criteria are included in the 207 subs meeting crit t1 and t2 (as expected)

tcs = subset(x, sample_id!='',select=c(gsmsid,sample_id,cohage,wave,efin,m16,gt1,gt2,gt3,flag,jsex,race3))
tid = (unique(as.character(tcs$gsmsid)))
tsid = (unique(as.character(tcs$sample_id)))

#########
###         Choose optimal t1, t2, t3 bloodspots
#########


## t1
tm1 = subset(tcs, gt1==1)
tm1[, op1 := as.integer(ifelse(cohage==max(cohage),1L,0L)), by = gsmsid]
tm1 = subset(tm1, op1==1, select = -op1)

## t2
tm2 = subset(tcs, gt2==1)
tm2[, op2 := as.integer(ifelse(efin==max(efin),1L,0L)), by = gsmsid]
tm2 = subset(tm2, op2==1)
tm2[, op2 := as.integer(ifelse(cohage==min(cohage),1L,0L)), by = gsmsid]
tm2 = subset(tm2, op2==1, select = -op2)
tm2 = unique(tm2)  ## drop one dup row

## t3
tm3 = subset(tcs, gt3==1)
tm3[, op3 := as.integer(ifelse(abs(cohage-28.1)==min(abs(cohage-28.1)),1L,0L)), by=gsmsid]
tm3 = subset(tm3, op3==1, select = -op3)


## Combine
sel = rbind(tm1,tm2,tm3)
sel = sel[order(sel$gsmsid, sel$cohage),]
table(sel$cohage[sel$gt3==1], useNA = 'ifany')
table(sel$flag, useNA = 'ifany')
table(tcs$flag, useNA = 'ifany')

sid = (unique(as.character(sel$gsmsid)))
ssid = (unique(as.character(sel$sample_id)))


#########
###         Read current stage 1 sample selection list and compare picks with new picks considering paulas new bloodspots 
#########

revl = read.table(paste0(nd,'revised_GSMS_trauma_selections_011516.csv'), sep=',', header=T)
names(revl)
rvl = subset(revl, include_in_stage1_01152016==1)
rid = (unique(as.character(rvl$gsmsid)))
rsid = (unique(as.character(rvl$sample_id)))

length(rid[rid %in% sid])      ### all subs in current list are also in new list incl paulas data
length(rsid[rsid %in% tsid])   ### all 468 bloodspots in current list are also in new list incl paulas data 
length(rsid[rsid %in% ssid])   ### 392 of 468 bloodspots in current list are optimal in new list incl paulas data 

###
###   Confirm that old spots are ok
###

names(rvl)
rvl$rvl = 1
rvl0 = subset(rvl, select=c(gsmsid,sample_id,cohage,wave,sex,race,rvl))
sel$sel = 1
sel0 = subset(sel, select=c(gsmsid,sample_id,sel))

x = merge(tcs, rvl0,  by=c('gsmsid','sample_id'), all=T)
x = merge(x, sel0,  by=c('gsmsid','sample_id'), all=T)
x = as.data.frame(x)
x = unique(x)
setDT(x)

x = subset(x, sel==1 | rvl==1)
x = x[order(x$gsmsid, x$cohage.x),]
table(x$flag,useNA='ifany')             ##### ???
x[, dup := .N, by = c('gsmsid')]
table(x$dup,useNA='ifany')
x = subset(x,dup==4 & (is.na(sel) |is.na(rvl)))   ## differences in opt picks are sensible -- always stem from new obs from pau


###
###   Make new Stage 1 datafile with additional samples identified from paula
###

nsel = subset(sel, gsmsid %in% sid[!(sid %in% rid)], select =-sel)
table(nsel$flag,useNA='ifany')       ## No flags


nsel[, race := as.character(ifelse(max(race3,na.rm=T)==2,'Cherokee','Anglo')), by=gsmsid]
nsel[, sex := as.character(ifelse(max(as.numeric(jsex),na.rm=T)==3,'M','F')), by=gsmsid]
nsel[, race0 := as.numeric(ifelse(max(race3,na.rm=T)==2,2,1)), by=gsmsid]
nsel[, sex0 := as.numeric(ifelse(max(as.numeric(jsex),na.rm=T)==3,1,2)), by=gsmsid]

nsel[, id := as.numeric(seq_len(.N)), by = c('race0','sex0')]
nsel[, id := id/max(id), by = c('race0','sex0')]
nsel[, id2 := min(id), by = gsmsid]
nsel = nsel[order(nsel$id2, nsel$gsmsid, nsel$cohage),]
nsel$prep_order_new = seq(1:length(nsel$id)) + max(revl$prep_order_new)

nsel[, c('flag','jsex','race3','m16','race0','sex0','id','id2') := NULL]
nsel$include_in_stage1_01152016 = 1
nsel$prep_order = NA
nsel$punched = NA
nsel$prep_order = NA
nsel$in1stDNAprep = NA
nsel$old = NA

setDT(revl)
revl[, sample_id := as.character(sample_id)]
revl[, sample_id := ifelse(prep_order==510,'O00289',sample_id)]  ### Next 3 lines replace Q00107 (DO NOT OWN) with O00289
revl[, wave := ifelse(prep_order==510,15,wave)]
revl[, cohage := ifelse(prep_order==510,26,cohage)]
x = rbind(revl,nsel)
#write.table(x, file = 'revised_GSMS_trauma_selections_021916.csv', sep=',', col.names = T, row.names = F)


###
### Calculate number of subjects for which events increase v stay the same betwen t2 and t3 selections
###

xx = subset(x,gt1!=1 &include_in_stage1_01152016==1)
xx[, chg := max(efin) - min(efin), by=gsmsid]
dim(subset(xx,chg!=0))[1]/2










#########
#########
###        Maltreatment selections and controls
#########
#########


##Check out age ranges for trauma selection as these dictate t1 and t2 for stage 2 according to grant
mean(x[gt1==1,cohage])
## 11.7
mean(x[gt2==1,cohage])
## 13.6
mean(x[gt3==1,cohage])
## 25.2


#########
###        Calculate new events count using adult data
#########

itm2 = c('jnvio1r','jnvio2r','jnl11r','jnl13r','jnl14r','jnl15r','jnl16r','jnl17r','jnl18r','jnl19r','jnl20r','jnl21r','jnl22r','jnl23r','jnl24r') 
itm1 = c("driskl73","driskl74","driskl77","drisk70","driskl80","driskl81","driskl82","driskl84","driskl85","drisk89","drisk15","schoolbully1","drisk57","drisk72","drisk76","j3los4","j3los5") 

yyy = subset(ev123p,cohage>=18)
for (i in itm1) {
  eval(parse(text = paste0('ev123p[, ',i,'y := as.numeric(ifelse(cohage>=18 , 0, ',i,')), by=gsmsid ]')))
  eval(parse(text = paste0('ev123p[, ',i,'y := as.numeric(ifelse(cumsum(',i,'y)>0, 1, 0)), by=gsmsid ]')))
#   print(i)
#   eval(parse(text = paste0('print(table(yyy$',i,'))')))
}

gudy = c("driskl73y","driskl74y","driskl77y","drisk70y","driskl80y","driskl81y","driskl82y","driskl84y","driskl85y","drisk89y",
        "drisk15y","schoolbully1y","drisk57y","drisk72y","drisk76y","j3los4y","j3los5y") 

ev123p$eyng = rowSums(subset(ev123p, select=gudy))
ev123p$eold = rowSums(subset(ev123p, select=itm2))
ev123p$eyng_old = rowSums(subset(ev123p, select=c(itm2,gudy)))
# ev123p[, c(gudy,itm2,itm1) := NULL]

table(ev123p$cohage,ev123p$eold, useNA='ifany')
table(ev123p$cohage,ev123p$eyng, useNA='ifany')
table(ev123p$cohage,ev123p$eyng_old, useNA='ifany')

table(ev123p$eyng,ev123p$efin, useNA='ifany')
table(ev123p$eold,ev123p$efin, useNA='ifany')
table(ev123p$eyng_old,ev123p$efin, useNA='ifany')
cor(ev123p$efin,ev123p$eyng_old, use = "complete.obs" )






#########
###        Maltreatment selection
#########

table(ev123p$group)
mal = subset(ev123p,m16id==1 & sample_id!='' & maltreatr ==1)
mal[, gt2 := as.integer(ifelse(max(cohage)>=18,1L,0L)), by=gsmsid]
mal[, gt1 := as.integer(ifelse(min(cohage)<18 ,1L,0L)), by=gsmsid]
mal = subset(mal,gt1==1 & gt2 == 1 )
length(unique(mal$gsmsid))
## N = 203

err = subset(revl, include_in_stage1_01152016==0 & gt1!=1)
setDT(err,key=c('gsmsid','cohage'))
# N = 49 maltreatment cases present in 'err' trauma selections

eid = (unique(as.character(err$gsmsid)))
esid = (unique(as.character(err$sample_id)))
es1id = (unique(as.character(err[gt2==1,sample_id])))
es2id = (unique(as.character(err[gt3==1,sample_id])))

malid = unique(as.character(mal$gsmsid))
malsid = unique(as.character(mal$sample_id))

length(eid[eid %in% malid])  
length(es1id[es1id %in% malsid])   
length(es2id[es2id %in% malsid])   

et1 = eid[es1id %in% malsid]  
et2 = eid[es2id %in% malsid]  
est1 = es1id[es1id %in% malsid]  
est2 = es2id[es2id %in% malsid]  


###
###  Choose opt t1 & t2
###

## t1
mal[, t1 := as.integer(ifelse(abs(cohage-13.6)==min(abs(cohage-13.6)),1L,0L)), by=gsmsid]
mal[, t1 := as.integer(ifelse(gsmsid %in% et1,0L,t1))]
mal[, t1 := as.integer(ifelse(sample_id %in% est1,1L,t1))]
mt2 = subset(mal, t1==1)
mt2$t2 = 0
mean(mt2[,cohage])
## 13.9   (good match to 13.6 mean in t2 of trauma selections)


## t2
mal[, t2 := as.integer(ifelse(abs(cohage-25.2)==min(abs(cohage-25.2)),1L,0L)), by=gsmsid]
mal[, t2 := as.integer(ifelse(gsmsid %in% et2,0L,t2))]
mal[, t2 := as.integer(ifelse(sample_id %in% est2,1L,t2))]

### Patch 4/11/16: missing original selected bloodspot C00200-N00100, replace with C00200-90364
mal[,t2 := ifelse(sample_id=='N00100',0,t2)]
mal[,t2 := ifelse(sample_id=='90364',1,t2)]

mt3 = subset(mal, t2==1)
mt3 = subset(mt3, sample_id!='P00067')  ## drop the 1 instances of mult sample_ids per assessment
mean(mt3[,cohage])
## 24.2  (acceptable match to 25.2 mean in t2 of trauma selections)

m = rbind(mt2,mt3)
msid = (unique(as.character(m$sample_id)))
length(esid[esid %in% msid])      ###  86 (of 98 t2+t3, or 147 total) mal samps in error list are also in current list


###
###  Polish, randomize and output list
###

m[, race := as.character(ifelse(max(race3,na.rm=T)==2,'Cherokee','Anglo')), by=gsmsid]
m[, sex := as.character(ifelse(max(as.numeric(jsex),na.rm=T)==3,'M','F')), by=gsmsid]
m[, race0 := as.numeric(ifelse(max(race3,na.rm=T)==2,2,1)), by=gsmsid]
m[, sex0 := as.numeric(ifelse(max(as.numeric(jsex),na.rm=T)==3,1,2)), by=gsmsid]

m[, id := as.numeric(seq_len(.N)), by = c('race0','sex0')]
m[, id := id/max(id), by = c('race0','sex0')]
m[, id2 := min(id), by = gsmsid]
m = m[order(m$id2, m$gsmsid, m$cohage),]
m$prep_order_maltreat = seq(1:length(m$id))
m[, prep_stage1 := as.numeric(ifelse(sample_id %in% esid,1L,0L))]

m = subset(m, select=c("gsmsid","sample_id","cohage","wave","maltreat","m16","eyng","eold","eyng_old","t1","t2","race","sex","prep_order_maltreat","prep_stage1"))
write.table(m, file = 'GSMS_maltreatment_selections_041116.csv', sep=',', col.names = T, row.names = F)

table(mt2$maltreat)

table(mt3$eyng)
mean(mt3$eyng,na.rm=T)

table(mt3$eold)
mean(mt3$eold,na.rm=T)





#########
###        Control selection
#########

ctl = subset(ev123p, eyng==0 & m16id== 0 & sample_id!='' )
length(as.character(unique(ctl$gsmsid)))
## ctl: sub = 543 ; obs = 1941

table(ctl$cohage, useNA = 'ifany')
table(ctl$cohage,ctl$eold, useNA = 'ifany')
table(ctl$cohage,ctl$maltreatr, useNA = 'ifany')
x = subset(ctl,is.na(maltreatr), select=sample_id)
setDT(ctl, key=c('gsmsid','cohage'))

ctl[, gt2 := as.integer(ifelse(max(cohage)>=18,1L,0L)), by=gsmsid]
ctl[, gt1 := as.integer(ifelse(min(cohage)<18 ,1L,0L)), by=gsmsid]
ctl = subset(ctl,gt1==1 & gt2 == 1 ) #, select = c("gsmsid","sample_id","cohage","wave","maltreat","maltreatr","m16","eyng","eold","eyng_old","gt1","gt2" ))
length(as.character(unique(ctl$gsmsid)))
## ctl: sub = 197 ; obs = 1205


### USe orig grant criteria
ctl1 = subset(ctl,cohage<18)
ctl1[, t1 := as.integer(ifelse(abs(cohage-13.6)==min(abs(cohage-13.6)),1L,0L)), by=gsmsid]
ctl1 = subset(ctl1, t1==1)
mean(ctl1[,cohage])  ### mean = 13.8, very close to t2 stage 1 mean = 13.6
ctl1$t2=0

ctl2 = subset(ctl,cohage>=18)

table(ctl2$cohage,ctl2$eold, useNA = 'ifany')
table(ctl2$cohage,ctl2$maltreatr, useNA = 'ifany')

ctl2[, maltreat := as.numeric(ifelse(is.na(maltreatr),0,maltreat))]
ctl2[, t2 := as.integer(ifelse(maltreat==min(maltreat),1L,0L)), by=gsmsid]
ctl2 = subset(ctl2,t2==1)
ctl2[, t2 := ifelse(eold==min(eold),1L,0L), by=gsmsid]
ctl2 = subset(ctl2,t2==1)
ctl2[, t2 := as.integer(ifelse(abs(cohage-25.2)==min(abs(cohage-25.2)),1L,0L)), by=gsmsid]
ctl2 = subset(ctl2, t2==1)
mean(ctl2[,cohage])  ### mean = 23.3, acceptably close to t2 stage 1 mean = 25.2


table(ctl2$cohage,ctl2$eold, useNA = 'ifany')
table(ctl2$cohage,ctl2$maltreatr, useNA = 'ifany')
table(ctl2$eold,ctl2$maltreatr, useNA = 'ifany')

ctl2$t1=0
c = rbind(ctl1,ctl2)



###
###  Polish, randomize and output list
###

c[, race := as.character(ifelse(max(race3,na.rm=T)==2,'Cherokee','Anglo')), by=gsmsid]
c[, sex := as.character(ifelse(max(as.numeric(jsex),na.rm=T)==3,'M','F')), by=gsmsid]
c[, race0 := as.numeric(ifelse(max(race3,na.rm=T)==2,2,1)), by=gsmsid]
c[, sex0 := as.numeric(ifelse(max(as.numeric(jsex),na.rm=T)==3,1,2)), by=gsmsid]

c[, id := as.numeric(seq_len(.N)), by = c('race0','sex0')]
c[, id := id/max(id), by = c('race0','sex0')]
c[, id2 := min(id), by = gsmsid]
c = c[order(c$id2, c$gsmsid, c$cohage),]
c$prep_order_maltreat = seq(1:length(c$id))
c[, prep_stage1 := as.numeric(ifelse(sample_id %in% esid,1L,0L))]

c = subset(c, select=c("gsmsid","sample_id","cohage","wave","maltreat","m16","eyng","eold","eyng_old","t1","t2","race","sex","prep_order_maltreat","prep_stage1"))
write.table(c, file = 'GSMS_control_selections_021516.csv', sep=',', col.names = T, row.names = F)


###
###  Curate stage 2 mal list to drop evt>2 and mal ==1. Combine stage 2 ctl and mal lists.
###



ctl = read.table(paste0(nd,'GSMS_control_selections_021516.csv'), sep=',', header=T)
ctl = subset(ctl, maltreat==0 & eold<2)
mal = read.table(paste0(nd,'GSMS_maltreatment_selections_021216.csv'), sep=',', header=T)

setDT(ctl)
setDT(mal)

ctl$case = 0
ctl[, id := as.numeric(seq_len(.N)), by = c('race','sex')]
ctl[, id := id/max(id), by = c('race','sex')]
ctl[, id2 := min(id), by = gsmsid]
ctl = ctl[order(ctl$id2, ctl$gsmsid, ctl$cohage),]

mal$case = 1
mal[, id := as.numeric(seq_len(.N)), by = c('race','sex')]
mal[, id := id/max(id), by = c('race','sex')]
mal[, id2 := min(id), by = gsmsid]
mal = mal[order(mal$id2, mal$gsmsid, mal$cohage),]

m_c = rbind(ctl,mal)
m_c = m_c[order(m_c$id2, m_c$gsmsid, m_c$cohage),]
m_c$prep_order_stage2 = seq(1:length(m_c$id))


m_c = subset(m_c, select=c("gsmsid","sample_id","cohage","wave","maltreat","eyng","eold","eyng_old","t1","t2","race","sex","case","prep_stage1","prep_order_stage2"))
write.table(m_c, file = 'GSMS_stage2_selections_032016.csv', sep=',', col.names = T, row.names = F)



###
###   Confirm no overlappping samples in stage 1 and stage 2
###


a = read.table(paste0(nd,'revised_GSMS_trauma_selections_020516.csv'), sep=',', header=T)
b = read.table(paste0(nd,'GSMS_stage2_selections_032016.csv'), sep=',', header=T)
c = subset(a, include_in_stage1_01152016==1)

x = as.character(a$sample_id)
y = as.character(b$sample_id)
z = as.character(c$sample_id)
x %in% y
z %in% y
