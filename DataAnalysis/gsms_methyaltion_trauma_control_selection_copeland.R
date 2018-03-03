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

nd = "C:\\Users\\u6006921\\Dropbox\\trauma_methyl\\data_anly\\"
sd = "C:\\Users\\u6006921\\Dropbox\\trauma_methyl\\data_anly\\trauma_data\\"
setwd(nd)




#########
#########
###        Select stage 2 controls
#########
#########

ev123p = read.table(paste0(nd,'managed_data_for_methyl_trauma_ctl_select.csv'), sep=',', header=T)
stage1 = read.table(paste0(nd,'revised_GSMS_trauma_selections_021916.csv'), sep=',', header=T)

setDT(ev123p)
setDT(stage1)

##Check out age ranges for trauma selection as these dictate t1 and t2 for stage 2 according to grant
mean(stage1[gt1==1,cohage])
## 11.7
mean(stage1[gt2==1,cohage])
## 13.6
mean(stage1[gt3==1,cohage])
## 25.2


#########
###        Calculate new events count using adult data
#########

itm2 = c('jnvio1r','jnvio2r','jnl11r','jnl13r','jnl14r','jnl15r','jnl16r','jnl17r','jnl18r','jnl19r','jnl20r','jnl21r','jnl22r','jnl23r','jnl24r') 
itm1 = c("driskl73","driskl74","driskl77","drisk70","driskl80","driskl81","driskl82","driskl84","driskl85","drisk89","drisk15","schoolbully1","drisk57","drisk72","drisk76","j3los4","j3los5") 


for (i in itm1) {
  eval(parse(text = paste0('ev123p[, ',i,'y := as.numeric(ifelse(cohage>=18 , 0, ',i,')), by=gsmsid ]')))
  eval(parse(text = paste0('ev123p[, ',i,'y := as.numeric(ifelse(cumsum(',i,'y)>0, 1, 0)), by=gsmsid ]')))
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
names(ctl)

for_bill = subset(ctl, select=c("gsmsid","sample_id","cohage","wave","maltreat","m16","eyng","eold","eyng_old","gt1","gt2","race3","jsex",itm1,itm2,gudy))

#write.table(for_bill, file = 'control_data_WORK_WITH_THIS_BILL.csv', sep=',', col.names = T, row.names = F)


#########
#########
###        Beyond this point, I winow the controls down to our curent selections. Bill, you will most likely want to use the file above
###         to identify additional selections. I include the code below in the event you are curious how we got where we are currently.
#########
#########


### Use orig grant criteria

ctl = subset(ctl,gt1==1 & gt2 == 1 ) #, select = c("gsmsid","sample_id","cohage","wave","maltreat","maltreatr","m16","eyng","eold","eyng_old","gt1","gt2" ))
length(as.character(unique(ctl$gsmsid)))
## ctl: sub = 197 ; obs = 1205
names(ctl)



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
###  For the current list of controls, Edwin requested that adult obs be excluded for subs experiencing maltreat or more than 1 
###    event in adulthood (as well as the earlier criteria of no youth events and no youth trauma) 
###     This knocked the # of controls w adult obs to n=180. the # of controls w childhood obs still n=197 
 


c2 = subset(c, maltreat==0 & eold<2)
names(c2)
table(c2$t1,c2$t2)


###
###  Remove samples already used
###


resid = subset(for_bill, !(gsmsid %in% c2$gsmsid))
resid2 = subset(for_bill, !(sample_id %in% c2$sample_id))

write.table(resid, file = 'control_data_WORK_WITH_THIS_BILL_minus_subs_previously_included.csv', sep=',', col.names = T, row.names = F)
write.table(resid2, file = 'control_data_WORK_WITH_THIS_BILL_minus_obs_previously_included.csv', sep=',', col.names = T, row.names = F)
