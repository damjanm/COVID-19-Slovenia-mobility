


#####packages needed
Sys.setenv(MAKEFLAGS = "-j8")


spline.df=4 #3,4,5
only.spline=FALSE



spline.type ="ns"
degree=1
 # spline.type =c("bs","ns","const") 
 #degree=c(1,2,3)  

mindate=NULL
maxdate=NULL

 

##this combination uses all data (2nd wave from 1.7., first wave up to 30.6)
type=1
datechange1="2020-06-30"  
datechange2="2020-06-01"

N2 =184


load("data.Rdata")

library(smooth)
library(splines)
library(gsheet)
library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(optparse)
library(zoo)




######parameters

 





#### OPTIMAL PARAMETERS:







tau1<-0.2  
tau2<-0.01  
tau3<-0.2
tau5<-0.2
tau7<-0.1 
tau6<-1
tau4<-1

mu1<-6 
mu2<-17  
mu3<-7  
mu4<-14  
mu5<-4
mu6<-11
mu7<-7  



nms<-paste("t1",tau1,"t2",tau2,"t3",tau3,"t4",tau4,"t5",tau5,"t6",tau6,"t7",tau7,"m1",mu1,"m2",mu2,"m3",mu3,"m4",mu4,"m5",mu5,"m6",mu6,"m7",mu7,sep="")

 
# Serial interval.  
 
 mu.serial <- 6.5
cv.serial <- 0.62



##code is later


####use either forecast, or N2, set the other to NULL!
 #it is overriden later by forecats=N-N2, where N is the number of data for Slo
 

d1_pop = 2078890 #population size of Slo, change if JiJs forcasts will someday come true.


# various distributions required for modeling
mean0 = 5.1
cv0 = 0.86; # infection to onset


x1 = rgammaAlt(1e7,mean0,cv0) # infection-to-onset distribution
cv1 =0.25
cv2 = 0.45  

cv3 = 0.45  

cv4 = 0.45  

cv5 = 0.45 

cv6=0.45
cv7=0.45

cvs<-paste("m0",mean0,"cv1",cv1,"cv2",cv2,"cv3",cv3,"cv4",cv4,"cv5",cv5,"cv6",cv6,"cv7",cv7,sep="")

x2 = rgammaAlt(1e7,mu1,cv1) 

pi1<-x1+x2


xd = rgammaAlt(1e7,mu2,cv2)

pi2=x1+xd

xh=rgammaAlt(1e7,mu3,cv3)

pi3 = x1+xh  
pi4 = rgammaAlt(1e7,mu4,cv4)  
pi5 = rgammaAlt(1e7,mu5,cv5)  
pi6 = rgammaAlt(1e7,mu6,cv6)  
pi7 = rgammaAlt(1e7,mu7,cv7)  


ecdf.saved.1 = ecdf(pi1)   
ecdf.saved.2 = ecdf(pi2) 
ecdf.saved.3 = ecdf(pi3) 
ecdf.saved.4 = ecdf(pi4) 
ecdf.saved.5 = ecdf(pi5) 
ecdf.saved.6 = ecdf(pi6) 
ecdf.saved.7 = ecdf(pi7)  



####data

idrow<-1:(nrow(dd)) #-which.max((diff(cumsum(is.na(rev(dd$datum)))))==0)
 
dd <- dd[idrow,]
 



hospitalizirani<-as.numeric(dd$hospitalizirani..trenutno.[-1]) #cumulative! prispevek je lahko negativen!
hospitalizirani[is.na(hospitalizirani)]<-0

hospitalizirani<-c(rep(0,50),hospitalizirani)



hospitaliziraniin<-dd$novi[-1]

hospitaliziraniin[is.na(hospitaliziraniin)]<-0

hospitaliziraniin<-c(rep(0,50),hospitaliziraniin)


hospitaliziraniout<-as.numeric(dd$iz.bol..oskrbe..vsi.[-1])
hospitaliziraniout[is.na(hospitaliziraniout)]<-0

hospitaliziraniout<-c(rep(0,50),hospitaliziraniout)
hospitaliziraniouti<-hospitaliziraniout
for (i in 2:length(hospitaliziraniout)){
  hospitaliziraniouti[i]<-hospitaliziraniout[i]-hospitaliziraniout[i-1]
}
hospitaliziraniout<-hospitaliziraniouti

dd$potrjeni.danes[-1] <- gsub('\\.','',dd$potrjeni.danes[-1])
cases<- c(rep(0,50),as.numeric(dd$potrjeni.danes[-1]))
cases[is.na(cases)]<-0


deaths<- as.numeric(dd$umrli..vsi.[-1]) #cumulative!
deaths[is.na(deaths)]<-0

deaths<-c(rep(0,50),deaths)
deathsi<-deaths
for (i in 2:length(deaths)){
  deathsi[i]<-deaths[i]-deaths[i-1]
}
deaths<-deathsi


icu<-c(rep(0,50), as.numeric(dd$intenzivna.enota..trenutno.[-1])) #cumulative! prispevek je lahko negativen!
icu[is.na(icu)] <- 0




icuin<-c(rep(0,50), as.numeric(dd$I.i[-1])) #cumulative! prispevek je lahko negativen!
icuin[is.na(icuin)] <- 0

icuout<-c(rep(0,50), as.numeric(dd$I.o[-1])) #cumulative! prispevek je lahko negativen!
icuout[is.na(icuout)] <- 0



deathsc<-c(rep(0,50), as.numeric(dd$D.u[-1])) #cumulative! prispevek je lahko negativen!
deathsc[is.na(deathsc)] <- 0

deathsh<-c(rep(0,50), as.numeric(dd$H.D[-1])) #cumulative! prispevek je lahko negativen!
deathsh[is.na(deathsh)] <- 0


day<-strsplit(dd[-1,2],split="\\.")

ff<-function(i,x){
  decembr <- which(unlist(lapply(x, function(iter) iter[2]=='12')))
  if(identical(integer(0), decembr)){
    leto <- '.2020'
  } else{
    if(max(december)<i){
      leto <- '.2021'
    } else{
      leto <- '.2020'
    }
  }
  paste(x[[i]][1],".",strsplit(day[[i]][2]," ")[[1]][2],leto,sep="")
  
}

day<-unlist(lapply(1:length(day),ff,day))
day<-as.Date(day,format="%d.%m.%Y")

day<-c((day[1]-1:50)[order(day[1]-1:50)],day)

contacts<-rep(-0.5,length(day))
contacts[day>=as.Date("2020-10-09")]<-0.5

contacts[day<=as.Date(datechange1)]<-0.5 #to da ok sliko

#contacts[day<as.Date("2020-03-04")]<-0.5 
#contacts[day>as.Date("2020-03-12")]<-0.5 

#####mobilty data


ddm<-read.csv("data.csv")



ddms<-ddm[ddm$sub_region_1=="",]
ddms$date<-as.Date(as.character(ddms$date))

#replace NAs with last value
mynareplace<-function(x){
  n<-length(x)
  for (i in 1:n){
    if (is.na(x[i])) x[i]<-x[i-1]
  }
  x
}

ddms$retail_and_recreation_percent_change_from_baseline<-mynareplace(ddms$retail_and_recreation_percent_change_from_baseline)
ddms$grocery_and_pharmacy_percent_change_from_baseline<-mynareplace(ddms$grocery_and_pharmacy_percent_change_from_baseline)
ddms$parks_percent_change_from_baseline<-mynareplace(ddms$parks_percent_change_from_baseline)

ddms$transit_stations_percent_change_from_baseline<-mynareplace(ddms$transit_stations_percent_change_from_baseline)
ddms$workplaces_percent_change_from_baseline<-mynareplace(ddms$workplaces_percent_change_from_baseline)
ddms$residential_percent_change_from_baseline<-mynareplace(ddms$residential_percent_change_from_baseline)

 




ddms$grocery_and_pharmacy_percent_change_from_baseline<-apply(cbind(
  ddms$grocery_and_pharmacy_percent_change_from_baseline,
  ddms$retail_and_recreation_percent_change_from_baseline,
  ddms$transit_stations_percent_change_from_baseline
),1,mean)

fg<-rollmean(ddms[,9:14],k=10,fill=NA,align="left")
fg<-apply(fg,2,function(x) {x[is.na(x)]<-x[which(is.na(x)==TRUE)[1]-1];x})

ddms[,9:14]<-fg
#ddms$grocery_and_pharmacy_percent_change_from_baseline<- ddms$grocery_and_pharmacy_percent_change_from_baseline
#ddms$parks_percent_change_from_baseline<-apply(cbind(ddms$parks_percent_change_from_baseline,ddms$residential_percent_change_from_baseline),1,mean)

#ddms$parks_percent_change_from_baseline<-ddms$residential_percent_change_from_baseline

ddms$scools<-rep(-0.5,nrow(ddms))
ddms$scools[which(ddms$date>="2020-03-16"&ddms$date<="2020-05-18")]<-0.5 
#ddms$scools[which(ddms$date>="2020-06-25"&ddms$date<="2020-09-1")]<-1 
ddms$scools[which(ddms$date>="2020-10-19")]<-0.5 


 



datemax<-which(day==ddms$date[nrow(ddms)])
datemin<-which(day==ddms$date[1])
contacts<-contacts[datemin:datemax]
hospitalizirani<-hospitalizirani[datemin:datemax]
hospitaliziraniin<-hospitaliziraniin[datemin:datemax]
hospitaliziraniout<-hospitaliziraniout[datemin:datemax]
cases<-cases[datemin:datemax]
deaths<-deaths[datemin:datemax]
icu<-icu[datemin:datemax]
icuin<-icuin[datemin:datemax]

icuout<-icuout[datemin:datemax]
deathsc<-deathsc[datemin:datemax]
deathsh<-deathsh[datemin:datemax]

day<-day[datemin:datemax]

if (!is.null(maxdate)){
  date.limit<-which(day==maxdate)
ddms<-ddms[1:date.limit,]
contacts<-contacts[1:date.limit]
hospitalizirani<-hospitalizirani[1:date.limit]
hospitaliziraniin<-hospitaliziraniin[1:date.limit]
hospitaliziraniout<-hospitaliziraniout[1:date.limit]
cases<-cases[1:date.limit]
deaths<-deaths[1:date.limit]
icu<-icu[1:date.limit]
icuin<-icuin[1:date.limit]

icuout<-icuout[1:date.limit]
deathsc<-deathsc[1:date.limit]
deathsh<-deathsh[1:date.limit]

day<-day[1:date.limit]
}


if (!is.null(mindate)){
  date.limit<-which(day==mindate)
  ddms<-ddms[date.limit:length(day),]
  contacts<-contacts[date.limit:length(day)]
  hospitalizirani<-hospitalizirani[date.limit:length(day)]
  hospitaliziraniin<-hospitaliziraniin[date.limit:length(day)]
  hospitaliziraniout<-hospitaliziraniout[date.limit:length(day)]
  cases<-cases[date.limit:length(day)]
  deaths<-deaths[date.limit:length(day)]
  icu<-icu[date.limit:length(day)]
  icuin<-icuin[date.limit:length(day)]
  
  icuout<-icuout[date.limit:length(day)]
  deathsc<-deathsc[date.limit:length(day)]
  deathsh<-deathsh[date.limit:length(day)]
  
  day<-day[date.limit:length(day)]
}


#####



###############Bayes part

Country="Slovenia"
countries <- Country




x.si<-rgammaAlt(1e7,mu.serial,cv.serial)
ecdf.saved.si = ecdf(x.si)

convolution.si = function(u) ( ecdf.saved.si(u))

f.si = rep(0,N2) # f is the probability of dying on day i given infection
f.si[1] = (convolution.si(1.5) - convolution.si(0))
for(i in 2:N2) {
  f.si[i] = (convolution.si(i+.5) - convolution.si(i-.5))
}

serial.interval.r <- data.frame(X=1:N2, fit=rep(NA,N2))
serial.interval.r[,"fit"]<-f.si
serial.interval<-serial.interval.r

# IFR is the overall probability of dying given infection
convolution.1 = function(u) (tau1 * ecdf.saved.1(u))

f1 = rep(0,N2) # f is the probability of dying on day i given infection
f1[1] = (convolution.1(1.5) - convolution.1(0))
for(i in 2:N2) {
  f1[i] = (convolution.1(i+.5) - convolution.1(i-.5))
}

convolution.2 = function(u) (tau2 * ecdf.saved.2(u))

f2 = rep(0,N2) # f is the probability of dying on day i given infection
f2[1] = (convolution.2(1.5) - convolution.2(0))
for(i in 2:N2) {
  f2[i] = (convolution.2(i+.5) - convolution.2(i-.5))
}


convolution.3 = function(u) (tau3 * ecdf.saved.3(u))


f3 = rep(0,N2) # f is the probability of dying on day i given infection
f3[1] = (convolution.3(1.5) - convolution.3(0))
for(i in 2:N2) {
  f3[i] = (convolution.3(i+.5) - convolution.3(i-.5))
}

convolution.4 = function(u) (tau4 * ecdf.saved.4(u))


f4 = rep(0,N2) # f is the probability of dying on day i given infection
f4[1] = (convolution.4(1.5) - convolution.4(0))
for(i in 2:N2) {
  f4[i] = (convolution.4(i+.5) - convolution.4(i-.5))
}

convolution.5 = function(u) (tau5 * ecdf.saved.5(u))


f5 = rep(0,N2) # f is the probability of dying on day i given infection
f5[1] = (convolution.5(1.5) - convolution.5(0))
for(i in 2:N2) {
  f5[i] = (convolution.5(i+.5) - convolution.5(i-.5))
}

convolution.6 = function(u) (tau6 * ecdf.saved.6(u))


f6 = rep(0,N2) # f is the probability of dying on day i given infection
f6[1] = (convolution.6(1.5) - convolution.6(0))
for(i in 2:N2) {
  f6[i] = (convolution.6(i+.5) - convolution.6(i-.5))
}

convolution.7 = function(u) (tau7 * ecdf.saved.7(u))


f7 = rep(0,N2) # f is the probability of dying on day i given infection
f7[1] = (convolution.7(1.5) - convolution.7(0))
for(i in 2:N2) {
  f7[i] = (convolution.7(i+.5) - convolution.7(i-.5))
}


stan_data = list(M=1,N=NULL,deaths=NULL,deathsh=NULL,deathsc=NULL,hosp=NULL,icu=NULL,hospin=NULL,hospout=NULL,icuin=NULL,icuout=NULL,
                 grocery=NULL,parks=NULL,transit=NULL,work=NULL, home=NULL, recreation=NULL,schools=NULL,
                 f1=NULL,f2=NULL,f3=NULL,f4=NULL,f5=NULL,f6=NULL,f7=NULL,
                 N0=6,cases=NULL,SI=serial.interval$fit[1:N2],
                 EpidemicStart = NULL, pop = NULL,contacts=NULL) # N0 = 6 to make it consistent with Rayleigh



dates = list()
reported_cases = list()

deathsc_by_country =deathsh_by_country =deaths_by_country =hosps_by_country=icus_by_country=hospsin_by_country=
  hospsout_by_country=icusin_by_country=icusout_by_country= 
    grocery_by_country=
  parks_by_country=
  transit_by_country=
  work_by_country=schools_by_country=
home_by_country=
recreation_by_country=contacts_by_country=
  list()



d1=data.frame(DateRep=day,Cases=cases,Deaths=deaths,Deathsh=deathsh,Deathsc=deathsc,Hosp=hospitalizirani,Icu=icu,Hospin=hospitaliziraniin,
              Hospout=hospitaliziraniout,Icuin=icuin,Icuout=icuout,
              grocery=ddms$grocery_and_pharmacy_percent_change_from_baseline,
              parks=ddms$parks_percent_change_from_baseline,
              transit=ddms$transit_stations_percent_change_from_baseline,
              work=ddms$workplaces_percent_change_from_baseline,schools=ddms$scools,
              home=ddms$residential_percent_change_from_baseline,
              recreation=ddms$retail_and_recreation_percent_change_from_baseline ,
              contacts=contacts
)
d1$date<-d1$DateRep
d1$t = decimal_date(d1$DateRep)
d1=d1[order(d1$t),]

d1.1<-d1[which(d1$date<=as.Date(datechange1)),]

d1.2<-d1[which(d1$date>as.Date(datechange2)),]

for (m in 1:2){
  stan_data$pop<-c(stan_data$pop,d1_pop)
  if (m==1) d1<-d1.1 else d1<-d1.2

 date_min <- as.Date("2019-12-31") 
if ( d1$DateRep[1]  >  date_min ){
  #print(paste(Country,'In padding'))
  pad_days <- d1$DateRep[1]  - date_min
  pad_dates <- date_min + days(1:pad_days[[1]]-1)
  padded_data <- data.frame("DateRep" =  pad_dates ,
                            "t" = decimal_date(as.Date(pad_dates,format='%d/%m/%Y')),
                            "date" = as.Date(pad_dates,format='%d/%m/%Y'),
                            "Cases" = as.integer(rep(0, pad_days)),
                            "Deaths" = as.integer(rep(0, pad_days)),
                            "Deathsh" = as.integer(rep(0, pad_days)),
                            "Deathsc" = as.integer(rep(0, pad_days)),
                            "Hosp" = as.integer(rep(0, pad_days)),
                            "Icu" = as.integer(rep(0, pad_days)),
                            "Hospin" = as.integer(rep(0, pad_days)),
                            "Hospout" = as.integer(rep(0, pad_days)),
                            "Icuin" = as.integer(rep(0, pad_days)),
                            "Icuout" = as.integer(rep(0, pad_days)),
                            "grocery"=as.integer(rep(0, pad_days)),
                            "parks"=as.integer(rep(0, pad_days)),
                            "transit"=as.integer(rep(0, pad_days)),
                            "work"=as.integer(rep(0, pad_days)),
                            "scools"=as.integer(rep(0, pad_days)),
                            "contacts"=as.integer(rep(0, pad_days)),
"home"=as.integer(rep(0, pad_days)),
"recreation"=as.integer(rep(0, pad_days)),
                              stringsAsFactors=F)
  
  d1 <- bind_rows(padded_data, d1)
 }
index = which(d1$Cases>0)[1]
index1 = which(cumsum(d1$Deaths)>=10)[1] # also 5
index2 = index1-30
if (type==2&m==2) index2<-which(d1$date==as.Date(datechange2))
#print(sprintf("First non-zero cases is on day %d, and 30 days before 10 deaths is day %d",index,index2))


d1=d1[index2:nrow(d1),]






dates[[m]] = d1$date
# hazard estimation
N = length(d1$Cases)
#print(sprintf("%s has %d days of data",Country,N))
forecast = N2 - N


###Roks serial interval code

stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)

#stan_data$timechange<-which(d1$date==datechange)



reported_cases[[m]] = as.vector(as.numeric(d1$Cases))
cases=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))

contacts=c(as.vector(as.numeric(d1$contacts)),rep(0.5,forecast))
contacts_by_country[[m]] = as.vector(as.numeric(d1$contacts))


deaths=c(as.vector(as.numeric(d1$Deaths)),rep(-1,forecast))
deaths_by_country[[m]] = as.vector(as.numeric(d1$Deaths))


deathsh=c(as.vector(as.numeric(d1$Deathsh)),rep(-1,forecast))
deathsh_by_country[[m]] = as.vector(as.numeric(d1$Deathsh))

deathsc=c(as.vector(as.numeric(d1$Deathsc)),rep(-1,forecast))
deathsc_by_country[[m]] = as.vector(as.numeric(d1$Deathsc))


hosps=c(as.vector(as.numeric(d1$Hosp)),rep(-1,forecast))
hosps_by_country[[m]] = as.vector(as.numeric(d1$Hosp))

hospsin=c(as.vector(as.numeric(d1$Hospin)),rep(-1,forecast))
hospsin_by_country[[m]] = as.vector(as.numeric(d1$Hospin))

hospsout=c(as.vector(as.numeric(d1$Hospout)),rep(-1,forecast))
hospsout_by_country[[m]] = as.vector(as.numeric(d1$Hospout))


icus=c(as.vector(as.numeric(d1$Icu)),rep(-1,forecast))
icus_by_country[[m]] = as.vector(as.numeric(d1$Icu))

icusin=c(as.vector(as.numeric(d1$Icuin)),rep(-1,forecast))
icusin_by_country[[m]] = as.vector(as.numeric(d1$Icuin))

icusout=c(as.vector(as.numeric(d1$Icuout)),rep(-1,forecast))
icusout_by_country[[m]] = as.vector(as.numeric(d1$Icuout))

grocery=c(as.vector(as.numeric(d1$grocery)),rep(as.numeric(d1$grocery)[nrow(d1)],forecast))
grocery_by_country[[m]] = as.vector(as.numeric(d1$grocery))

schools=c(as.vector(as.numeric(d1$schools)),rep(as.numeric(d1$schools)[nrow(d1)],forecast))
schools_by_country[[m]] = as.vector(as.numeric(d1$schools))



work=c(as.vector(as.numeric(d1$work)),rep(as.numeric(d1$work)[nrow(d1)],forecast))
work_by_country[[m]] = as.vector(as.numeric(d1$work))


transit=c(as.vector(as.numeric(d1$transit)),rep(as.numeric(d1$transit)[nrow(d1)],forecast))
transit_by_country[[m]] = as.vector(as.numeric(d1$transit))

parks=c(as.vector(as.numeric(d1$parks)),rep(as.numeric(d1$parks)[nrow(d1)],forecast))
parks_by_country[[m]] = as.vector(as.numeric(d1$parks))

home=c(as.vector(as.numeric(d1$home)),rep(as.numeric(d1$home)[nrow(d1)],forecast))
home_by_country[[m]] = as.vector(as.numeric(d1$home))


recreation=c(as.vector(as.numeric(d1$recreation)),rep(as.numeric(d1$recreation)[nrow(d1)],forecast))
recreation_by_country[[m]] = as.vector(as.numeric(d1$recreation))


stan_data$f1 = cbind(stan_data$f1,f1)
stan_data$f2 = cbind(stan_data$f2,f2)
stan_data$f3 = cbind(stan_data$f3,f3)
stan_data$f4 = cbind(stan_data$f4,f4)
stan_data$f5 = cbind(stan_data$f5,f5)
stan_data$f6 = cbind(stan_data$f6,f6)
stan_data$f7 = cbind(stan_data$f7,f7)

stan_data$deaths = cbind(stan_data$deaths,deaths)
stan_data$contacts = cbind(stan_data$contacts,contacts)

stan_data$deathsh = cbind(stan_data$deathsh,deathsh)
stan_data$deathsc = cbind(stan_data$deathsc,deathsc)


stan_data$cases = cbind(stan_data$cases,cases)

stan_data$hosp = cbind(stan_data$hosp,hosps)
stan_data$hospin = cbind(stan_data$hospin,hospsin)
stan_data$hospout = cbind(stan_data$hospout,hospsout)

stan_data$icu = cbind(stan_data$icu,icus)
stan_data$icuin = cbind(stan_data$icuin,icusin)
stan_data$icuout = cbind(stan_data$icuout,icusout)

stan_data$grocery = cbind(stan_data$grocery,grocery/100)

stan_data$parks = cbind(stan_data$parks,parks/100)

stan_data$transit = cbind(stan_data$transit,transit/100)
stan_data$work = cbind(stan_data$work,work/100)

stan_data$home = cbind(stan_data$home,home/100)
stan_data$recreation = cbind(stan_data$recreation,recreation/100)



stan_data$schools = cbind(stan_data$schools,schools)

stan_data$N2=N2
stan_data$x=1:N2

stan_data$N = c(stan_data$N,N)



}




if (spline.type=="bs"){
spline1<-bs(1:stan_data$N2,df=spline.df,degree=degree,Boundary.knots =c(1,stan_data$N[1]))
spline2<-bs(1:stan_data$N2,df=spline.df,degree=degree,Boundary.knots =c(1,stan_data$N[2]))
}
if (spline.type=="ns"){
spline1<-ns(1:stan_data$N2,df=spline.df,Boundary.knots =c(1,stan_data$N[1]))
spline2<-ns(1:stan_data$N2,df=spline.df,Boundary.knots =c(1,stan_data$N[2]))
}


if (spline.df==3){
  stan_data$spline1<-cbind(spline1[,1],spline2[,1])
  stan_data$spline2<-cbind(spline1[,2],spline2[,2])
  
  stan_data$spline3<-cbind(spline1[,3],spline2[,3])
  
   
}


if (spline.df==4){
stan_data$spline1<-cbind(spline1[,1],spline2[,1])
stan_data$spline2<-cbind(spline1[,2],spline2[,2])

stan_data$spline3<-cbind(spline1[,3],spline2[,3])

stan_data$spline4<-cbind(spline1[,4],spline2[,4])
}

if (spline.df==5){
  stan_data$spline1<-cbind(spline1[,1],spline2[,1])
  stan_data$spline2<-cbind(spline1[,2],spline2[,2])
  
  stan_data$spline3<-cbind(spline1[,3],spline2[,3])
  
  stan_data$spline4<-cbind(spline1[,4],spline2[,4])
  stan_data$spline5<-cbind(spline1[,5],spline2[,5])
}

#stan_data$spline<-bs(1:stan_data$N2,df=4)
  
  stan_data$M=2



 

###lock down would happen on 


stan_data$grocery_1<-stan_data$grocery
stan_data$work_1<-stan_data$work
stan_data$home_1<-stan_data$home

ld1<-as.Date("2020-08-10")
idld<-which(dates[[2]]==ld1)
true.ld<-as.Date("2020-10-26")
dif.day<-true.ld-ld1

stan_data$grocery_1[idld:(stan_data$N2-dif.day),2]<-stan_data$grocery[(idld+dif.day):stan_data$N2,2]
stan_data$grocery_1[(stan_data$N2-dif.day+1):stan_data$N2,2]<-stan_data$grocery[stan_data$N2,2]


stan_data$work_1[idld:(stan_data$N2-dif.day),2]<-stan_data$work[(idld+dif.day):stan_data$N2,2]
stan_data$work_1[(stan_data$N2-dif.day+1):stan_data$N2,2]<-stan_data$work[stan_data$N2,2]
 
stan_data$home_1[idld:(stan_data$N2-dif.day),2]<-stan_data$home[(idld+dif.day):stan_data$N2,2]
stan_data$home_1[(stan_data$N2-dif.day+1):stan_data$N2,2]<-stan_data$home[stan_data$N2,2]



stan_data$grocery_2<-stan_data$grocery
stan_data$work_2<-stan_data$work
stan_data$home_2<-stan_data$home

 


stan_data$grocery_2[118:stan_data$N2,2]<-min(stan_data$grocery[,1])
stan_data$work_2[118:stan_data$N2,2]<-min(stan_data$work[,1])
stan_data$home_2[118:stan_data$N2,2]<-max(stan_data$home[,1])


##learn the day at which there was lockdown

lock_down<- which(dates[[1]]=="2020-03-20")

#stan_data$hosp[,2][which(stan_data$hosp[,2]>=stan_data$hosp[,1][lock_down])]
#stan_data$icu[,2][which(stan_data$icu[,2]<=stan_data$icu[,1][lock_down])]
#76 dan bi morali ukrepati to je 14.9. bil pa je 24.10.


ld1<-as.Date("2020-09-14")
idld<-which(dates[[2]]==ld1)
true.ld<-as.Date("2020-10-26")
dif.day<-true.ld-ld1


stan_data$grocery_3<-stan_data$grocery
stan_data$work_3<-stan_data$work
stan_data$home_3<-stan_data$home



stan_data$grocery_3[idld:(stan_data$N2-dif.day),2]<-stan_data$grocery[(idld+dif.day):stan_data$N2,2]
stan_data$grocery_3[(stan_data$N2-dif.day+1):stan_data$N2,2]<-stan_data$grocery[stan_data$N2,2]


stan_data$work_3[idld:(stan_data$N2-dif.day),2]<-stan_data$work[(idld+dif.day):stan_data$N2,2]
stan_data$work_3[(stan_data$N2-dif.day+1):stan_data$N2,2]<-stan_data$work[stan_data$N2,2]

stan_data$home_3[idld:(stan_data$N2-dif.day),2]<-stan_data$home[(idld+dif.day):stan_data$N2,2]
stan_data$home_3[(stan_data$N2-dif.day+1):stan_data$N2,2]<-stan_data$home[stan_data$N2,2]




stan_data$grocery_4<-stan_data$grocery
stan_data$work_4<-stan_data$work
stan_data$home_4<-stan_data$home

stan_data$grocery_4[which(dates[[2]]>=as.Date("2020-10-01"))[1]:stan_data$N2,2]<--0.5
stan_data$work_4[which(dates[[2]]>=as.Date("2020-10-01"))[1]:stan_data$N2,2]<--0.7
stan_data$home_4[which(dates[[2]]>=as.Date("2020-10-01"))[1]:stan_data$N2,2]<-0.3


stan_data$grocery_5<-stan_data$grocery
stan_data$work_5<-stan_data$work
stan_data$home_5<-stan_data$home

stan_data$grocery_5[,2]<-c(stan_data$grocery[76:stan_data$N2,2],rep(stan_data$grocery[stan_data$N2,2],75))
stan_data$work_5[,2]<-c(stan_data$work[76:stan_data$N2,2],rep(stan_data$work[stan_data$N2,2],75))
stan_data$home_5[,2]<-c(stan_data$home[76:stan_data$N2,2],rep(stan_data$home[stan_data$N2,2],75))

premik_vrste<-as.Date("2020-10-26")-as.Date("2020-09-14")



stan_data$grocery_6<-stan_data$grocery
stan_data$work_6<-stan_data$work
stan_data$home_6<-stan_data$home

stan_data$grocery_6[,2]<-c(stan_data$grocery[premik_vrste:stan_data$N2,2],rep(stan_data$grocery[stan_data$N2,2],premik_vrste-1))
stan_data$work_6[,2]<-c(stan_data$work[premik_vrste:stan_data$N2,2],rep(stan_data$work[stan_data$N2,2],premik_vrste-1))
stan_data$home_6[,2]<-c(stan_data$home[premik_vrste:stan_data$N2,2],rep(stan_data$home[stan_data$N2,2],premik_vrste-1))



#######start analysis






 




options(mc.cores = 4)

 if (only.spline==FALSE){
   if (spline.df==3) m = stan_model('stan-models/FinalCounterContactsDF3.stan')
   if (spline.df==4) m = stan_model('stan-models/FinalCounterContactsDF4.stan')
   if (spline.df==5) m = stan_model('stan-models/FinalCounterContactsDF5.stan')  
 } else {
   if (spline.df==3) m = stan_model('stan-models/FinalCounterOnlySplineContactDF3.stan')  
   if (spline.df==4) m = stan_model('stan-models/FinalCounterOnlySplineContactDF4.stan')
   if (spline.df==5) m = stan_model('stan-models/FinalCounterOnlySplineContactDF5.stan')
}



fit = sampling(m,data=stan_data,iter=1600,warmup=800,chains=4,thin=1,control = list(adapt_delta = 0.98, max_treedepth = 15))  

out = rstan::extract(fit)
 