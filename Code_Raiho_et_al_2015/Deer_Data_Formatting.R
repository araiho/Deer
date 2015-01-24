
#####
#####
##### Data Editing File
##### "Forecasting the Effects of Fertility Control on Overabundant Ungulates" 
##### Created by: Ann Raiho 
##### Last Edited on: 26 November 2014
#####
#####

library(reshape) #package for "melting" and "casting" data

#####
##### Importing data sheets #####
#####

anti=(read.csv("ANTI.csv"))
cato=(read.csv("CATO.csv"))
choh=(read.csv("CHOH.csv"))
gwmp=(read.csv("GWMP.csv"))
mana=(read.csv("MANA.csv"))
mono=(read.csv("MONO.csv"))
nace=(read.csv("NACE.csv"))
prwi=(read.csv("PRWI.csv"))
rocr=(read.csv("ROCR.csv"))

#####
##### Combining into one sheet #####
#####

raw.pop.data=rbind(anti,cato,choh,gwmp,mana,mono,nace,prwi,rocr)

#####
##### Deleting distance sampling data not used #####
#####

raw.pop.data$Event_Notes<-NULL
raw.pop.data$Temperature<-NULL
raw.pop.data$Moon_Phase<-NULL
raw.pop.data$Sky<-NULL
raw.pop.data$Wind.Speed<-NULL
raw.pop.data$Verified<-NULL
raw.pop.data$Date_Time_Entered<-NULL
raw.pop.data$Time<-NULL
raw.pop.data$Mileage<-NULL
raw.pop.data$Import_Event_ID<-NULL
raw.pop.data$DistanceFromRoad<-NULL
raw.pop.data$Angle<-NULL
raw.pop.data$DistanceToDeer<-NULL
raw.pop.data$Route_Length<-NULL
raw.pop.data$Route_Length<-NULL
raw.pop.data$Route.Section<-NULL
raw.pop.data$Data_ID<-NULL
raw.pop.data$Event_ID<-NULL
raw.pop.data$Event_Group_ID<-NULL
raw.pop.data$Start_Date<-NULL
raw.pop.data$Park_Code<-NULL

#####
##### Categorical Data #####
#####

raw.pop.data=as.data.frame(raw.pop.data) #coercing to correct format
raw.noFOWA=raw.pop.data[raw.pop.data[,1]!="FOWA",] #removing park without enough data
raw.noFOWA=raw.noFOWA[raw.noFOWA[,1]!="PISC",] #removing park with hunting in vicinity
raw.noFOWA=raw.noFOWA[raw.noFOWA[,1]!="CATO",] #removing park with culling
raw=as.data.frame(raw.noFOWA)

#creating factors for parks and years to do indexing in model
raw[,1]<-as.numeric(factor(raw[,1],levels=unique(raw.noFOWA[,1]),ordered=TRUE))
raw[,2]<-as.numeric(as.factor(raw[,2]))
raw[is.na(raw)]<-0

#creating dataset of deer in categorized each year from raw distance sampling data
raw.melt=melt(raw,id.var=1:2,na.rm=TRUE)
y.alpha2=as.data.frame(cast(raw.melt,Route+Year~variable,sum))
y.alpha2[,c(4,6)]<-y.alpha2[,c(6,4)] #switch Bucks and fawns
colnames(y.alpha2)=c("Route","Year","Group_Index","fawns","Doe","bucks","unknown","Total_Group") 
y.alpha2[,8]=rowSums(y.alpha2[,4:6])
y.alpha2=y.alpha2[y.alpha2[,8]!=0,]

#creating dataset of deer in categories from raw distance sampling data NOT SUMMED for each year
y.alpha=raw[raw[,4]!=0 | raw[,5]!=0 | raw[,6]!=0,]
y.alpha[,c(4,6)]<-y.alpha[,c(6,4)] #switch Bucks and fawns
colnames(y.alpha)=c("Route","Year","Group_Index","fawns","Doe","bucks","unknown","Total_Group") 
y.alpha[,8]=rowSums(y.alpha[,4:6])
y.alpha=y.alpha[y.alpha[,8]!=0,]

#####
##### Density Estimates from Program DISTANCE #####
#####

pop.data=(read.csv("deer-brmd data request4.csv")) 
names(pop.data)=c("park","year","y.park","y.yr","density", "95lower","95upper","stderr","cvs","cvt","#ofgroups","%unclassified","%bucks","%does","%fawns","fawn/doe","buck/doe") 

N.obs=as.data.frame(pop.data)#coercing to correct type of object
N.obs = N.obs[N.obs[,1]!="CATO",]#removing park with culling
N.obs = N.obs[N.obs[,1]!="PISC",]#removing park with hunting in the area
N.obs[,1]=as.numeric(factor(N.obs[,1],levels=unique(raw.noFOWA[,1]),ordered=TRUE))#making factors for parks for indexing
N.obs[,c(1,3)]<-N.obs[,c(3,1)]
N.obs[is.na(N.obs)]<-0
N.obs=as.matrix(N.obs)#coercing to correct type of object

#area=as.list(c(13.15,23.35,81.91,28.92,18.3,6.07,4.76,17.06,65.09,11.41)) #areas for all parks 
area=as.list(c(13.15,81.91,28.92,18.3,6.07,4.76,65.09,11.41)) #areas for 8 parks used.

#formatting additional data for years 2012 and 2013. Contact ann.raiho@gmail.com to add more data.
pop.data1=(read.csv("BRMD Deer REQUEST 2012 2013.csv"))
pop.data1=pop.data1[pop.data1[,1]!="FTWA",]#removing park without enough data
pop.data1=pop.data1[pop.data1[,1]!="CATO",]#removing park with culling
pop.data1=pop.data1[pop.data1[,1]!="PISC",]#removing park with hunting in area
pop.data1[,1]=as.numeric(factor(pop.data1[,1],levels=unique(raw.noFOWA[,1]),ordered=TRUE))#making factors for parks for indexing
new.year.factor=c(rep(13,Park),rep(12,Park))#making factors for years for indexing

#adding the additional density data
pop.data1=as.matrix(rbind(pop.data1[3:10,],pop.data1[12:nrow(pop.data1),]))
new.pop=matrix(0,nrow=nrow(pop.data1),ncol=ncol(N.obs))
new.pop[,3]=as.numeric(pop.data1[,1])
new.pop[,4]=as.numeric(new.year.factor)
new.pop[,5]=as.numeric(pop.data1[,3])
new.pop[,8]=as.numeric(pop.data1[,5])

#creating final density matrices
N.obs=as.data.frame(rbind(N.obs,new.pop))#has year 1
N.obsno1s=as.matrix(N.obs[N.obs[,4]!=1,])#does not have year 1 #for estimating 

if(length(unique(N.obs[,3])) == Park) {
  print("parks stated equals parks in density data")
  } else {
    print("need to reformat data to match number of parks")
    }#checking to see if we have the right number of parks in data


#adding the additional categorical data
new.alpha=matrix(0,nrow=nrow(pop.data1),ncol=ncol(y.alpha2))
new.alpha[,1]=as.factor(pop.data1[,1])
new.alpha[,2]=as.numeric(new.year.factor)
new.alpha[,4]=as.numeric(pop.data1[,10])
new.alpha[,5]=as.numeric(pop.data1[,9])
new.alpha[,6]=as.numeric(pop.data1[,8])
new.alpha[,8]=rowSums(new.alpha[,4:6])
new.alpha=as.data.frame(new.alpha)
names(new.alpha)<-names(y.alpha2)

#creating final categorical matrices
y.alpha2=as.data.frame(rbind(y.alpha2,new.alpha))

if(length(unique(y.alpha2[,1])) == Park) {
  print("parks stated equals parks in categorical data")
} else {
  print("need to reformat data to match number of parks")
}##checking to see if we have the right number of parks in data


