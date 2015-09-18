rm(list=ls())
setwd("~/Documents/Deer/Data")

Park = 8 #Change if more parks are added or removed. #MUST DEFINE
source("~/Documents/Deer/Code_Raiho_et_al_2015/Deer_Data_Formatting.R")


DRAW=TRUE #set to TRUE if you want to save figures as pdf automatically
SAVE=TRUE #set to TRUE if you want to save data to file automatically

n.adapt = 4000
n.update = 8000
n.iter = 40000

model.dir = "~/Documents/Deer/Code_Raiho_et_al_2015/"
dump.dir = "~/Documents/Deer/Final_Products/"

source("~/Documents/Deer/Code_Raiho_et_al_2015/Survival.R")
source(paste(model.dir,"Deer_Run_Model.R",sep=""))
source(paste(model.dir,"Deer_Draw_Plots.R",sep=""))
