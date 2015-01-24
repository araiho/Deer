rm(list=ls())
setwd("~/Documents/MultinomialLikelihood/NPSwtdeer_update_code")

Park = 8 #Change if more parks are added or removed. #MUST DEFINE
source("~/Documents/Deer/Code_Raiho_et_al_2015/Deer_Data_Formatting.R")

DRAW=TRUE #set to TRUE if you want to save figures as pdf automatically
SAVE=TRUE #set to TRUE if you want to save data to file automatically

n.adapt = 2
n.update = 4
n.iter = 20

model.dir = "~/Documents/Deer/Code_Raiho_et_al_2015/"
dump.dir = "~/Documents/Deer/Dump/"

source(paste(model.dir,"Deer_Run_Model.R",sep=""))
source(paste(model.dir,"Deer_Draw_Plots.R",sep=""))
