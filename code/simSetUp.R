#Simulation Setup Codes
simSetUp <- function()
{
	#Load All Necessary sources
	source("code/REML-EM.R")


	#Load All Necessary libraries
  suppressPackageStartupMessages({
      library(MASS)
  })


	###################################################################################################
	#Read Data
	###################################################################################################
	ped <- read.table("code/ped.ped", sep= "\t");
	info <- read.table("code/info.info", head=T);

	return(list(ped=ped, info=info));
}
