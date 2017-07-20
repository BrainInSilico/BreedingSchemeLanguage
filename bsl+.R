######################## bsl+ ######################
#
# creation date : 09/06/2017
# last update : 26/06/2017
# author : Nicolas Beaume (nicolasbeaume.consultancy@gmail.com)
#
# description : run a simulation of population using BSL (Shiori Yabe, Hiroyoshi Iwata, and Jean-Luc Jannink)
# with extra features
#
#####################################################

#******************** main **************************
source("installBSL.R")
# test if the package is installed, in case not, install it
if(!isTRUE("BreedingSchemeLanguage" %in% .packages(all.available=TRUE))) {
  installBSL()
}
library(BreedingSchemeLanguage)
library(miscTools)
#** set some parameters **
# args = commandArgs(trailingOnly=TRUE)
args <- "bslParam.txt"
args <- readParameters(args)
trainingPop <- 0
currentPop <- 0
#** run simulation **

#### simulate phenotipic selection
# create initial population
simEnv <- defineSpecies(nSim = args$nSim, nMarkers = args$nMarkers, nQTL = args$nQTL)
if(any("qtlEffects"%in%names(args))) {
  setQTLeffect(effects = args$qtlEffects)
}
defineCosts(phenoCost = args$phenoCost, crossCost = args$crossCost, genoCost = args$genoCost)
initializePopulation(nInd = args$initPop)
phenotype()
# prediction made on the same population, the only one available for now
predictValue(popID = currentPop, trainingPopID = trainingPop)
# select best parents
select(nSelect = args$parentSelected)
# F2
cross(nProgeny = args$progenySize)
# predict breeding value of F2 using F1 as training
currentPop <- currentPop+1
# genomic prediction to make sure this function works too
genotype()
predictValue(popID = currentPop, trainingPopID = trainingPop)
# select best individual on F2
select(nSelect = args$parentSelected)
# F3
cross(nProgeny = args$progenySize)
# predict breeding values on F3 using F1 and F2 as training
trainingPop <- c(trainingPop, currentPop)
currentPop <- currentPop+1
phenotype()
predictValue(popID = currentPop, trainingPopID = trainingPop)
# select the best F3
select(nSelect = args$parentSelected)
pdf("testBSL.pdf")
plotData()
dev.off()
selectionRate()
geneticGain()
computeHeterozygosity(popID = unique(simEnv$sims[[1]]$genoRec$popID))
# write total cost, selection rate, genetic gate, heterozygosity and EBV
totalCost <- data.frame(simulation=paste("simulation", 1:simEnv$nSim), totalCost=rep(-1, simEnv$nSim))
ebv <- data.frame(lineID=simEnv$sims[[1]]$predRec$predGID,
                  generationID=simEnv$sims[[1]]$predRec$predNo)
heterozygosity <- data.frame(lineID=simEnv$sims[[1]]$genoRec$GID)
for(i in 1:simEnv$nSim) {
  totalCost[i,2] <- simEnv$sims[[i]]$totalCost
  ebv <- data.frame(ebv, simEnv$sims[[i]]$predRec$predict)
  heterozygosity <- data.frame(heterozygosity, simEnv$sims[[i]]$genoRec$heterozygosityRate)
  write.csv(selRate[[i]], file = paste("selectionRate_sim",i,".csv", sep=""), row.names = T)
  write.csv(genGain[[i]], file = paste("geneticGain_sim",i,".csv", sep=""), row.names = T)
}
write.csv(totalCost, file = "totalCost.csv", row.names = F)
colnames(ebv) <- c("lineID", "generationID", paste("sim", 1:simEnv$nSim, "_EBV", sep=""))
write.csv(ebv, file = "breedingValues.csv", row.names = F)
colnames(heterozygosity) <- c("lineID", paste("sim", 1:simEnv$nSim, "_heterozygosity", sep=""))
write.csv(heterozygosity, file = "heterozygosity.csv", row.names = F)

