source("script/KEGGLoader.R")
source("script/CausalityAnalysisWorkflow.R")

#load data for gene database and matrix for logFC and expression value
load("data/database.dat")
load("data/DESeq2resultBDC.dat")

#get pathway network from KEGG
wntPathway <- getPathwayRelationshipTable("hsa04310")
newWNTPathway <- getEnsemblBasedRelTable(wntPathway)
WNTadj <- createAdjacencyMatrix(newWNTPathway,TRUE)
separatedWNTadj <- separateAdjMatrix(WNTadj) #KEGG network contain subnetwork, separate to individual network
canonPathway <- separatedWNTadj[[3]] #manually choose the canon WNT pathway

logFCData <- FC_BDC_Normal

#check for one gene, LEF1
#ProcessTF is the function that calculate the regression between target and TFs gene
convertGeneID("LEF1","symbol","ensembl")
getTF("ENSG00000138795",canonPathway) #TF for LEF1
resultAnalysis <- ProcessTF("ENSG00000138795",canonPathway)
resultAnalysis[[1]] #inhibitor
resultAnalysis[[1]][1] #significant inhibitor that is found
resultAnalysis[[1]][2] #not significant inhibitor
resultAnalysis[2] #activator
resultAnalysis[[2]][1] #significant activator that is found
resultAnalysis[[2]][3] #regression result
resultAnalysis[[2]][2] #not significant activator


#run for all nodes
allResult <- list()
for(n in colnames(canonPathway)){
  allResult[[n]] <- ProcessTF(n,canonPathway)
}




pi3Pathway <- getPathwayRelationshipTable("hsa04151")
newpi3Pathway <- getEnsemblBasedRelTable(pi3Pathway)
pi3adj <- createAdjacencyMatrix(newpi3Pathway,TRUE)
separatedpi3adj <- separateAdjMatrix(pi3adj)
for(i in (1:length(separatedpi3adj))){
  print(nrow(separatedpi3adj[[i]]))
}

allResultPI3 <- list()
for(n in colnames(separatedpi3adj[[1]])){
  allResultPI3[[n]] <- ProcessTF(n,separatedpi3adj[[1]])
}