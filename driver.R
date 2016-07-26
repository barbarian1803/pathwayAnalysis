source("script/KEGGLoader.R")
load("data/database.dat")
load("data/DESeq2resultBDC.dat")


wntPathway <- getPathwayRelationshipTable("hsa04310")
newWNTPathway <- getEnsemblBasedRelTable(wntPathway)
WNTadj <- createAdjacencyMatrix(newWNTPathway,TRUE)
separatedWNTadj <- separateAdjMatrix(WNTadj)


logFCData <- FC_BDC_Normal
source("script/CausalityAnalysisWorkflow.R")
convertGeneID("LEF1","symbol","ensembl")
resultAnalysis <- ProcessTF("ENSG00000138795",canonPathway)


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