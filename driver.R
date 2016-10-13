source("script/KEGGLoader.R")
source("script/CausalityAnalysisWorkflow.R")

#load data for gene database and matrix for logFC and expression value
load("data/database.dat")
load("data/DESeq2resultBDC.dat")

#get pathway network from KEGG
wntPathway <- getPathwayRelationshipTable("hsa04310")
newWNTPathway <- getEnsemblBasedRelTable(wntPathway)
WNTadj <- createAdjacencyMatrix(newWNTPathway,TRUE)
separatedWNTadj <- separateAdjMatrix(WNTadj)
canonPathway <- separatedWNTadj[[3]]

logFCData <- FC_BDC_Normal
#check for one gene, LEF1
#ProcessTF is the function that calculate the regression between target and TFs gene
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