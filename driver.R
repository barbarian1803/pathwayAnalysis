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
