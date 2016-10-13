library(DESeq2)

#calculate TFs for target gene
ProcessTF <- function(Target,graph){
  message(paste("Processing gene ",Target))
  allTF <- getTF(Target,graph)
  if(is.null(allTF)){
    print("No TF found")
    return(NULL)
  }
  inhibitor <- c()
  activator <- c()
  unknown <- c()
  
  for(tf.id in allTF){
	#filter out data that is not exist in the expression table
	#it is possible data from KEGG contain gene/compound that is not defined in ensembl database
    if(identical(tf.id,character(0))){
      next
    }
    if(!tf.id %in% rownames(logFCData)){
      next
    }
	#filter TF that is not differentially expressed
    if(abs(logFCData[tf.id,"log2FoldChange"])<0.65){
      next
    }
	#check if activator or inhibitor
    type <- getTFType(tf.id,Target,graph)
    
    if(type%in%c(3,5,9)){
      activator <- c(activator,tf.id)
    }else if(type%in%c(4,6,10)){
      inhibitor <- c(inhibitor,tf.id)
    }else{
      unknown <- c(unknown,tf.id)
    }
  }
  
  #inhibitor variable will contain with TF that is categorized as inhibitor
  #activator variable will contain with TF that is categorized as activator
  inhibitorAnalysis <- checkInhibitor(inhibitor,Target)  
  activatorAnalysis <- checkActivator(activator,Target)
  output <- list(inhibitorAnalysis,activatorAnalysis)
  output
}

checkInhibitor <- function(inhibitor,Target){
	#check the inputted inhibitor and splitted the inputted inhibitor into 2, significant and not
  logTarget <- logFCData[Target,"log2FoldChange"][1]
  
  significant <- c()
  nonsignificant<-c()
  
  for(Gene in inhibitor){
    corr <- cor(AssayNormalCancer[Gene,],AssayNormalCancer[Target,])
    if(corr < (-0.6)){
      significant <- c(significant,Gene)
    }else{
      nonsignificant <- c(nonsignificant,Gene)      
    }
  }
  
  #return value of this function is a list contains 2 vectors
  #1st vector is the significant inhibitor that is found after calculation
  #2nd vector is the remaining inhibitor that we found not significant
  cause <- list()
  cause[[1]]<-significant
  cause[[2]]<-nonsignificant
  cause
}

checkActivator <- function(activator,Target){
  #do the regression of the inputted activator
  #the outpur is a list contains 3 vectors
  #1st vector is the gene names that are found significant (used in the linear regression model)
  #2nd vector is the gene names that are found not significant (not used in the linear regression model)
  #3rd vector is the summary of the linear regression of the most significant model (genes from 1st vector)
  
  
  oriactivator <- activator
  #filter the inputted activator
  activator <- filterActivator(activator,Target)
  
  if(length(activator)<1){
    res <- list()
    res[[1]]<-c()
    res[[2]]<-oriactivator
    res[[3]]<-NULL
    return(res)
  }
  
  #use ensemblID, used to be changing the symbol but since already use the ensemblID, just ignore
  activatorEnsembl <- c()
  for (a in activator){
    activatorEnsembl <- c(activatorEnsembl,a)
  }
  
  activator <- activatorEnsembl
  max <- -1
  selected <- c()
  LMSummary <- NULL
  for(i in length(activator):1){
	#brute force to get all possible combination of TFs
    list <- enum.choose(activator,i)
    for(j in 1:length(list)){
      formula <- paste(Target,paste(list[[j]],collapse="+"),sep="~")
      formula <- as.formula(formula)
      LM <- summary(lm(formula,as.data.frame(t(AssayNormalCancer))))
      rsqr <- LM$adj.r.squared
      coeff <- LM$coefficients
      if( rsqr>max & !AnyNegativeCoeff(coeff[-1,1]) & AnySignificant(coeff[-1,4]) ){
        max <- rsqr
        selected <- list[[j]]
        LMSummary <- LM
      }
    }
  }
  selectedSymbol <- c()
  for(n in selected){
    selectedSymbol <- c(selectedSymbol,n)
  }
  res <- list()
  res[[1]] <- selectedSymbol
  res[[2]]<-setdiff(oriactivator,selectedSymbol)
  res[[3]]<-LMSummary
  res
}

#remove activator which has negative correlation with its target
filterActivator <- function(activator,Target){
  out <- c()
  for( act in activator){
    corr <- cor(AssayNormalCancer[act,],AssayNormalCancer[Target,])
    if(corr > 0){
      out <- c(out,act)
    }
  }
  out
}

enum.choose <- function(x, k) {
  if(k > length(x)) stop('k > length(x)')
  if(choose(length(x), k)==1){
    list(as.vector(combn(x, k)))
  } else {
    cbn <- combn(x, k)
    lapply(seq(ncol(cbn)), function(i) cbn[,i])
  }
}

AnyNegativeCoeff <- function(vector){
  out <- FALSE
  for(v in vector){
    if(v<0){
      out <- TRUE
      break
    }
  }
  out
}

AnySignificant <- function(vector){
  out <- FALSE
  for(v in vector){
    
    if(!is.nan(v)&v<=0.05){
      out <- TRUE
      break
    }
  }
  out
}