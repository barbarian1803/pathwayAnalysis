library(DESeq2)


ProcessTF <- function(Target,graph){
  message(paste("Processing gene ",Target))
  allTF <- getTF(Target,graph)
  if(is.null(allTF)){
    print("masuk null")
    return(NULL)
  }
  inhibitor <- c()
  activator <- c()
  unknown <- c()
  
  for(tf.id in allTF){
    if(identical(tf.id,character(0))){
      next
    }
    if(!tf.id %in% rownames(logFCData)){
      next
    }
    if(abs(logFCData[tf.id,"log2FoldChange"])<0.65){
      next
    }
    type <- getTFType(tf.id,Target,graph)
    
    if(type%in%c(3,5,9)){
      activator <- c(activator,tf.id)
    }else if(type%in%c(4,6,10)){
      inhibitor <- c(inhibitor,tf.id)
    }else{
      unknown <- c(unknown,tf.id)
    }
    
  }
  inhibitorAnalysis <- checkInhibitor(inhibitor,Target)
  activatorAnalysis <- checkActivator(activator,Target)
  output <- list(inhibitorAnalysis,activatorAnalysis)
  output
}

checkInhibitor <- function(inhibitor,Target){
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
  
  cause <- list()
  cause[[1]]<-significant
  cause[[2]]<-nonsignificant
  cause
}

checkActivator <- function(activator,Target){
  
  oriactivator <- activator
  activator <- filterActivator(activator,Target)
  
  if(length(activator)<1){
    res <- list()
    res[[1]]<-c()
    res[[2]]<-oriactivator
    res[[3]]<-NULL
    return(res)
  }
  activatorEnsembl <- c()
  for (a in activator){
    activatorEnsembl <- c(activatorEnsembl,a)
  }
  activator <- activatorEnsembl
  max <- -1
  selected <- c()
  LMSummary <- NULL
  for(i in length(activator):1){
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