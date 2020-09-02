
#Function to get a score for each kinase for each iteration
iterationScore <- function(sampleList, kinaseSep, uniqueKin){
  
  splitKinase=NULL
  for(i in 1:length(sampleList)){
    splitKinase=paste0(splitKinase, sampleList[i], sep=kinaseSep)
  }
  splitKinase2=unlist(strsplit(splitKinase, split=kinaseSep))
  returnVector=c()
  for(i in 1:length(uniqueKin)){
    temp=paste("^",uniqueKin[i],"$", sep="")
    returnVector[i]=length(grep(temp, splitKinase2))
  }
  return(returnVector)
  
}