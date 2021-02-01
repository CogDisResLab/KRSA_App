# krsa 

krsa <-
  function (file ,iterations, numSigSubs, sigSubs) {
    withProgress(message = 'KRSA', value = 0, {
      #incProgress(1/5, detail = "Reading file")
      subKinase <- read.table(file, header = TRUE, sep = "\t", as.is=T, quote = "\"", check.names=FALSE, stringsAsFactors=FALSE)
      kinaseSep=" "
      subKinaseOnly=toupper(subKinase[,2])
      
      splitKinase=NULL
      #incProgress(2/5, detail = "Finding unique kinases")
      for(i in 1:length(subKinaseOnly)){
        splitKinase=paste0(splitKinase, subKinaseOnly[i], sep=kinaseSep )
      }
      uniqueKin=unique(unlist(strsplit(splitKinase, split=kinaseSep)))
      scores=data.frame(matrix(ncol = iterations, nrow = length(uniqueKin)), row.names=uniqueKin)
      
      for(i in 1:iterations){
        incProgress((1/iterations), detail = paste0(i, "/", iterations))
        sampleList=sample(subKinaseOnly, numSigSubs) 
        scores[,i]=iterationScore(sampleList, kinaseSep, uniqueKin)
      }
      
      #incProgress(4/5, detail = "Averaging")
      scores=cbind(scores, rowMeans(scores)) #get the average kinase amount
      yourSubs=subset(subKinase, subKinase$Substrates %in% sigSubs)
      yourScores=data.frame(matrix(ncol = 1, nrow = length(uniqueKin)), row.names=uniqueKin)
      yourScores[,1]=iterationScore(yourSubs[,2], kinaseSep, uniqueKin)
      
      
      #section to make overall table
      #incProgress(5/5, detail = "Creating results")
      overallTable=as.data.frame(cbind(Kinase=uniqueKin, DxAverage=yourScores[,1], ReSamplingAverage=scores[,iterations+1]))
      overallTable[,2]=as.numeric(as.character(overallTable[,2]))
      overallTable[,3]=as.numeric(as.character(overallTable[,3]))
      overallTable=cbind(overallTable, AbsoluteDiference=abs(overallTable[,2]-overallTable[,3]))
      overallTable=cbind(overallTable, StandardDeviation=unlist(apply(scores[,1:iterations], 1, sd)))
      overallTable=cbind(overallTable, Minus2Sd=overallTable[,3]-2*(overallTable[,5]) , Plus2SD=overallTable[,3]+2*(overallTable[,5]))
      overallTable=cbind(overallTable, Down=(ifelse(overallTable[,2]<overallTable[,6], "YES", "x")), Up=(ifelse(overallTable[,2]>overallTable[,7], "YES", "x")))
      overallTable=cbind(overallTable, ZScore=abs((overallTable[,2]-overallTable[,3])/overallTable[,5]))
      overallTable2=overallTable[order(overallTable$ZScore, decreasing=TRUE),]
      
      #kittyYES<<-overallTable2
      print("Done")
    })
    return(list(overallTable2, scores))
  }

# new krsa

samplingPep <- function(x,CovFile,sum_num) {
  CovFile %>%
    group_by(Kin) %>%
    summarise(
      counts = sum(Substrates %in% sample(CovFile$Substrates %>% unique(),sum_num)), .groups = 'drop'
    ) %>% mutate(itr = x) -> res
  
  return(res)
}


krsa_2 <- function(x,CovFile,peps) {
  
  # file %>% separate_rows(Kinases) %>%
  #   rename(Kin = Kinases) -> CovFile
  
  purrr::map_df(1:x,samplingPep,CovFile,length(peps)) -> res
  
  res %>% group_by(Kin) %>% summarise(SamplingAvg = mean(counts), SD= sd(counts)) -> avrages
  
  obs <- CovFile %>%
    group_by(Kin) %>%
    summarise(
      Observed = sum(Substrates %in% peps)
    )
  
  left_join(avrages, obs) %>% mutate(Z = (Observed-SamplingAvg)/SD) %>% arrange(desc(abs(Z))) %>% 
    filter(!Kin %in% c("BARK1", "VRK2")) %>% 
    select(Kin, Observed, SamplingAvg, SD, Z) %>% rename(Kinase = Kin) -> fin
  
  
  return(list(
    res = res,
    krsa_table = fin
  ))
}
