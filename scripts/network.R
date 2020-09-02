# Network code

plotInputN <- function(KinHits, hit_col, int_col, n_size, n_text, layout){
  HitsColor=hit_col #user input
  InterColor=int_col #user input
  NodeSize=n_size #user input
  NodeTextSize= abs(n_text - 11) #user input
  nodes2 <- readRDS("./data/sup/FinNodes_Final_FF.rds")
  edges <- readRDS("./data/sup/FinEdges_Final_FF.rds")
  
  nodes2 %>% dplyr::filter(FinName %in% KinHits) %>% pull(FinName) -> sigHITS
  edges %>% dplyr::filter(Source %in% sigHITS | Target %in% sigHITS) %>% dplyr::filter(Source != Target) -> modEdges 
  
  modsources <- pull(modEdges, Source)
  modtargets <- pull(modEdges, Target)
  
  modALL <- unique(c(modsources,modtargets))
  
  nodes2 %>% dplyr::filter(FinName %in% modALL) -> nodesF
  
  listoFkin <- pull(nodesF, FinName)
  
  edges %>% dplyr::filter(Source %in% nodesF$FinName & Target %in% nodesF$FinName) %>% dplyr::filter(Source != Target) -> modEdges
  
  modEdges %>% mutate(line = ifelse(Source %in% sigHITS | Target %in% sigHITS, 2,1 )) -> modEdges
  
  modsources <- pull(modEdges, Source)
  modtargets <- pull(modEdges, Target)
  
  modALL <- c(modsources,modtargets)
  as.data.frame(table(modALL)) -> concts 
  
  concts %>% dplyr::rename(FinName = modALL) -> concts
  
  concts$FinName <-  as.character(concts$FinName)
  
  right_join(nodesF,concts) -> nodesF
  nodesF %>% mutate(cl = ifelse(FinName %in% sigHITS, HitsColor, InterColor)) -> nodesF
  # filter low freqs 
  
  nodesF %>% dplyr::filter(Freq>=1|cl==HitsColor) %>% 
    pull(FinName) -> FinKinases
  
  modEdges %>% dplyr::filter(Source %in% FinKinases & Target %in% FinKinases) -> modEdges
  nodesF %>% dplyr::filter(FinName %in% FinKinases) %>% mutate(Freq = ifelse(Freq < 4, 4, Freq)) -> nodesF
  
  net <- graph_from_data_frame(d=modEdges, vertices=nodesF, directed=T) 
  net <- igraph::simplify(net,remove.loops = F, remove.multiple = F)
  
  V(net)$size = log2(V(net)$Freq)*NodeSize
  colrs <- c(HitsColor, InterColor)
  V(net)$color <- V(net)$cl
  
  colrs2 <- c("gray", "black")
  E(net)$color <- colrs2[E(net)$line]
  plot(net, edge.arrow.size=.05,vertex.label=V(net)$FinName,vertex.label.color = "black",
       vertex.label.cex=log2(V(net)$Freq)/NodeTextSize, layout = layout_in_circle)
  
  # We can even set the network layout:
  if(layout=="Circle"){
    plot(net, edge.arrow.size=.05,vertex.label=V(net)$FinName,vertex.label.color = "black",
         vertex.label.cex=log2(V(net)$Freq)/NodeTextSize, layout = layout_in_circle)
    
  }else if(layout=="Fruchterman-Reingold"){
    plot(net, edge.arrow.size=.05,vertex.label=V(net)$FinName,vertex.label.color = "black",
         vertex.label.cex=log2(V(net)$Freq)/NodeTextSize, layout = layout_with_fr)
    
  }else if(layout=="Kamada Kawai"){
    plot(net, edge.arrow.size=.05,vertex.label=V(net)$FinName,vertex.label.color = "black",
         vertex.label.cex=log2(V(net)$Freq)/NodeTextSize, layout = layout_with_kk)
    
  }else if(layout=="LGL"){
    plot(net, edge.arrow.size=.05,vertex.label=V(net)$FinName,vertex.label.color = "black",
         vertex.label.cex=log2(V(net)$Freq)/NodeTextSize, layout = layout_with_lgl)
    
  }
  
}