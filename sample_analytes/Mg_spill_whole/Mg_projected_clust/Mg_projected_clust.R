wrapper_projected<-function(chunk_ID){
  ## Loading required libraries
  library(geosphere)
  library(network)
  library(igraph)
  library(mapdata)
  library(intergraph)
  require(sna)
  require(maps)
  library(GGally)
  library(MASS)
  library(foreach)
  library(doParallel)
  library(data.table)
  
  graph_subset_indicator<-F
  
  file_path<-"/storage/home/aua257/work/spill_datasets_run/Mg_spill_whole/"

  ## Loading the dataframe "df_polluter_nodeID_aggregated"
  load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
  
  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))
  ## Loading the "total_edgelist_character_modified"
  load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))
  
  ## Defining a wrapper function to calculate flow distance of all polluter nodes in parallel
  wrapper_flow_dist_cal<-function(df_polluter,anpoll_edgelist, shortest_path_anpoll_edgelist, total_edgelist_character_modified, from_indicator, to_indicator, file_path){
    ## Defining a function to calculate between a test node and a current node. The current nodes are usually polluter nodes.
    flow_dist_cal<-function(anpoll_edgelist, shortest_path_anpoll_edgelist, total_edgelist_character_modified, test_node_ID, current_node_ID, from_indicator, to_indicator, file_path){
      ## Getting the edge row ID
      if(from_indicator){
        edge_row_ID<-which((anpoll_edgelist[,1]==test_node_ID)&(anpoll_edgelist[,2]==current_node_ID))
      }else if(to_indicator){
        edge_row_ID<-which((anpoll_edgelist[,1]==current_node_ID)&(anpoll_edgelist[,2]==test_node_ID))
      }
      ## Initializing and loop to fill up the matrix of edges that are formed by the nodeIDs in the path of test and current node ID
      flow_path_nodeIDs<-matrix(NA_character_,nrow = (length(shortest_path_anpoll_edgelist[[edge_row_ID]])-1),ncol = 2)
      for (i in 1:(length(shortest_path_anpoll_edgelist[[edge_row_ID]])-1)){
        flow_path_nodeIDs[i,]<-c(shortest_path_anpoll_edgelist[[edge_row_ID]][i],shortest_path_anpoll_edgelist[[edge_row_ID]][i+1])
      }
      ## Loading the total_edgelist_character_modified
      load(file = paste0(file_path,"common_files_modified/stream_path_dist_vec.RData"))
      ## Initializing the flow distance in meters
      flow_dist_m<-0
      ## Loop to add up the flow distance of all edges in the flow_path_nodeIDs
      for(i in 1:nrow(flow_path_nodeIDs)){
        flow_dist_m<-flow_dist_m+stream_path_dist_vec[which((total_edgelist_character_modified[,1]==flow_path_nodeIDs[i,1])&(total_edgelist_character_modified[,2]==flow_path_nodeIDs[i,2]))]
      }
      ## flow distance in km
      flow_dist_km<-flow_dist_m/1000
      return(flow_dist_km)
    }
    ## Extracting all connected nodeIDs
    if(from_indicator){
      connected_nodeIDs<-anpoll_edgelist[which(anpoll_edgelist[,2]==df_polluter$nodeID[index]),1]
    }else if(to_indicator){
      connected_nodeIDs<-anpoll_edgelist[which(anpoll_edgelist[,1]==df_polluter$nodeID[index]),2]
    }
    ## Initializing the flow_dist_vec as NA. If there are no connected nodes, this would remain NA
    flow_dist_vec<-NA
    if(length(connected_nodeIDs)>0){
      flow_dist_vec<-rep(NA_integer_,length(connected_nodeIDs))
      #print(length(connected_nodeIDs))
      for(j in 1:length(connected_nodeIDs)){
        flow_dist_vec[j]<-flow_dist_cal(anpoll_edgelist = anpoll_edgelist,shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,total_edgelist_character_modified = total_edgelist_character_modified,test_node_ID = connected_nodeIDs[j],current_node_ID = df_polluter$nodeID[index],from_indicator = from_indicator,to_indicator = to_indicator,file_path = file_path)
        #print(j)
      }
    }
    return(flow_dist_vec)
  }
  
  ## Loading the dataframe "df_anpoll_processed"
  load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
  ## Loading the vector of projected node IDs
  load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
  
  ## Writing a function to generate chunk_list from projected_nodeIDs_vec
  projected_chunks_generator<-function(projected_nodeIDs_vec){
    n_total<-length(projected_nodeIDs_vec)
    n_chunks<-80
    chunk_size<-ceiling(n_total/n_chunks)
    indices_all<-1:n_total
    indices_chunk_list<-split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))
    return(indices_chunk_list)
  }
  
  chunk_list<-projected_chunks_generator(projected_nodeIDs_vec = projected_nodeIDs_vec)
  
  ## Getting the flow dist from list for projected node IDs of polluters
  df_projected_nodeIDs<-data.frame("nodeID"=projected_nodeIDs_vec,stringsAsFactors = F)
  flow_dist_from_list_projected<-list()
  cl <- makeCluster(5)
  registerDoParallel(cl)
  flow_dist_from_list_projected<-foreach (index = chunk_list[[chunk_ID]])%dopar% {
    wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
  }
  stopCluster(cl)
  
  str(flow_dist_from_list_projected)
  names(flow_dist_from_list_projected)<-projected_nodeIDs_vec[chunk_list[[chunk_ID]]]
  save(flow_dist_from_list_projected,file = paste0(file_path,"inference/projected_results/flow_dist_from_list_projected_",chunk_ID,".RData"))
  print("from_list_projected is done")
}

#wrapper_projected(chunk_ID = 100)

