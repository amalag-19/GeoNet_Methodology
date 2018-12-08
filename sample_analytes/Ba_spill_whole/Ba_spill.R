########################################################################################################
############################################ Ba Spill Analysis ##########################################
########################################################################################################
########################################################################################################

file_path<-"/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/River network project/Analysis for Spill Paper/Ba_spill_whole/"

## Sourcing the modular functions for analysis
source(file = "/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/River network project/BaSu_network_v15_modular_functions.R")

dir.create(path = paste0(file_path,"common_files"))
dir.create(path = paste0(file_path,"common_files_modified"))
dir.create(path = paste0(file_path,"analyte_files"))
dir.create(path = paste0(file_path,"polluter_files"))
dir.create(path = paste0(file_path,"anpoll_files"))
dir.create(path = paste0(file_path,"inference"))
dir.create(path = paste0(file_path,"inference/projected_results"))

########################################################################################################
## Loading the river net shape file
load(file = "~/Box Sync/PSU/Fall 2018/Geoscience_Research/River network project/NetStreams1998/shape.RData")

## Running the function "shapefile_pre_processing" to save the first three important files
shapefile_pre_processing_list<-shapefile_pre_processing(shape_obj = shape,output_path = file_path)

####################################################
## Getting analyte dataframe
## Analyte preprocessing; Reading and cleaning the analyte dataset
df_analyte_raw<-read.csv(file = paste0(file_path,"data/Ba_data.csv"),stringsAsFactors = F)
str(df_analyte_raw)

## Preprocessing the analyte dataframe
df_analyte_preprocessed<-df_analyte_raw[,c(1,2,3)]
df_analyte_preprocessed$date<-as.Date(df_analyte_raw$obsDate,"%Y-%m-%d")
names(df_analyte_preprocessed)<-c("conc","lon","lat","date")
str(df_analyte_preprocessed)
## saving this dataframe
save(df_analyte_preprocessed,file = paste0(file_path,"analyte_files/df_analyte_preprocessed.RData"))

####################################################
## Preprocessing the polluter dataframe
df_polluter_raw<-read.csv(file = paste0(file_path,"data/Spill_data.csv"),stringsAsFactors = F)
str(df_polluter_raw)

## Preprocessing the analyte dataframe
df_polluter_preprocessed<-df_polluter_raw[,c(1,2)]
df_polluter_preprocessed$date<-as.Date(df_polluter_raw$Date,format="%m/%d/%y")
str(df_polluter_preprocessed)
head(df_polluter_preprocessed)
## saving this dataframe
save(df_polluter_preprocessed,file = paste0(file_path,"polluter_files/df_polluter_preprocessed.RData"))

####################################################
## Combining the analyte and polluter dataframes
df_anpoll_preprocessed<-data.frame(matrix(NA,nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed),4))
df_anpoll_preprocessed<-data.frame(df_analyte_preprocessed,"anpoll_indicator"=rep("analyte",nrow(df_analyte_preprocessed)),stringsAsFactors = F)
df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),c(2,3,4)]<-df_polluter_preprocessed
df_anpoll_preprocessed[(nrow(df_analyte_preprocessed)+1):(nrow(df_analyte_preprocessed)+nrow(df_polluter_preprocessed)),"anpoll_indicator"]<-"polluter"

str(df_anpoll_preprocessed)
head(df_anpoll_preprocessed)
tail(df_anpoll_preprocessed)
df_anpoll_preprocessed[(nrow(df_anpoll_preprocessed)-nrow(df_polluter_preprocessed)):nrow(df_anpoll_preprocessed),]

## saving this dataframe
save(df_anpoll_preprocessed,file = paste0(file_path,"anpoll_files/df_anpoll_preprocessed.RData"))

########################################################################################################
## Loading the dataframe "df_anpoll_preprocessed"
load(file = paste0(file_path,"anpoll_files/df_anpoll_preprocessed.RData"))

## load(file = paste0(file_path,"anpoll_files/dist_mapped.RData"))

## Performing the two step dynamic mapping procedure that can map all C-PP locations to the base river map
two_step_anpoll_mapper_list<-two_step_anpoll_mapper(df_anpoll_preprocessed = df_anpoll_preprocessed,output_path = file_path)

## Calculating the path distances of all streams inside "stream_list_modified"
stream_path_dist_vec<-stream_path_dist_cal(output_path = file_path)

## Loading the dataframe "df_anpoll_processed"
load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))

## Getting the 2nd MOST IMPORTANT dataframe and other lists for analyte
df_list_analyte_obj<-df_list_analyte_generator(df_analyte_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="analyte"),],output_path = file_path)

## Getting the 3rd MOST IMPORTANT dataframe for polluter
df_polluter_nodeID_aggregated<-df_polluter_generator(df_polluter_processed = df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="polluter"),c("nodeID","lon_mapped","lat_mapped","date")],output_path = file_path)

########################################################################################################
## Subsetting the river net and getting 1) "igraph_river_whole" 2) "graph order" and 3) "igraph_river_decomposed_list"
river_net_subset_list<-river_net_subset(output_path = file_path)

## Defining graph_subset_indicator and rank_subgraph
graph_subset_indicator<-F
#rank_subgraph<-2

## Loading the df_analyte_nodeID_aggregated and df_polluter_nodeID_aggregated
load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))

## Getting the analyte_vertex_ids
analyte_vertex_ids<-vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_analyte_nodeID_aggregated,output_path = file_path,analyte=T,polluter=F,graph_subset=graph_subset_indicator)

## Getting the polluter_vertex_ids
polluter_vertex_ids<-vertex_IDs_generator(df_anpoll_nodeID_aggregated = df_polluter_nodeID_aggregated,output_path = file_path,analyte=F,polluter=T, graph_subset=graph_subset_indicator)

## Getting the projected_nodeIDs_list
#undebug(projected_nodeIDs_list_generator)
projected_nodeIDs_list<-projected_nodeIDs_list_generator(file_path=file_path,projected_threshold_dist_km = 50)

#debug(flow_dist_polluter_projected_cal)
flow_dist_polluter_projected_list<-flow_dist_polluter_projected_cal(file_path = file_path)
#load(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/River network project/Analysis for Spill Paper/old_results/v3/Cl_spill_whole/polluter_files/projected_nodeIDs_vec.RData")
#load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
#df_polluter_nodeID_aggregated

########################################################################################################
## Creating shortest path edgelist by defining the segments and epochs and running in parallel 6 segments in each epoch with foreach

## Loading the analyte and polluter vertex ids
load(file = paste0(file_path,"analyte_files/analyte_vertex_ids.RData"))
load(file = paste0(file_path,"polluter_files/polluter_vertex_ids.RData"))
## Loading the vector of projected node IDs
load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
## Loading the igraph object for whole river network
load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
## Getting the vertex IDs for the projected node IDs
projected_vertex_ids<-which(V(igraph_river_whole)$name%in%projected_nodeIDs_vec)

## Combining the analyte, polluter and projected vertex ids
analyte_polluter_projected_vertex_ids<-sort(unique(c(analyte_vertex_ids,polluter_vertex_ids,projected_vertex_ids)))

## Defining number of chunks to be parallelized, chunk size and indices inside each chunk
n_chunks<-6
chunk_size<-ceiling(length(analyte_polluter_projected_vertex_ids)/n_chunks)
indices_all<-1:length(analyte_polluter_projected_vertex_ids)
indices_chunk_list<-split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))

cl <- makeCluster(spec = n_chunks)
registerDoParallel(cl)
shortest_path_anpoll_edgelist_chunk<-foreach (index = 1:n_chunks)%dopar% {
  shortest_path_edgelist_creator_parallelized_2(n_chunks=n_chunks,file_path = file_path)
}
stopCluster(cl)

shortest_path_anpoll_edgelist<-unlist(x = shortest_path_anpoll_edgelist_chunk,recursive = F)

## Removing duplicate paths
shortest_path_anpoll_edgelist<-unique(shortest_path_anpoll_edgelist)

## Removing the length 1 paths
shortest_path_length_vec<-sapply(X = shortest_path_anpoll_edgelist,FUN = length)
if(length(which(shortest_path_length_vec==1))>=1){
  shortest_path_anpoll_edgelist[which(shortest_path_length_vec==1)]<-NULL
}

## saving this path edgelist
save(shortest_path_anpoll_edgelist,file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))

####################################################
## Creating the analyte polluter network
anpoll_network_list<-anpoll_network_creator(output_path=file_path)

head(anpoll_network_list[[1]])
########################################################################################################
############################################# Inference ################################################
########################################################################################################
graph_subset_indicator<-F
## Loading the dataframe "df_analyte_nodeID_aggregated"
load(file = paste0(file_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
## Loading the dataframe "df_polluter_nodeID_aggregated"
load(file = paste0(file_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))

## Loading the analyte-polluter network edgelist "anpoll_edgelist"
load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
## Loading the analyte-polluter network edgelist "anpoll_edgelist"
load(file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))
## Loading the "total_edgelist_character_modified"
load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))

## Loading the dataframe "df_anpoll_processed"
load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
## Loading the vector of projected node IDs
load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))

## Subsetting the df_polluter_processed from df_anpoll_processed
df_polluter_processed<-df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="polluter"),]
nrow(df_polluter_processed)

## Subsetting the dataframe "df_polluter_nodeID_aggregated" for the case when river network is subsetted geographically, calculating the flow distance for each polluter from and to all connecting nodes in the anpoll network and storing them in a list
if(graph_subset_indicator){
  df_polluter_nodeID_aggregated_rank_subgaph<-df_polluter_nodeID_aggregated[which(df_polluter_nodeID_aggregated$rank_subgraph),]
  row.names(df_polluter_nodeID_aggregated_rank_subgaph)<-NULL
  ## Getting the flow dist from list
  flow_dist_from_list<-list()
  cl <- makeCluster(6)
  registerDoParallel(cl)
  flow_dist_from_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
    wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
  }
  stopCluster(cl)
  names(flow_dist_from_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
  save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))
  print("from_list is done")
  ## Getting the flow dist to list
  flow_dist_to_list<-list()
  cl <- makeCluster(6)
  registerDoParallel(cl)
  flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated_rank_subgaph))%dopar% {
    wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated_rank_subgaph,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = F, to_indicator = T, file_path = file_path)
  }
  stopCluster(cl)
  names(flow_dist_to_list)<-df_polluter_nodeID_aggregated_rank_subgaph$nodeID
  save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
}else{
  ## Getting the flow dist from list for polluters
  flow_dist_from_list_polluters<-list()
  cl <- makeCluster(6)
  registerDoParallel(cl)
  flow_dist_from_list_polluters<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated))%dopar% {
    wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
  }
  stopCluster(cl)
  names(flow_dist_from_list_polluters)<-df_polluter_nodeID_aggregated$nodeID
  print("from_list_polluter is done")
  
  ## Getting the flow dist from list for projected node IDs of polluters
  df_projected_nodeIDs<-data.frame("nodeID"=projected_nodeIDs_vec,stringsAsFactors = F)
  flow_dist_from_list_projected<-list()
  cl <- makeCluster(6)
  registerDoParallel(cl)
  #flow_dist_from_list_projected<-foreach (index = 1:nrow(df_projected_nodeIDs))%dopar% {
  flow_dist_from_list_projected<-foreach (index = 1:12)%dopar% {
    wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
  }
  stopCluster(cl)
  
  
  undebug(wrapper_flow_dist_cal)
  ptm<-proc.time()
  index<<-2
  temp<- wrapper_flow_dist_cal(df_polluter=df_projected_nodeIDs,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = T, to_indicator = F, file_path = file_path)
  print(proc.time()-ptm)
  
  
  
  
  names(flow_dist_from_list_projected)<-projected_nodeIDs_vec
  print("from_list_projected is done")
  
  ## Appending the from list for polluters and projected
  flow_dist_from_list<-append(flow_dist_from_list_polluters,flow_dist_from_list_projected)
  ## Saving the flow distance from list
  save(flow_dist_from_list,file = paste0(file_path,"inference/flow_dist_from_list.RData"))
  
  ## Getting the flow dist to list
  flow_dist_to_list<-list()
  cl <- makeCluster(6)
  registerDoParallel(cl)
  flow_dist_to_list<-foreach (index = 1:nrow(df_polluter_nodeID_aggregated))%dopar% {
    wrapper_flow_dist_cal(df_polluter=df_polluter_nodeID_aggregated,anpoll_edgelist = anpoll_edgelist, shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist, total_edgelist_character_modified = total_edgelist_character_modified, from_indicator = F, to_indicator = T, file_path = file_path)
  }
  stopCluster(cl)
  names(flow_dist_to_list)<-df_polluter_nodeID_aggregated$nodeID
  save(flow_dist_to_list,file = paste0(file_path,"inference/flow_dist_to_list.RData"))
}

########################################################################################################
## Summarizing the overall downstream distances for each polluter
flow_dist_to_summary_list<-lapply(X = flow_dist_to_list, FUN = function(x){if (!is.na(x)) return (summary(x)) else return(NULL)})
flow_dist_to_summary_list[sapply(flow_dist_to_summary_list, is.null)]<-NULL
df_to_summary<-do.call(rbind,flow_dist_to_summary_list)
colMeans(df_to_summary)
mean(df_to_summary[,"Median"])

########################################################################################################
## Getting the final matrices and dataframes with p values
graph_subset_indicator<-F
## Loading the dataframe "df_anpoll_processed"
load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
## Loading the dataframe "df_polluter_processed"
load(file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))
## Defining the downstream_threshold_dist_km lower and upper
#downstream_threshold_dist_km<-matrix(c(0,10,10,50),2,2,byrow=TRUE)
df_threshold_dist_km<-data.frame("polluter_intersection"=numeric(),"upstream"=numeric(),"downstream_lower"=numeric(),"downstream_upper"=numeric())

df_threshold_dist_km[1,]<-c(5,5,0,10)
df_threshold_dist_km[2,]<-c(45,5,10,50)

## Reading the spill_data_processed to get the affected water body and county
spill_data_processed<-read.csv(file = "Box Sync/PSU/Summer 2018/Geoscience_Research/River network project/Analysis for Spill Paper/Spill_data_processed.csv",stringsAsFactors = F)
spill_data_processed$Date<-as.Date(spill_data_processed$Date,format ="%m/%d/%y")
## Loop to get test results for spills for downstream close and far
for (j in 1:nrow(df_threshold_dist_km)){
  polluter_test_matrix<-data.frame(matrix(NA,nrow(df_polluter_processed),15))
  upstream_downstream_obs_list<-list()
  ## Loop to get test results for spills in rank_subgraph
  for(i in 1:nrow(polluter_test_matrix)){
    if(graph_subset_indicator){
      test_result<-polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"],date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
    }else{
      test_result<-polluter_test_dist_time(df_polluter = df_polluter_processed, polluter_lon = df_polluter_processed$lon_mapped[i], polluter_lat = df_polluter_processed$lat_mapped[i], polluter_projected_dist_km = df_threshold_dist_km[j,"polluter_intersection"],upstream_threshold_dist_km = df_threshold_dist_km[j,"upstream"], downstream_threshold_lower_dist_km = df_threshold_dist_km[j,"downstream_lower"], downstream_threshold_upper_dist_km = df_threshold_dist_km[j,"downstream_upper"], date_start=as.Date("1920-01-01"), date_end=as.Date("2017-12-31"), spill_date = df_polluter_processed$date[i], file_path = file_path)
    }
    print(i)
    #assign("last.warning", NULL, envir = baseenv())
    ## Getting mean, median and no. of observations for upstream and downstream in the first 6 columns
    polluter_test_matrix[i,c(1:6)]<-test_result[[3]][c(1:6)]
    ## Storing the p values and test indicator for upstream temporal test
    if(as.numeric(test_result[[1]][3])==1){
      polluter_test_matrix[i,7]<-as.numeric(test_result[[1]][1])
      polluter_test_matrix[i,8]<-as.numeric(test_result[[1]][2])
    }
    polluter_test_matrix[i,9]<-as.numeric(test_result[[1]][3])
    ## Storing the p values and test indicator for downstream temporal test
    if(as.numeric(test_result[[2]][3])==1){
      polluter_test_matrix[i,10]<-as.numeric(test_result[[2]][1])
      polluter_test_matrix[i,11]<-as.numeric(test_result[[2]][2])
    }
    polluter_test_matrix[i,12]<-as.numeric(test_result[[2]][3])
    ## Storing the p values and test indicator for final upstream vs. downstream spatio-temporal test
    if(as.numeric(test_result[[3]][9])==1){
      polluter_test_matrix[i,13]<-as.numeric(test_result[[3]][7])
      polluter_test_matrix[i,14]<-as.numeric(test_result[[3]][8])
    }
    polluter_test_matrix[i,15]<-as.numeric(test_result[[3]][9])
    ## Initializing and storing the observations for upstream and downstream
    upstream_downstream_obs_list[[i]]<-list()
    upstream_downstream_obs_list[[i]][[1]]<-test_result[[4]]
    upstream_downstream_obs_list[[i]][[2]]<-test_result[[5]]
    #print(i)
  }
  
  colnames(polluter_test_matrix)<-c("upstream mean","downstream mean","upstream median","downstream median","number of obs upstream","number of obs downstream", "t test (up) p value", "Wilcoxon test (up) p value","test result indicator (up)","t test (down) p value", "Wilcoxon test (down) p value","test result indicator (down)","t test (updown) p value", "Wilcoxon test (updown) p value","test result indicator (updown)")
  save(polluter_test_matrix,file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
  
  #####################################################
  ## Doing the fdr analysis for >=15 pvalues in eadch test
  load(file = paste0(file_path,"inference/polluter_test_matrix_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
  df_polluter_test<-fdr_analysis_wrapper(polluter_test_matrix = polluter_test_matrix,alpha = 0.1,county=spill_data_processed[,1],water_body=spill_data_processed[,2],file_path = file_path)
  df_polluter_test_mean<-df_polluter_test[[1]]
  df_polluter_test_median<-df_polluter_test[[2]]
  save(df_polluter_test_mean,file = paste0(file_path,"inference/df_polluter_test_mean_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
  save(df_polluter_test_median,file = paste0(file_path,"inference/df_polluter_test_median_",df_threshold_dist_km[j,"downstream_upper"],".RData"))
}
########################################################################################################