########################################################################################################
################ FINAL MODULARIZED CODE FOR RIVER NETWORK POLLUTER DETECTION ANALYSIS ##################
########################################################################################################
######################## A series of modular functions (properly commented) ############################
########################################################################################################

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
# library(lubridate)

########################################################################################################
########################################################################################################
## Writing a function which takes inputs as shape file (output of the function readOGR) and directory path where output files should be stored at and gives output as 1)  1st MOST IMPORTANT Dataframe "df_node_latlong" 2) character edgleist as "total_edgelist_character" 3) list with path information as "stream_list"
shapefile_pre_processing<-function(shape_obj,output_path){
  ## Extracting the attribute data frame
  df_shape<-shape_obj@data
  ## Extracting the lines list
  list_lines<-shape_obj@lines
  n_streams<-nrow(df_shape)
  ## Getting the edgelist in character form
  total_edgelist_character<-cbind(as.character(df_shape[,1]),as.character(df_shape[,2]),df_shape[,"STRAHLER"])
  ## Getting the total_edgelist_character
  save(total_edgelist_character,file = paste0(output_path,"common_files/total_edgelist_character.RData"))
  ## Getting the latlong aray for initial and final nodes of all streams (used later in creating the dataframe df_node_latlong) and stream_list
  network_latlong<-array(NA_real_,dim=c(2,2,n_streams))
  stream_list<-list()
  for(i in 1:n_streams){
    stream_list[[i]]<-list_lines[[i]]@Lines[[1]]@coords
    network_latlong[1,,i]<-stream_list[[i]][1,]
    network_latlong[2,,i]<-stream_list[[i]][nrow(stream_list[[i]]),]
  }
  ## Saving the stream_list
  save(stream_list,file = paste0(output_path,"common_files/stream_list.RData"))
  ## Putting all latlong info in one dataframe
  df_network_latlong<-data.frame(t(apply(X = network_latlong,MARGIN = 3,FUN = function(x){
    latlong_vec<-c(x[1,],x[2,])
    return(latlong_vec)
  })))
  names(df_network_latlong)<-c("start_lon","start_lat","end_lon","end_lat")
  ## combining the latlong dataframe with original dataframe
  df_shape_latlong<-cbind(df_shape,df_network_latlong)
  ######################################################
  ## Getting the edgelist in numeric form
  total_edgelist_numeric<-cbind(as.numeric(as.character(df_shape[,1])),as.numeric(as.character(df_shape[,2])))
  ## Creating a data frame that consists of unique node IDs and corresponding latlong
  df_node_latlong_whole<-data.frame(cbind(c(total_edgelist_numeric[,1],total_edgelist_numeric[,2]),c(df_shape_latlong$start_lon,df_shape_latlong$end_lon),c(df_shape_latlong$start_lat,df_shape_latlong$end_lat)))
  names(df_node_latlong_whole)<-c("nodeID","lon","lat")
  ## Subsetting for unique nodes only
  df_node_latlong<-unique(df_node_latlong_whole)
  df_node_latlong$nodeID<-as.character(df_node_latlong$nodeID)
  ## Saving the dataframe df_node_latlong
  save(df_node_latlong,file = paste0(output_path,"common_files/df_node_latlong.RData"))
  return(list(df_node_latlong,total_edgelist_character,stream_list))
}

########################################################################################################
## Writing a function for a two step mapping procedure of all C-PP locations by mapping them to nearest nodes in the base network or dividing the streams and creating new nodes. The 1st main goal is to add these mapped node IDs to the dataframe "df_anpoll_preprocessed" and return that as first argument. The 2nd main goal of this function is to modify all the preprocessed river network files ("df_node_latlong", "total_edgelist_character" and "stream_list") based on this finalized mapping. Inputs are 1) preprocessed dataframe for anpoll "df_anpoll_preprocessed" and 2) directory path where output file should be stored at and outputs 1) df_anpoll_preprocessed with additional column node ID and 2) all outputs same as shapefile_pre_processing but all of them modified
two_step_anpoll_mapper<-function(df_anpoll_preprocessed,output_path){
  ## Loading the dataframe "df_node_latlong" from common files
  load(file = paste0(output_path,"common_files/df_node_latlong.RData"))
  ## Loading the "total_edgelist_character" from common files
  load(file = paste0(output_path,"common_files/total_edgelist_character.RData"))
  ## Loading the "stream_list" from common files
  load(file = paste0(output_path,"common_files/stream_list.RData"))
  ###################################################
  ## Defining a function to create node mid stream and update everything in the global environment
  mid_stream_node_creator<-function(){
    flag_mid<<-flag_mid+1
    #print(paste0("flag_mid = ",flag_mid))
    ## Creating new node correspoding to this mapping over stream
    new_node_ID<<-as.character(max(as.numeric(df_node_latlong_modified[,1]))+1)
    ## Updating new node counter
    new_node_counter<<-new_node_counter+1
    ## Appending the "df_anpoll_preprocessed" with node_ID, lon_mapped and lat_mapped
    df_anpoll_preprocessed[i,"nodeID"]<<-new_node_ID
    df_anpoll_preprocessed[i,c("lon_mapped","lat_mapped")]<<-as.numeric(dist_to_stream[which.min(dist_to_stream[,1]),c(2,3)])
    ## Appending the total character edgelist
    total_edgelist_character_modified[nrow(total_edgelist_character_modified)+1,]<<-c(total_edgelist_character_modified[stream_ID_anpoll_selected,1],new_node_ID,total_edgelist_character_modified[stream_ID_anpoll_selected,3])
    total_edgelist_character_modified[nrow(total_edgelist_character_modified)+1,]<<-c(new_node_ID,total_edgelist_character_modified[stream_ID_anpoll_selected,2],total_edgelist_character_modified[stream_ID_anpoll_selected,3])
    ## Deleting the mapped stream from total character edgelist
    total_edgelist_character_modified<<-total_edgelist_character_modified[-stream_ID_anpoll_selected,]
    row.names(total_edgelist_character_modified)<<-NULL
    ## Storing new node ID latlong in the dataframe df_node_latlong_modified
    df_node_latlong_modified[nrow(df_node_latlong_modified)+1,1]<<-new_node_ID
    df_node_latlong_modified[nrow(df_node_latlong_modified),c(2,3)]<<-as.numeric(dist_to_stream[which.min(dist_to_stream[,1]),c(2,3)])
    ## Appending the new stream_list with path information
    stream_list_modified[[length(stream_list_modified)+1]]<<-as.matrix(rbind(as.matrix(stream_sub_df[stream_selected_row_IDs[1]:stream_selected_row_IDs[min(order(dist_to_points)[c(1,2)])],c(1,2)]),as.numeric(df_node_latlong_modified[nrow(df_node_latlong_modified),c(2,3)])))
    row.names(stream_list_modified[[length(stream_list_modified)]])<<-NULL
    colnames(stream_list_modified[[length(stream_list_modified)]])<<-NULL
    stream_list_modified[[length(stream_list_modified)+1]]<<-as.matrix(rbind(as.numeric(df_node_latlong_modified[nrow(df_node_latlong_modified),c(2,3)]),as.matrix(stream_sub_df[stream_selected_row_IDs[max(order(dist_to_points)[c(1,2)])]:stream_selected_row_IDs[length(stream_selected_row_IDs)],c(1,2)])))
    row.names(stream_list_modified[[length(stream_list_modified)]])<<-NULL
    colnames(stream_list_modified[[length(stream_list_modified)]])<<-NULL
    ## Deleting the mapped stream path from stream list
    stream_list_modified[[stream_ID_anpoll_selected]]<<-NULL
  }
  ###################################################
  ## Defining a function to map to node end of the stream and just update df_anpoll_preprocessed in the global environment
  end_stream_node_mapper<-function(){
    ## Appending the "df_anpoll_preprocessed" with node_ID, lon_mapped and lat_mapped
    column_ID<<-if(order(dist_to_points)[1]==1) 1 else 2
    df_anpoll_preprocessed[i,"nodeID"]<<-as.matrix(total_edgelist_character_modified)[stream_ID_anpoll_selected,column_ID]
    temp_ID<<-if(order(dist_to_points)[1]==1) 1 else length(dist_to_points)
    df_anpoll_preprocessed[i,c("lon_mapped","lat_mapped")]<<-stream_sub_df[stream_selected_row_IDs[temp_ID],c(1,2)]
  }
  ###################################################
  ## Initializing a dataframe for storing edgelist correspoding to new streams created after mapping
  total_edgelist_character_modified<-data.frame(total_edgelist_character,stringsAsFactors = F)
  ## Initializing a new stream_list that contains path information (latlong) of the two new streams created for all sample sites after mapping
  stream_list_modified<-stream_list
  ## Initializing the dataframe "df_node_latlong_modified"
  df_node_latlong_modified<-df_node_latlong
  ## Initializing the nodeID column for finally mapped nodes in the dataframe "df_anpoll_preprocessed"
  df_anpoll_preprocessed$nodeID<-NA
  df_anpoll_preprocessed$lon_mapped<-NA
  df_anpoll_preprocessed$lat_mapped<-NA
  
  ## Initializing the vector dist_mapped to calculate distances between original latlong of C-PP locations and mapped latlong on river streams
  dist_mapped<-rep(NA_real_,nrow(df_anpoll_preprocessed))
  ## Initializing the new node counter to count the number of new nodes created
  new_node_counter<-0
  ## Initializing the flag to count how many times a stream was broken
  flag_mid<-0
  ## Two step mapping loop for each row in df_anpoll_preprocessed; 1st step: Get top ten nearest node IDs, 2nd step: Create new nodes and streams by mapping to the river streams corresponding to top ten nearest node IDs
  for(i in 1:nrow(df_anpoll_preprocessed)){
    dist_vec_step1<-distm(x = as.matrix(df_anpoll_preprocessed[i,c("lon","lat")]),y = as.matrix(df_node_latlong_modified[,c("lon","lat")]))
    ## Getting the top ten nearest node IDs for sampling site i 
    top_ten_node_IDs_sampling_site_i<-df_node_latlong_modified$nodeID[order(as.vector(dist_vec_step1))[1:10]]
    ## Getting the stream IDs in the original edgelist which orginate or end in these top ten node IDs
    stream_IDs_sampling_site_i<-sort(unique(c(which(total_edgelist_character_modified[,1]%in%top_ten_node_IDs_sampling_site_i),which(total_edgelist_character_modified[,2]%in%top_ten_node_IDs_sampling_site_i))))
    ## Getting a smaller list of streams corresponding to these stream_IDs where each element stores the path information (sequence of lat long) of that stream 
    stream_sub_list<-stream_list_modified[stream_IDs_sampling_site_i]
    dist_to_stream<-matrix(NA_real_,nrow = length(stream_sub_list),3)
    for (j in 1:length(stream_sub_list)){
      stream_sub_list[[j]]<-data.frame(stream_sub_list[[j]])
      stream_sub_list[[j]][,3]<-j
      stream_sub_list[[j]][,4]<-0 ## Initializing start and end flag of stream
      stream_sub_list[[j]][1,4]<-1 ## Marking the start flag of stream as 1
      stream_sub_list[[j]][nrow(stream_sub_list[[j]]),4]<-2 ## Marking the start flag of stream as 2
      dist_to_stream[j,]<-dist2Line(p = as.matrix(df_anpoll_preprocessed[i,c("lon","lat")]),line = stream_sub_list[[j]][,c(1,2)])
    }
    ## Combining the path information for all these streams in one big dataframe with 3rd column depicting the local stream_ID
    stream_sub_df<-do.call(what = rbind,args = stream_sub_list)
    ## stream ID selected corresponding to closest stream from i^th sampling site
    stream_ID_anpoll_selected<-stream_IDs_sampling_site_i[which.min(dist_to_stream[,1])]
    ## Row ids in stream_sub_df which corresponds to selected stream
    stream_selected_row_IDs<-which(stream_sub_df[,3]==which.min(dist_to_stream[,1]))
    ## Calculating the distance of the i^th C-PP location to all the path latlong in this dataframe to judge whether the shortest distance to the selected stream occurs at midpoint or one of the endpoints
    dist_to_points<-as.vector(distm(x = as.matrix(df_anpoll_preprocessed[i,c("lon","lat")]), y = as.matrix(stream_sub_df[stream_selected_row_IDs,c(1,2)])))
    if(order(dist_to_points)[1]%in%c(1,length(dist_to_points))){
      if(order(dist_to_points)[1]==1){
        dist_12<-as.numeric(distm(x = as.matrix(stream_sub_df[stream_selected_row_IDs[1],c(1,2)]), y = as.matrix(stream_sub_df[stream_selected_row_IDs[2],c(1,2)])))
        angle_CPP_1_2<-acos(x = ((dist_to_points[1]^2+dist_12^2-dist_to_points[2]^2)/(2*dist_to_points[1]*dist_12)))*(180/pi)
        coincide_condition<-prod((as.numeric(stream_sub_df[stream_selected_row_IDs[1],c(1,2)])==as.numeric(dist_to_stream[which.min(dist_to_stream[,1]),c(2,3)])))
        if(!is.nan(angle_CPP_1_2)){
          if((angle_CPP_1_2>=90)|(coincide_condition)){
            end_stream_node_mapper()
          }else{
            mid_stream_node_creator()
          }
        }else{
          if(coincide_condition){
            end_stream_node_mapper()
          }else{
            mid_stream_node_creator()
          }
        }
      }else if(order(dist_to_points)[1]==length(dist_to_points)){
        dist_end_12<-as.numeric(distm(x = as.matrix(stream_sub_df[stream_selected_row_IDs[length(stream_selected_row_IDs)-1],c(1,2)]), y = as.matrix(stream_sub_df[stream_selected_row_IDs[length(stream_selected_row_IDs)],c(1,2)])))
        angle_CPP_end_1_2<-acos(x = ((dist_to_points[length(dist_to_points)]^2+dist_end_12^2-dist_to_points[length(dist_to_points)-1]^2)/(2*dist_to_points[length(dist_to_points)]*dist_end_12)))*(180/pi)
        coincide_condition<-prod((as.numeric(stream_sub_df[stream_selected_row_IDs[length(stream_selected_row_IDs)],c(1,2)])==as.numeric(dist_to_stream[which.min(dist_to_stream[,1]),c(2,3)])))
        if(!is.nan(angle_CPP_end_1_2)){
          if((angle_CPP_end_1_2>=90)|(coincide_condition)){
            end_stream_node_mapper()
          }else{
            mid_stream_node_creator()
          }
        }else{
          if(coincide_condition){
            end_stream_node_mapper()
          }else{
            mid_stream_node_creator()
          }
        }
      }
    }else{
      mid_stream_node_creator()
    }
    if((i%%1000)==0){
      print(i)
      #river_overlay_map<-river_overlay_map %>% fitBounds(~min_lon, ~min_lat, ~max_lon, ~max_lat) %>% addCircleMarkers(lng = df_anpoll_preprocessed[i,"lon"],lat = df_anpoll_preprocessed[i,"lat"], radius = 6,color = "red") %>%addCircleMarkers(lng = stream_sub_df[stream_selected_row_IDs[1],1],lat = stream_sub_df[stream_selected_row_IDs[1],2], radius = 4,color = "green")%>% addCircleMarkers(lng = df_anpoll_preprocessed[i,"lon_mapped"],lat = df_anpoll_preprocessed[i,"lat_mapped"], radius = 6,color = "black")
    }
    ##%>%addCircleMarkers(lng = dist_to_stream[,2],lat = dist_to_stream[,3], radius = 4,color = "yellow")
    dist_mapped[i]<-distm(x = df_anpoll_preprocessed[i,c("lon","lat")],y = df_anpoll_preprocessed[i,c("lon_mapped","lat_mapped")])
    #print(i)
  }
  
  ## Saving the mapped distances vector
  save(dist_mapped,file=paste0(output_path,"anpoll_files/dist_mapped.RData"))
  ## Renaming the modified dataframe of anpoll
  df_anpoll_processed<-df_anpoll_preprocessed
  row.names(df_anpoll_processed)<-NULL
  ## Saving the dataframe "df_anpoll_processed"
  save(df_anpoll_processed,file=paste0(output_path,"anpoll_files/df_anpoll_processed.RData"))
  ## Subsetting the df_polluter_processed from df_anpoll_processed
  df_polluter_processed<-df_anpoll_processed[which(df_anpoll_processed$anpoll_indicator=="polluter"),]
  row.names(df_polluter_processed)<-NULL
  nrow(df_polluter_processed)
  save(df_polluter_processed,file=paste0(file_path,"polluter_files/df_polluter_processed.RData"))
  ###################################################
  ## Saving the modified 1st MOST IMPORTANT dataframe df_node_latlong, modified total_character_edgelist and modified stream_list
  save(df_node_latlong_modified,file=paste0(output_path,"common_files_modified/df_node_latlong_modified.RData"))
  total_edgelist_character_modified<-as.matrix(total_edgelist_character_modified)
  save(total_edgelist_character_modified,file=paste0(output_path,"common_files_modified/total_edgelist_character_modified.RData"))
  save(stream_list_modified,file=paste0(output_path,"common_files_modified/stream_list_modified.RData"))
  return(list(df_anpoll_processed,df_node_latlong_modified,total_edgelist_character_modified,stream_list_modified))
}

########################################################################################################
## Writing a function to calculate the path distances of all streams in the modified common files. Input is directory path where output file should be stored at and output is path distance vector
stream_path_dist_cal<-function(output_path){
  ## Loading the "total_edgelist_character_modified" from common files modified
  load(file = paste0(output_path,"common_files_modified/total_edgelist_character_modified.RData"))
  ## Loading the "stream_list_modified" from common files modified
  load(file = paste0(output_path,"common_files_modified/stream_list_modified.RData"))
  stream_path_dist_vec<-rep(0,nrow(total_edgelist_character_modified))
  for(i in 1:nrow(total_edgelist_character_modified)){
    for (j in 1:(nrow(stream_list_modified[[i]])-1)){
      stream_path_dist_vec[i]<-stream_path_dist_vec[i]+as.numeric(distm(x = stream_list_modified[[i]][j,],y = stream_list_modified[[i]][j+1,]))
    }
    if((i%%1000)==0){
      print(i)
    }
  }
  save(stream_path_dist_vec,file = paste0(output_path,"common_files_modified/stream_path_dist_vec.RData"))
  return(stream_path_dist_vec)
}

########################################################################################################
########################################################################################################
## Writing a function which takes inputs as preprocessed dataframe for analyte "df_analyte_preprocessed" and directory path where output file should be stored at and gives output as the 2nd MOST IMPORTANT dataframe
df_list_analyte_generator<-function(df_analyte_processed,output_path){
  ## Taking mean of observations with duplicate node IDs
  df_analyte_nodeID_aggregated<-aggregate(df_analyte_processed[,c("lon_mapped","lat_mapped","conc")], by=list(df_analyte_processed$nodeID),FUN=mean, na.rm=TRUE)
  names(df_analyte_nodeID_aggregated)[1]<-"nodeID"
  ## Saving the dataframe "df_analyte_nodeID_aggregated"
  save(df_analyte_nodeID_aggregated,file = paste0(output_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
  ## List with each element corresponding to dataframe of all time information of ne node ID
  list_analyte_time<-split(df_analyte_processed, as.factor(df_analyte_processed$nodeID))
  list_analyte_time<-lapply(X = list_analyte_time,FUN = function(x){
    row.names(x)<-NULL
    return(x)
  })
  listID_nodeID_matrix<-matrix(NA,length(list_analyte_time),2)
  for(i in 1:length(list_analyte_time)){
    listID_nodeID_matrix[i,1]<-i
    listID_nodeID_matrix[i,2]<-list_analyte_time[[i]][1,"nodeID"]
  }
  save(list_analyte_time,file = paste0(output_path,"analyte_files/list_analyte_time.RData"))
  save(listID_nodeID_matrix,file = paste0(output_path,"analyte_files/listID_nodeID_matrix.RData"))
  return(list(df_analyte_nodeID_aggregated,list_analyte_time,listID_nodeID_matrix))
}

########################################################################################################
## Writing a function which takes inputs as preprocessed dataframe for polluter "df_polluter_preprocessed" and directory path where output file should be stored at and gives output as the 3rd MOST IMPORTANT dataframe
df_polluter_generator<-function(df_polluter_processed,output_path){
  ## Taking mean of observations with duplicate node IDs
  df_polluter_nodeID_aggregated<-aggregate(df_polluter_processed[,c("lon_mapped","lat_mapped")], by=list(df_polluter_processed$nodeID),FUN=mean, na.rm=TRUE)
  names(df_polluter_nodeID_aggregated)[1]<-"nodeID"
  names(df_polluter_nodeID_aggregated)<-c("nodeID","lon_mapped","lat_mapped")
  ## Saving the dataframe "df_polluter_nodeID_aggregated"
  save(df_polluter_nodeID_aggregated,file = paste0(output_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
  return(df_polluter_nodeID_aggregated)
}

########################################################################################################
########################################################################################################
## Writing a function that must be used after analyte and polluter mapping is done. Inputs are total_edgelist_character and directory path where output files should be stored at and gives output as 1) igraph object for whole network "igraph_river_whole" 2) "graph order" giving the vector of indices with decreasing size of most connected clusters 3) "igraph_river_decomposed_list" as a list of igraph objects corresponding to connected subgraphs
river_net_subset<-function(output_path){
  ## Loading the "total_edgelist_character_modified"
  load(file = paste0(output_path,"common_files_modified/total_edgelist_character_modified.RData"))
  ## Creating the igraph object from "total_edgelist_character_modified"
  igraph_river_whole<-graph_from_edgelist(total_edgelist_character_modified[,c(1,2)], directed = TRUE)
  save(igraph_river_whole,file = paste0(output_path,"common_files_modified/igraph_river_whole.RData"))
  ## Finding connected clusters in whole graph and getting "graph_order"
  connected_clusters<-clusters(igraph_river_whole)
  graph_order<-order(connected_clusters$csize,decreasing = TRUE)
  save(graph_order,file = paste0(output_path,"common_files_modified/graph_order.RData"))
  ## Decomposing graph into connected subgraphs
  igraph_river_decomposed_list <- decompose.graph(igraph_river_whole)
  save(igraph_river_decomposed_list,file = paste0(output_path,"common_files_modified/igraph_river_decomposed_list.RData"))
  return(list(igraph_river_whole,graph_order,igraph_river_decomposed_list))
}

########################################################################################################
########################################################################################################
## Writing a function which takes inputs as  1) analyte/polluter nodeID aggregated dataframe "df_anpoll_nodeID_aggregated 2) directory path where output files should be stored at and 3) rank_subgraph and gives output as vertex ids corresponding to analyte/polluter for the rank_subgraph
vertex_IDs_generator<-function(df_anpoll_nodeID_aggregated,output_path,analyte=T,polluter=F,graph_subset=T,rank_subgraph=2){
  ## Loading 1) igraph object for whole network "igraph_river_whole" 2) "graph order" giving the vector of indices with decreasing size of most connected clusters 3) "igraph_river_decomposed_list" as a list of igraph objects corresponding to connected subgraphs
  load(file = paste0(output_path,"common_files_modified/igraph_river_whole.RData"))
  load(file = paste0(output_path,"common_files_modified/graph_order.RData"))
  load(file = paste0(output_path,"common_files_modified/igraph_river_decomposed_list.RData"))
  ## Adding Observation vertex attributes in rank_subgraph
  ## checking how many observation node IDs overlap over overall nodeIDs in whole graph
  sum(V(igraph_river_whole)$name%in%df_anpoll_nodeID_aggregated$nodeID)
  ## checking how many observation node IDs overlap over overall nodeIDs in rank_subgraph biggest subgraph
  sum(V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$name%in%df_anpoll_nodeID_aggregated$nodeID)
  if(graph_subset){
    if(analyte){
      ## setting vertex attributes as observation being true/false
      igraph_river_decomposed_list[[graph_order[rank_subgraph]]]<-set_vertex_attr(graph = igraph_river_decomposed_list[[graph_order[rank_subgraph]]],name = "analyte",index = V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]]),value = V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$name%in%df_anpoll_nodeID_aggregated$nodeID)
      ## extracting vertex ids corresponding to observation locations
      analyte_vertex_ids<-which(vertex_attr(graph = igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$`analyte`==TRUE)
      ## saving the vertex ids in rank subgraph corresponding to analyte observations unique sites
      save(analyte_vertex_ids,file = paste0(output_path,"analyte_files/analyte_vertex_ids.RData"))
      anpoll_vertex_ids<-analyte_vertex_ids
      ## analyte node IDs for the rank subgraph
      analyte_nodeIDs_rank_subgraph<-V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$name[analyte_vertex_ids]
      ## Appending the "df_anpoll_nodeID_aggregated" with a column rank_subgraph_indicator that is TRUE or FALSE depending on whether the corresponding nodes are in rank_subgraph
      df_anpoll_nodeID_aggregated$rank_subgraph_indicator<-df_anpoll_nodeID_aggregated$nodeID%in%analyte_nodeIDs_rank_subgraph
      df_analyte_nodeID_aggregated<-df_anpoll_nodeID_aggregated
      ## Saving the updated dataframe "df_analyte_nodeID_aggregated"
      save(df_analyte_nodeID_aggregated,file = paste0(output_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
    }else if(polluter){
      ## setting vertex attributes as observation being true/false
      igraph_river_decomposed_list[[graph_order[rank_subgraph]]]<-set_vertex_attr(graph = igraph_river_decomposed_list[[graph_order[rank_subgraph]]],name = "polluter",index = V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]]),value = V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$name%in%df_anpoll_nodeID_aggregated$nodeID)
      ## extracting vertex ids corresponding to observation locations
      polluter_vertex_ids<-which(vertex_attr(graph = igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$`polluter`==TRUE)
      ## saving the vertex ids in rank subgraph corresponding to polluter observations unique sites
      save(polluter_vertex_ids,file = paste0(output_path,"polluter_files/polluter_vertex_ids.RData"))
      anpoll_vertex_ids<-polluter_vertex_ids
      ## polluter node IDs for the rank subgraph
      polluter_nodeIDs_rank_subgraph<-V(igraph_river_decomposed_list[[graph_order[rank_subgraph]]])$name[polluter_vertex_ids]
      ## Appending the "df_anpoll_nodeID_aggregated" with a column rank_subgraph_indicator that is TRUE or FALSE depending on whether the corresponding nodes are in rank_subgraph
      df_anpoll_nodeID_aggregated$rank_subgraph_indicator<-df_anpoll_nodeID_aggregated$nodeID%in%polluter_nodeIDs_rank_subgraph
      df_polluter_nodeID_aggregated<-df_anpoll_nodeID_aggregated
      ## Saving the updated dataframe "df_anpoll_nodeID_aggregated"
      save(df_polluter_nodeID_aggregated,file = paste0(output_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
    }
  }else{
    if(analyte){
      ## setting vertex attributes as observation being true/false
      igraph_river_whole<-set_vertex_attr(graph = igraph_river_whole,name = "analyte",index = V(igraph_river_whole),value = V(igraph_river_whole)$name%in%df_anpoll_nodeID_aggregated$nodeID)
      ## extracting vertex ids corresponding to observation locations
      analyte_vertex_ids<-which(vertex_attr(graph = igraph_river_whole)$`analyte`==TRUE)
      ## saving the vertex ids in rank subgraph corresponding to analyte observations unique sites
      save(analyte_vertex_ids,file = paste0(output_path,"analyte_files/analyte_vertex_ids.RData"))
      anpoll_vertex_ids<-analyte_vertex_ids
    }else if(polluter){
      ## setting vertex attributes as observation being true/false
      igraph_river_whole<-set_vertex_attr(graph = igraph_river_whole,name = "polluter",index = V(igraph_river_whole),value = V(igraph_river_whole)$name%in%df_anpoll_nodeID_aggregated$nodeID)
      ## extracting vertex ids corresponding to observation locations
      polluter_vertex_ids<-which(vertex_attr(graph = igraph_river_whole)$`polluter`==TRUE)
      ## saving the vertex ids in rank subgraph corresponding to polluter observations unique sites
      save(polluter_vertex_ids,file = paste0(output_path,"polluter_files/polluter_vertex_ids.RData"))
      anpoll_vertex_ids<-polluter_vertex_ids
    }
  }
  return(anpoll_vertex_ids)
}

########################################################################################################
########################################################################################################
## Writing a function to extract a list of projected node IDs in the base river network for each polluter node
projected_nodeIDs_list_generator<-function(file_path,projected_threshold_dist_km){
  ## Loading the dataframe df_polluter_nodeID_aggregated
  load(file = paste0(file_path,"anpoll_files/df_anpoll_processed.RData"))
  ## Loading the dataframe df_polluter_processed
  load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
  ## Loading the original base river network dataframe df_node_latlong
  load(file = paste0(file_path,"common_files/df_node_latlong.RData"))
  ## Loading the original base river network edgelist total_edgelist_character
  load(file = paste0(file_path,"common_files/total_edgelist_character.RData"))
  ## Loading the modified base river network edgelist total_edgelist_character_modified
  load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))
  ## Loading igraph object for whole network "igraph_river_whole" 
  load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
  ## Initializing the projected node IDs list
  projected_nodeIDs_list<-list()
  ## Loop to get the node ids of top 1000 nearest nodes, getting the corresponding ids and then calculating shortest path from original node to each of these vertex ids and then selecting those which have a higher Strahler number and storing them in the projected nodeIDs list
  #vertex_IDs_all<-which(V(igraph_river_whole)$name%in%df_node_latlong$nodeID)
  for(i in 1:nrow(df_polluter_processed)){
    ## Getting the node ids of top 1000 nearest nodes
    dist_vec<-distm(x = as.matrix(df_polluter_processed[i,c("lon","lat")]),y = as.matrix(df_node_latlong[,c("lon","lat")]))
    ## Getting the node IDs which lie within the projected_distance_threshold_km 
    node_IDs_site_i<-df_node_latlong$nodeID[which(dist_vec<=(projected_threshold_dist_km*1000))]
    
    ## Extracting vertex ids corresponding to these subsetted nodeIDs
    vertex_IDs_site_i<-which(V(igraph_river_whole)$name%in%node_IDs_site_i)
    #vertex_IDs_site_i<-vertex_IDs_all
    ## Getting the shortest path from polluter_i to all the subsetted nodes
    #ptm<-proc.time()
    short_path<-shortest_paths(graph = igraph_river_whole,from=which(V(igraph_river_whole)$name%in%df_polluter_processed$nodeID[i]), to= vertex_IDs_site_i,mode = "out",output = "vpath")
    #print(proc.time()-ptm)
    ## Trimming the paths with length 0 and storing all paths (among nodeIDs) in a list
    short_path_list<-lapply(X = short_path$vpath,FUN = function(x){
      if(length(as_ids(x))>0){
        y<-as_ids(x)
        return(y)
      }
    })
    short_path_list[sapply(short_path_list, is.null)] <- NULL
    ## Removing duplicate paths
    short_path_list<-unique(short_path_list)
    ## Removing the length 1 paths
    short_path_length_vec<-sapply(X = short_path_list,FUN = length)
    if(length(which(short_path_length_vec==1))>=1){
      short_path_list[which(short_path_length_vec==1)]<-NULL
    }
    
    
    ## Deleting all nodes that occur on same stream (same Strahler number for the stream downstream) but are farther apart.
    ## Initializing the list IDs to set NULL later
    # delete_IDs<-c()
    # if(length(short_path_list)>1){
    #   for (j in 1:length(short_path_list)){
    #     for (k in 1:length(short_path_list)){
    #       if(j!=k){
    #         if(short_path_list[[j]][length(short_path_list[[j]])]%in%short_path_list[[k]][-c(1,length(short_path_list[[k]]))]){
    #           if(length(which(total_edgelist_character[,2]==short_path_list[[k]][length(short_path_list[[k]])]))==1){
    #             delete_IDs<-c(delete_IDs,k)
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    # short_path_list[sort(unique(delete_IDs))]<-NULL
    ## Extracting the end node_IDs and declaring them as projected node IDs
    projected_nodeIDs_list[[i]]<-sapply(X = short_path_list,FUN = function(x) return(x[length(x)]))
    print(i)
  }
  names(projected_nodeIDs_list)<-df_polluter_processed$nodeID
  ## Saving the projected nodeIDs list
  save(projected_nodeIDs_list,file = paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))
  ## Extracting the vector of projected node_IDs
  projected_nodeIDs_vec<-unlist(projected_nodeIDs_list)
  attr(projected_nodeIDs_vec,"names")<-NULL
  ## Saving the projected nodeIDs vector
  save(projected_nodeIDs_vec,file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
  return(list(projected_nodeIDs_list,projected_nodeIDs_vec))
}

#######################################################################################################
########################################################################################################
## Defining a function to calculate flow distance of all polluter nodes to their projected nodes
flow_dist_polluter_projected_cal<-function(file_path){
  ## Loading the projected node IDs list
  load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
  ## Loading the total_edgelist_character_modified  
  load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))
  ## Loading the total_edgelist_character_modified
  load(file = paste0(file_path,"common_files_modified/stream_path_dist_vec.RData"))
  ## Loading the projected node IDs list
  load(file=paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))
  
  ## Defining a function to calculate between a test node and a current node. The current nodes are usually polluter nodes.
  flow_dist_cal<-function(short_path_vec, test_node_ID, current_node_ID, file_path){
    ## Initializing and loop to fill up the matrix of edges that are formed by the nodeIDs in the path of test and current node ID
    flow_path_nodeIDs_mat<-matrix(NA_character_,nrow = (length(short_path_vec)-1),ncol = 2)
    for (k in 1:(length(short_path_vec)-1)){
      flow_path_nodeIDs_mat[k,]<-c(short_path_vec[k],short_path_vec[k+1])
    }
    ## Initializing the flow distance in meters
    flow_dist_m<-0
    ## Loop to add up the flow distance of all edges in the flow_path_nodeIDs
    for(k in 1:nrow(flow_path_nodeIDs_mat)){
      flow_dist_m<-flow_dist_m+stream_path_dist_vec[which((total_edgelist_character_modified[,1]==flow_path_nodeIDs_mat[k,1])&(total_edgelist_character_modified[,2]==flow_path_nodeIDs_mat[k,2]))]
    }
    ## flow distance in km
    flow_dist_km<-flow_dist_m/1000
    return(flow_dist_km)
  }
  
  ## Initializing the list to store flow distances from polluter to its projected nodes (intersections)
  flow_dist_polluter_projected_list<-list()
  for (i in 1:length(projected_nodeIDs_list)){
    current_node_ID<-names(projected_nodeIDs_list)[i]
    ## Initializing the flow_dist_vec as NA. If there are no connected nodes, this would remain NA
    flow_dist_vec<-NA
    ## if condition to check if there is atleast one projected node for ith polluter
    if(length(projected_nodeIDs_list[[i]])>0){
      ## Extracting vertex ids corresponding to these subsetted nodeIDs
      vertex_IDs_projected_i<-which(V(igraph_river_whole)$name%in%projected_nodeIDs_list[[i]])
      #vertex_IDs_site_i<-vertex_IDs_all
      ## Getting the shortest path from polluter_i to all the subsetted nodes
      #ptm<-proc.time()
      short_path<-shortest_paths(graph = igraph_river_whole,from=which(V(igraph_river_whole)$name%in%current_node_ID), to= vertex_IDs_projected_i,mode = "out",output = "vpath")
      #print(proc.time()-ptm)
      ## Trimming the paths with length 0 and storing all paths (among nodeIDs) in a list
      short_path_list<-lapply(X = short_path$vpath,FUN = function(x){
        if(length(as_ids(x))>0){
          y<-as_ids(x)
          return(y)
        }
      })
      ## Initializing the flow_dist_vec
      flow_dist_vec<-rep(NA_real_,length(projected_nodeIDs_list[[i]]))
      for(j in 1:length(projected_nodeIDs_list[[i]])){
        flow_dist_vec[j]<-flow_dist_cal(short_path_vec=short_path_list[[j]],test_node_ID = projected_nodeIDs_list[[i]][j],current_node_ID = current_node_ID,file_path = file_path)
        #print(j)
      }
    }
    flow_dist_polluter_projected_list[[i]]<-flow_dist_vec
    print(i)
  }
  ## saving the flow_dist_polluter_projected_list
  save(flow_dist_polluter_projected_list, file=paste0(file_path,"polluter_files/flow_dist_polluter_projected_list.RData"))
  return(flow_dist_polluter_projected_list)
}

########################################################################################################
########################################################################################################
## Writing a function to create shortest path edgelist and creating the network for analyte-polluter. Inputs are  1) analyte and polluter nodeID aggregated dataframes "df_analyte_nodeID_aggregated" and "df_polluter_nodeID_aggregated" 2) directory path where output files should be stored at. Output shortest path edgelist
shortest_path_edgelist_creator<-function(graph_subset=T,rank_subgraph=2,file_path){
  library(igraph)
  ## Loading the analyte and polluter vertex ids
  load(file = paste0(file_path,"analyte_files/analyte_vertex_ids.RData"))
  load(file = paste0(file_path,"polluter_files/polluter_vertex_ids.RData"))
  
  ## Loading the vector of projected node IDs
  load(file = paste0(file_path,"polluter_files/projected_nodeIDs_vec.RData"))
  ## Loading the igraph object for whole river network
  load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
  ## Getting the vertex IDs for the projected node IDs
  projected_vertex_ids<-which(V(igraph_river_whole)$name%in%projected_nodeIDs_vec)
  
  ## Loading 1) igraph object for whole network "igraph_river_whole" 2) "graph order" giving the vector of indices with decreasing size of most connected clusters 3) "igraph_river_decomposed_list" as a list of igraph objects corresponding to connected subgraphs
  #load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
  #load(file = paste0(file_path,"common_files_modified/graph_order.RData"))
  #load(file = paste0(file_path,"common_files_modified/igraph_river_decomposed_list.RData"))
  ## Combining the analyte and polluter vertex ids
  analyte_polluter_projected_vertex_ids<-sort(unique(c(analyte_vertex_ids,polluter_vertex_ids,projected_vertex_ids)))
  ## Initializing a list to store shortest paths of all vertex ids in the current segment corresponding to index
  edge_list_index<-list()
  ## Defining the igraph_obj used to create short paths
  # if(graph_subset){
  #   igraph_obj<-igraph_river_decomposed_list[[graph_order[rank_subgraph]]]
  # }else{
  #   igraph_obj<-igraph_river_whole
  # }
  ## Loop to create the shortest path objects and shortest_path_edgelist for the 20 vertex ids corresponding to given segment to all other vertex ids
  for(i in 1:length(analyte_polluter_projected_vertex_ids)){
    ptm<-proc.time()
    short_path<-shortest_paths(graph = igraph_river_whole,from=analyte_polluter_projected_vertex_ids[i], to= analyte_polluter_projected_vertex_ids[-i],mode = "out",output = "vpath")
    ## Trimming the paths with length 0 and storing all paths (among nodeIDs) in a list
    edge_list_index[[i]]<-lapply(X = short_path$vpath,FUN = function(x){
      if(length(as_ids(x))>0){
        y<-as_ids(x)
        return(y)
      }
    })
    edge_list_index[[i]][sapply(edge_list_index[[i]], is.null)] <- NULL
    print(i)
    print(proc.time()-ptm)
  }
  
  edge_list_index<-unlist(x = edge_list_index,recursive = F)
  return(edge_list_index)
}

## Writing a function (used later in foreach) to create shortest path edgelist and creating the network for analyte-polluter. Inputs are  1) analyte and polluter nodeID aggregated dataframes "df_analyte_nodeID_aggregated" and "df_polluter_nodeID_aggregated" 2) directory path where output files should be stored at. Output shortest path edgelist; This parallelization function has some problems on cluster
shortest_path_edgelist_creator_parallelized<-function(segment_index,file_path){
  library(igraph)
  ## Loading the analyte and polluter vertex ids
  load(file = paste0(file_path,"analyte_files/analyte_vertex_ids.RData"))
  load(file = paste0(file_path,"polluter_files/polluter_vertex_ids.RData"))
  ## Loading 1) igraph object for whole network "igraph_river_whole" 2) "graph order" giving the vector of indices with decreasing size of most connected clusters 3) "igraph_river_decomposed_list" as a list of igraph objects corresponding to connected subgraphs
  load(file = paste0(file_path,"common_files_modified/igraph_river_whole.RData"))
  #load(file = paste0(file_path,"common_files_modified/graph_order.RData"))
  #load(file = paste0(file_path,"common_files_modified/igraph_river_decomposed_list.RData"))
  ## Getting the vertex IDs for the projected node IDs
  projected_vertex_ids<-which(V(igraph_river_whole)$name%in%projected_nodeIDs_vec)
  ## Combining the analyte and polluter vertex ids
  analyte_polluter_projected_vertex_ids<-sort(unique(c(analyte_vertex_ids,polluter_vertex_ids,projected_vertex_ids)))
  ## Defining the segment length and vertex_ids_segments
  segment_length<-10
  segment_start_indices<-seq(from=1,to = length(analyte_polluter_projected_vertex_ids),by = segment_length)
  ## Initializing a list to store shortest paths of all vertex ids in the current segment corresponding to index
  edge_list_index<-list()
  ## Defining the igraph_obj used to create short paths
  # if(graph_subset){
  #   igraph_obj<-igraph_river_decomposed_list[[graph_order[rank_subgraph]]]
  # }else{
  #   igraph_obj<-igraph_river_whole
  # }
  ## Loop to create the shortest path objects and shortest_path_edgelist for the 20 vertex ids corresponding to given segment to all other vertex ids
  if(segment_index!=length(segment_start_indices)){
    for(i in 1:segment_length){
      short_path<-shortest_paths(graph = igraph_river_whole,from=analyte_polluter_projected_vertex_ids[segment_start_indices[segment_index]+(i-1)], to= analyte_polluter_projected_vertex_ids[-(segment_start_indices[segment_index]+(i-1))],mode = "out",output = "vpath")
      ## Trimming the paths with length 0 and storing all paths (among nodeIDs) in a list
      edge_list_index[[i]]<-lapply(X = short_path$vpath,FUN = function(x){
        if(length(as_ids(x))>0){
          y<-as_ids(x)
          return(y)
        }
      })
      edge_list_index[[i]][sapply(edge_list_index[[i]], is.null)] <- NULL
    }
  }else{
    if((length(analyte_polluter_projected_vertex_ids)%%segment_length)>0){
      for(i in 1:(length(analyte_polluter_projected_vertex_ids)%%segment_length)){
        short_path<-shortest_paths(graph = igraph_obj,from=analyte_polluter_projected_vertex_ids[segment_start_indices[segment_index]+(i-1)], to= analyte_polluter_projected_vertex_ids[-(segment_start_indices[segment_index]+(i-1))],mode = "out",output = "vpath")
        ## Trimming the paths with length 0 and storing all paths (among nodeIDs) in a list
        edge_list_index[[i]]<-lapply(X = short_path$vpath,FUN = function(x){
          if(length(as_ids(x))>0){
            y<-as_ids(x)
            return(y)
          }
        })
        edge_list_index[[i]][sapply(edge_list_index[[i]], is.null)] <- NULL
      }
    }
  }
  edge_list_index<-unlist(x = edge_list_index,recursive = F)
  return(edge_list_index)
}

## Writing a function (used later in foreach) to create shortest path edgelist and creating the network for analyte-polluter. Inputs are  1) analyte and polluter nodeID aggregated dataframes "df_analyte_nodeID_aggregated" and "df_polluter_nodeID_aggregated" 2) directory path where output files should be stored at. Output shortest path edgelist
shortest_path_edgelist_creator_parallelized_2<-function(n_chunks,file_path){
  library(igraph)
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
  save(analyte_polluter_projected_vertex_ids,file = paste0(file_path,"anpoll_files/analyte_polluter_projected_vertex_ids.RData"))
  
  chunk_size<-ceiling(length(analyte_polluter_projected_vertex_ids)/n_chunks)
  indices_all<-1:length(analyte_polluter_projected_vertex_ids)
  indices_chunk_list<-split(x = indices_all, ceiling(seq_along(indices_all)/chunk_size))
  
  current_indices<-indices_chunk_list[[index]]
  
  ## Initializing a list to store shortest paths of all vertex ids in the current segment corresponding to index
  edge_list_index<-list()
  ## Loop to create the shortest path objects and shortest_path_edgelist for the 20 vertex ids corresponding to given segment to all other vertex ids
  for(i in 1:length(current_indices)){
    #ptm<-proc.time()
    short_path<-shortest_paths(graph = igraph_river_whole,from=analyte_polluter_projected_vertex_ids[current_indices[i]], to= analyte_polluter_projected_vertex_ids[-current_indices[i]],mode = "out",output = "vpath")
    ## Trimming the paths with length 0 and storing all paths (among nodeIDs) in a list
    edge_list_index[[i]]<-lapply(X = short_path$vpath,FUN = function(x){
      if(length(as_ids(x))>0){
        y<-as_ids(x)
        return(y)
      }
    })
    edge_list_index[[i]][sapply(edge_list_index[[i]], is.null)] <- NULL
    #print(i)
    #print(proc.time()-ptm)
  }
  edge_list_index<-unlist(x = edge_list_index,recursive = F)
  save(edge_list_index, file = paste0(file_path,"anpoll_files/shortest_path_chunks/shortest_path_anpoll_edgelist_chunk_",index,".RData"))
  return(edge_list_index)
}

########################################################################################################
########################################################################################################
## Writing a function to create edgelist and node dataframe for analyte-polluter network using shortest_path_edgelist. Inputs is directory path where output files should be stored at. Outputs are 1) "anpoll_edgelist" and 2) df_node_latlong_anpoll
anpoll_network_creator<-function(output_path){
  ## Loading the path edgelist
  load(file = paste0(output_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))
  ## Loop to create analyte_polluter_edgelist from shortest path edgelist
  anpoll_edgelist<-matrix(NA_character_,nrow = length(shortest_path_anpoll_edgelist),ncol = 2)
  for (i in 1:length(shortest_path_anpoll_edgelist)){
    n_i<-length(shortest_path_anpoll_edgelist[[i]])
    anpoll_edgelist[i,]<-shortest_path_anpoll_edgelist[[i]][c(1,n_i)]
    if((i%%1000)==0){
      print(i)
    }
  }
  ## Saving the "anpoll_edgelist"
  save(anpoll_edgelist,file = paste0(output_path,"anpoll_files/anpoll_edgelist.RData"))
  ## Loading the dataframes "df_analyte_nodeID_aggregated" and "df_polluter_nodeID_aggregated"
  load(file = paste0(output_path,"analyte_files/df_analyte_nodeID_aggregated.RData"))
  load(file = paste0(output_path,"polluter_files/df_polluter_nodeID_aggregated.RData"))
  ## Combining the unique nodeIDs for analyte and polluter in a single dataframe
  df_node_latlong_anpoll<-data.frame(rbind(df_analyte_nodeID_aggregated[,c(1,2,3)],df_polluter_nodeID_aggregated[,c(1,2,3)]),"anpoll_indicator"=c(rep("analyte",nrow(df_analyte_nodeID_aggregated)),rep("polluter",nrow(df_polluter_nodeID_aggregated))),stringsAsFactors = F)
  rownames(df_node_latlong_anpoll)<-NULL
  save(df_node_latlong_anpoll,file = paste0(output_path,"anpoll_files/df_node_latlong_anpoll.RData"))
  return(list(anpoll_edgelist,df_node_latlong_anpoll))
}

#######################################################################################################
########################################################################################################
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
    for(j in 1:length(connected_nodeIDs)){
      flow_dist_vec[j]<-flow_dist_cal(anpoll_edgelist = anpoll_edgelist,shortest_path_anpoll_edgelist = shortest_path_anpoll_edgelist,total_edgelist_character_modified = total_edgelist_character_modified,test_node_ID = connected_nodeIDs[j],current_node_ID = df_polluter$nodeID[index],from_indicator = from_indicator,to_indicator = to_indicator,file_path = file_path)
      #print(j)
    }
  }
  return(flow_dist_vec)
}

########################################################################################################
############################################# Inference ################################################
########################################################################################################

## Writing a function to collect all from_node_IDs within the threshold flow distance with respect to given node ID
from_nodeIDs_crawler<-function(nodeID,upstream_threshold_dist_km,file_path){
  ## Loading the common edgelist total_edgelist_character_modified
  load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))
  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
  ## Loading the "flow_dist_from_list"
  load(file = paste0(file_path,"inference/flow_dist_from_list.RData"))
  ## Collecting all from nodes that point to the polluter_nodeID
  from_edgelist_ids_total<-which(anpoll_edgelist[,2]==nodeID)
  ## Intializing from_counter_flag that becomes 1 if we could find atleast one node within the specified threshold flow distance
  from_counter_flag<-0
  ## Intializing the vector of from_nodeIDs
  from_nodeIDs<-c()
  ## If loop to check if there is atleast one node that is directed towards the specified polluter
  if(length(from_edgelist_ids_total)>0){
    ## Extracting all nodeIDs that are directed towards polluter_node_ID
    from_nodeIDs_total<-anpoll_edgelist[from_edgelist_ids_total,1]
    ## Computing the flow distance indicator vector to indicate which of these nodeIDs are within the specified threshold flow distance
    flow_dist_indicator_vec<-(flow_dist_from_list[[which(names(flow_dist_from_list)==nodeID)[1]]])<=(rep(upstream_threshold_dist_km,length(flow_dist_from_list[[which(names(flow_dist_from_list)==nodeID)[1]]])))
    ## Indexes within the "from_nodeIDs_total" vector which correspond to those nodeIDs that are within the specified threshold flow distance
    flow_dist_within_ids<-which(flow_dist_indicator_vec)
    ## Extracting the from_nodeIDs corresponding to the nodes within the specified threshold flow distance and from_edgelist_ids corresponding to row IDs in the "anpoll_edgelist" that are within specified threshold flow distance. Also updating "from_counter_flag"
    if(length(flow_dist_within_ids)>0){
      from_nodeIDs<-from_nodeIDs_total[flow_dist_within_ids]
      from_counter_flag<-1
    }else{
      from_counter_flag<-0
    }
  }
  return(list(from_nodeIDs, from_counter_flag))
}

## Writing a function to do inference that is testing the significance of polluter for given polluter coordinates, threshold flow distance and time interval, to collect all observations upstream and downstream within the specified threshold flow distance and time and do a two sample two sided t test and Wilcoxon rank sum test. Inputs are 1) latlong of polluter 2) threshold flow distance and 3) time interval. Outputs are mean, median and length of observations upstream and downstream within the specified threshold flow distance and year interval and the corresponding p values
polluter_test_dist_time<-function(df_polluter, polluter_lon, polluter_lat, polluter_projected_dist_km, upstream_threshold_dist_km, downstream_threshold_lower_dist_km, downstream_threshold_upper_dist_km, date_start, date_end, spill_date, file_path){
  ####################################################
  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path,"anpoll_files/anpoll_edgelist.RData"))
  ## Loading the analyte-polluter network edgelist "anpoll_edgelist"
  load(file = paste0(file_path,"anpoll_files/shortest_path_anpoll_edgelist.RData"))
  ## Loading the "total_edgelist_character_modified"
  load(file = paste0(file_path,"common_files_modified/total_edgelist_character_modified.RData"))
  ## Loading the "flow_dist_from_list"
  load(file = paste0(file_path,"inference/flow_dist_from_list.RData"))
  ## Loading the "flow_dist_to_list"
  load(file = paste0(file_path,"inference/flow_dist_to_list.RData"))
  ## Loading the projected_nodeIDs_list
  load(file=paste0(file_path,"polluter_files/projected_nodeIDs_list.RData"))
  ## Loading the distances of polluters to projected nodes
  load(file=paste0(file_path,"polluter_files/flow_dist_polluter_projected_list.RData"))
  
  ## Extracting polluter node ID from the dataframe given the polluter latlong
  polluter_node_ID<-unique(df_polluter$nodeID[which((df_polluter$lon_mapped==polluter_lon)&(df_polluter$lat_mapped==polluter_lat))])
  
  ####################################################
  ## Collecting all to nodes that are directed from the polluter_nodeID
  to_edgelist_ids_total<-which(anpoll_edgelist[,1]==polluter_node_ID)
  ## Intializing to_counter_flag that becomes 1 if we could find atleast one node within the specified threshold flow distance
  to_counter_flag<-0
  if(length(to_edgelist_ids_total)>0){
    ## Extracting all nodeIDs that are directed from the polluter_node_ID
    to_nodeIDs_total<-anpoll_edgelist[to_edgelist_ids_total,2]
    ## Computing the flow distance indicator vector to indicate which of these nodeIDs are within the specified threshold flow distance
    flow_dist_indicator_vec<-((flow_dist_to_list[[which(names(flow_dist_to_list)==polluter_node_ID)]])<=(rep(downstream_threshold_upper_dist_km,length(flow_dist_to_list[[which(names(flow_dist_to_list)==polluter_node_ID)]]))))&((flow_dist_to_list[[which(names(flow_dist_to_list)==polluter_node_ID)]])>(rep(downstream_threshold_lower_dist_km,length(flow_dist_to_list[[which(names(flow_dist_to_list)==polluter_node_ID)]]))))
    ## Indexes within the "to_nodeIDs_total" vector which correspond to those nodeIDs that are within the specified threshold flow distance
    flow_dist_within_ids<-which(flow_dist_indicator_vec)
    ## Extracting the to_nodeIDs corresponding to the nodes within the specified threshold flow distance and to_edgelist_ids corresponding to row IDs in the "anpoll_edgelist" that are within specified threshold flow distance. Also updating "to_counter_flag"
    if(length(flow_dist_within_ids)>0){
      to_nodeIDs<-to_nodeIDs_total[flow_dist_within_ids]
      to_edgelist_ids<-to_edgelist_ids_total[flow_dist_within_ids]
      to_counter_flag<-1
    }else{
      to_counter_flag<-0
    }
  }
  to_nodeIDs<-unique(to_nodeIDs)
  
  ####################################################
  ## Collecting all from nodes that point to the polluter_nodeID and corresponding projected node_IDs
  projected_nodeIDs_total<-projected_nodeIDs_list[[which(names(projected_nodeIDs_list)==polluter_node_ID)[1]]]
  
  ## Subsetting the projected node IDs based on the distances of polluter to projected
  projected_vec_ids<-which(flow_dist_polluter_projected_list[[which(names(projected_nodeIDs_list)==polluter_node_ID)[1]]]<=rep(polluter_projected_dist_km,length(projected_nodeIDs_total)))
  projected_nodeIDs_subset<-projected_nodeIDs_total[projected_vec_ids]
  
  ## Initializing the vector of from_nodeIDs containing upstream sampling nodes for the polluter node and it's projected nodes
  from_nodeIDs<-c()
  ## Intializing the from_counter_flag to know later whether there is atleast one upstream node
  from_count<-0
  
  ## Getting the upstream nodes for polluter node ID 
  from_nodeIDs_polluter_list<-from_nodeIDs_crawler(nodeID = polluter_node_ID,upstream_threshold_dist_km = upstream_threshold_dist_km, file_path = file_path)
  from_nodeIDs<-c(from_nodeIDs,from_nodeIDs_polluter_list[[1]])
  from_count<-from_count+from_nodeIDs_polluter_list[[2]]
  
  ## Getting the upstream nodes for projected node IDs
  for(i in 1:length(projected_nodeIDs_subset)){
    from_nodeIDs_projected_list<-from_nodeIDs_crawler(nodeID = projected_nodeIDs_subset[i],upstream_threshold_dist_km = upstream_threshold_dist_km, file_path = file_path)
    if(length(to_edgelist_ids_total)>0){
      from_nodeIDs_subset<-setdiff(from_nodeIDs_projected_list[[1]],to_nodeIDs_total)
    }else{
      from_nodeIDs_subset<-from_nodeIDs_projected_list[[1]]
    }
    from_nodeIDs<-c(from_nodeIDs,from_nodeIDs_subset)
    from_count<-from_count+from_nodeIDs_projected_list[[2]]
  }
  from_nodeIDs<-unique(from_nodeIDs)
  
  ####################################################
  ## Collecting all from and to observations within the time interval specified
  
  ## Loading the dataframe "df_node_latlong_anpoll.RData"
  load(file = paste0(file_path,"anpoll_files/df_node_latlong_anpoll.RData"))
  ## Loading the  "listID_nodeID_matrix"
  load(file = paste0(file_path,"analyte_files/listID_nodeID_matrix.RData"))
  ## Loading the  "list_analyte_time"
  load(file = paste0(file_path,"analyte_files/list_analyte_time.RData"))
  
  ## Initializing the from observation vectors upstream of given polluter node ID
  from_obs_total<-c()
  from_obs_before_spill<-c()
  from_obs_after_spill<-c()
  
  ## Initializing the to observation vectors upstream of given polluter node ID
  to_obs_total<-c()
  to_obs_before_spill<-c()
  to_obs_after_spill<-c()
  
  ## Loop to collect all from obervations within the time interval specified
  if(from_count>=1){
    for (i in 1:length(from_nodeIDs)){
      analyte_indicator<- c("analyte")%in%(df_node_latlong_anpoll[which(df_node_latlong_anpoll$nodeID==from_nodeIDs[i]),"anpoll_indicator"])
      if(analyte_indicator){
        listID_i<-which(listID_nodeID_matrix[,2]==from_nodeIDs[i])
        ## Getting all upstream observations
        row_ids_date_interval_total<-which(between((list_analyte_time[[listID_i]]$date), date_start, date_end))
        if(length(row_ids_date_interval_total)>0){
          from_obs_total<-c(from_obs_total,list_analyte_time[[listID_i]]$conc[row_ids_date_interval_total])
        }
        ## Getting upstream observations before the spill date
        row_ids_date_interval_before_spill<-which(between((list_analyte_time[[listID_i]]$date), date_start, spill_date))
        if(length(row_ids_date_interval_before_spill)>0){
          from_obs_before_spill<-c(from_obs_before_spill,list_analyte_time[[listID_i]]$conc[row_ids_date_interval_before_spill])
        }
        ## Getting upstream observations till one year after the spill date
        row_ids_date_interval_after_spill<-which(between((list_analyte_time[[listID_i]]$date), spill_date, spill_date %m+% years(1)))
        if(length(row_ids_date_interval_after_spill)>0){
          from_obs_after_spill<-c(from_obs_after_spill,list_analyte_time[[listID_i]]$conc[row_ids_date_interval_after_spill])
        }
      }
    }
  }
  
  ## Loop to collect all to obervations within the time interval specified
  if(to_counter_flag==1){
    for (i in 1:length(to_nodeIDs)){
      analyte_indicator<-c("analyte")%in%(df_node_latlong_anpoll[which(df_node_latlong_anpoll$nodeID==to_nodeIDs[i]),"anpoll_indicator"])
      if(analyte_indicator){
        listID_i<-which(listID_nodeID_matrix[,2]==to_nodeIDs[i])
        ## Getting all downstream observations with after one year spill date condition
        row_ids_date_interval_total<-which((between((list_analyte_time[[listID_i]]$date), date_start, date_end))&(between((list_analyte_time[[listID_i]]$date), spill_date, spill_date %m+% years(1))))
        if(length(row_ids_date_interval_total)>0){
          to_obs_total<-c(to_obs_total,list_analyte_time[[listID_i]]$conc[row_ids_date_interval_total])
        }
        ## Getting all downstream observations before the spill date
        row_ids_date_interval_before_spill<-which(between((list_analyte_time[[listID_i]]$date), date_start, spill_date))
        if(length(row_ids_date_interval_before_spill)>0){
          to_obs_before_spill<-c(to_obs_before_spill,list_analyte_time[[listID_i]]$conc[row_ids_date_interval_before_spill])
        }
        ## Getting downstream observations till one year after the spill date
        row_ids_date_interval_after_spill<-which(between((list_analyte_time[[listID_i]]$date), spill_date, spill_date %m+% years(1)))
        if(length(row_ids_date_interval_after_spill)>0){
          to_obs_after_spill<-c(to_obs_after_spill,list_analyte_time[[listID_i]]$conc[row_ids_date_interval_after_spill])
        }
      }
    }
  }
  
  ####################################################
  if(c("analyte")%in%(df_node_latlong_anpoll[which(df_node_latlong_anpoll$nodeID==polluter_node_ID),"anpoll_indicator"])){
    #polluter_use_flag<-0
    listID<-which(listID_nodeID_matrix[,2]==polluter_node_ID)
    ## Getting all downstream observations with after one year spill date condition
    row_ids_date_interval_total<-which((between((list_analyte_time[[listID]]$date), date_start, date_end))&((between((list_analyte_time[[listID]]$date), spill_date, spill_date %m+% years(1)))))
    if(length(row_ids_date_interval_total)>0){
      to_obs_total<-c(to_obs_total,list_analyte_time[[listID]]$conc[row_ids_date_interval_total])
      #polluter_use_flag<-1
    }
    ## Getting all downstream observations before the spill date
    row_ids_date_interval_before_spill<-which(between((list_analyte_time[[listID]]$date), date_start, spill_date))
    if(length(row_ids_date_interval_before_spill)>0){
      to_obs_before_spill<-c(to_obs_before_spill,list_analyte_time[[listID]]$conc[row_ids_date_interval_before_spill])
    }
    ## Getting downstream observations till one year after the spill date
    row_ids_date_interval_after_spill<-which(between((list_analyte_time[[listID]]$date), spill_date, spill_date %m+% years(1)))
    if(length(row_ids_date_interval_after_spill)>0){
      to_obs_after_spill<-c(to_obs_after_spill,list_analyte_time[[listID]]$conc[row_ids_date_interval_after_spill])
    }
  }
  
  ########################################################################################################
  ## First temporal test before vs. after spill date for upstream 
  if(length(from_obs_before_spill)==0){
    from_mean_before_spill<-NA
    from_median_before_spill<-NA
  }else{
    from_mean_before_spill<-mean(from_obs_before_spill)
    from_median_before_spill<-median(from_obs_before_spill)
  }
  from_n_before_spill<-length(from_obs_before_spill)
  
  ####################################################
  if(length(from_obs_after_spill)==0){
    from_mean_after_spill<-NA
    from_median_after_spill<-NA
  }else{
    from_mean_after_spill<-mean(from_obs_after_spill)
    from_median_after_spill<-median(from_obs_after_spill)
  }
  from_n_after_spill<-length(from_obs_after_spill)
  ####################################################
  if((from_n_before_spill>0)&(from_n_after_spill>0)){
    mean_diff<-from_mean_after_spill-from_mean_before_spill
    if((from_n_before_spill>1)&(from_n_after_spill>1)){
      p_value_t.test_1sided_upstream<-t.test(x = from_obs_before_spill,y = from_obs_after_spill,alternative = "less")$p.value
      #p_value_wilcox.test_2sided<-wilcox.test(x = from_obs_total,y = to_obs_total)$p.value
      p_value_wilcox.test_1sided_upstream<-wilcox.test(x = from_obs_before_spill,y = from_obs_after_spill,alternative = "less")$p.value
      test_pass_upstream<-1
    }else{
      p_value_t.test_1sided_upstream<-NA
      #p_value_wilcox.test_2sided<-NA
      p_value_wilcox.test_1sided_upstream<-NA
      test_pass_upstream<-0
    }
  }else{
    p_value_t.test_1sided_upstream<-NA
    #p_value_wilcox.test_2sided<-NA
    p_value_wilcox.test_1sided_upstream<-NA
    test_pass_upstream<-0
  }
  summary_upstream<-c(p_value_t.test_1sided_upstream,p_value_wilcox.test_1sided_upstream,test_pass_upstream)
  
  ########################################################################################################
  ## Second temporal test before vs. after spill date for downstream 
  if(length(to_obs_before_spill)==0){
    to_mean_before_spill<-NA
    to_median_before_spill<-NA
  }else{
    to_mean_before_spill<-mean(to_obs_before_spill)
    to_median_before_spill<-median(to_obs_before_spill)
  }
  to_n_before_spill<-length(to_obs_before_spill)
  
  ####################################################
  if(length(to_obs_after_spill)==0){
    to_mean_after_spill<-NA
    to_median_after_spill<-NA
  }else{
    to_mean_after_spill<-mean(to_obs_after_spill)
    to_median_after_spill<-median(to_obs_after_spill)
  }
  to_n_after_spill<-length(to_obs_after_spill)
  ####################################################
  if((to_n_before_spill>0)&(to_n_after_spill>0)){
    mean_diff<-to_mean_after_spill-to_mean_before_spill
    if((to_n_before_spill>1)&(to_n_after_spill>1)){
      p_value_t.test_1sided_downstream<-t.test(x = to_obs_before_spill,y = to_obs_after_spill,alternative = "less")$p.value
      #p_value_wilcox.test_2sided<-wilcox.test(x = from_obs_total,y = to_obs_total)$p.value
      p_value_wilcox.test_1sided_downstream<-wilcox.test(x = to_obs_before_spill,y = to_obs_after_spill,alternative = "less")$p.value
      test_pass_downstream<-1
    }else{
      p_value_t.test_1sided_downstream<-NA
      #p_value_wilcox.test_2sided<-NA
      p_value_wilcox.test_1sided_downstream<-NA
      test_pass_downstream<-0
    }
  }else{
    p_value_t.test_1sided_downstream<-NA
    #p_value_wilcox.test_2sided<-NA
    p_value_wilcox.test_1sided_downstream<-NA
    test_pass_downstream<-0
  }
  summary_downstream<-c(p_value_t.test_1sided_downstream,p_value_wilcox.test_1sided_downstream,test_pass_downstream)
  
  ########################################################################################################
  ## Third spatio-temporal test for before and after spill date upstream observations vs after spill date downstream observations.
  if(length(from_obs_total)==0){
    from_mean_total<-NA
    from_median_total<-NA
  }else{
    from_mean_total<-mean(from_obs_total)
    from_median_total<-median(from_obs_total)
  }
  from_n_total<-length(from_obs_total)
  
  ####################################################
  if(length(to_obs_total)==0){
    to_mean_total<-NA
    to_median_total<-NA
  }else{
    to_mean_total<-mean(to_obs_total)
    to_median_total<-median(to_obs_total)
  }
  to_n_total<-length(to_obs_total)
  ####################################################
  if((from_n_total>0)&(to_n_total>0)){
    mean_diff<-to_mean_total-from_mean_total
    if((from_n_total>1)&(to_n_total>1)){
      p_value_t.test_1sided_total<-t.test(x = from_obs_total,y = to_obs_total,alternative = "less")$p.value
      #p_value_wilcox.test_2sided<-wilcox.test(x = from_obs_total,y = to_obs_total)$p.value
      p_value_wilcox.test_1sided_total<-wilcox.test(x = from_obs_total,y = to_obs_total,alternative = "less")$p.value
      test_pass_total<-1
    }else{
      p_value_t.test_1sided_total<-NA
      #p_value_wilcox.test_2sided<-NA
      p_value_wilcox.test_1sided_total<-NA
      test_pass_total<-0
    }
  }else{
    p_value_t.test_1sided_total<-NA
    #p_value_wilcox.test_2sided<-NA
    p_value_wilcox.test_1sided_total<-NA
    test_pass_total<-0
  }
  summary_total<-c(from_mean_total,to_mean_total,from_median_total,to_median_total,from_n_total,to_n_total,p_value_t.test_1sided_total,p_value_wilcox.test_1sided_total,test_pass_total)
  
  return(list(summary_upstream,summary_downstream,summary_total,from_obs_total,to_obs_total))
}

########################################################################################################
## Defining a function for FDR control using BH procedure, Inputs are 1) a vector of p values and 2) significance level alpha (usually set as 0.1). Outputs are 1) a vector of decisions (reject null hypothesis is 1) and 2) a flag indicating if we made atleast one discovery
fdr_decision_cal<-function(p_val,alpha){
  p_val_ordered<-sort(p_val)
  fdr_decision<-rep(0,length(p_val_ordered))
  flag_atleast_1_discovery<-0
  for (i in length(p_val_ordered):1){
    if(p_val_ordered[i]<=((i/length(p_val_ordered))*alpha)){
      i_max<-i
      flag_atleast_1_discovery<-1
      fdr_decision[order(p_val)[1:i_max]]<-1
      break()
    }
  }
  return(list(fdr_decision,flag_atleast_1_discovery))
}

## Defining a function to genereate final test result dataframe. Inputs are 1) polluter test matrix, 2) significance level alpha (usually set as 0.1) and 3) file_path. Output is the final dataframe with starred p values with or withour fdr control (depending on whether the number of tests is greater than or equal to 15)
fdr_analysis_wrapper<-function(polluter_test_matrix,alpha,county,water_body,file_path){
  
  ####################################################
  p_val_starred_generator<-function(p_val_mat,file_path){
    ## FDR Analysis and saving the results
    ## Getting the row IDs in polluter_test_matrix for which we have some test results
    test_pass_ids<-which(p_val_mat[,3]==1)
    ## Loading the dataframe df_polluter_processed
    load(file = paste0(file_path,"polluter_files/df_polluter_processed.RData"))
    ## Getting the p values for one sided tests
    p_values_t_test_one_sided_starred<-p_val_mat[test_pass_ids,1]
    p_values_wilcox_test_one_sided_starred<-p_val_mat[test_pass_ids,2]
    
    if(length(test_pass_ids)>=15){
      ## Getting the FDR decision for Wilcoxon two sided and one side tests
      fdr_decision_t_test_1_sided<-fdr_decision_cal(p_val = p_val_mat[test_pass_ids,1],alpha = alpha)[[1]]
      fdr_decision_wilcox_test_1_sided<-fdr_decision_cal(p_val = p_val_mat[test_pass_ids,2],alpha = alpha)[[1]]
      
      p_values_t_test_one_sided_starred<-round(as.numeric(p_values_t_test_one_sided_starred),2)
      p_values_wilcox_test_one_sided_starred<-round(as.numeric(p_values_wilcox_test_one_sided_starred),2)
      ## Appending the p-values with * for fdr decision one cases
      p_values_t_test_one_sided_starred[which(fdr_decision_t_test_1_sided==1)]<-paste0(p_values_t_test_one_sided_starred[which(fdr_decision_t_test_1_sided==1)],"*")
      p_values_wilcox_test_one_sided_starred[which(fdr_decision_wilcox_test_1_sided==1)]<-paste0(p_values_wilcox_test_one_sided_starred[which(fdr_decision_wilcox_test_1_sided==1)],"*")
    }else{
      p_values_t_test_one_sided_starred[which(p_values_t_test_one_sided_starred>0.05)]<-round(as.numeric(p_values_t_test_one_sided_starred[which(p_values_t_test_one_sided_starred>0.05)]),2)
      p_values_wilcox_test_one_sided_starred[which(p_values_wilcox_test_one_sided_starred>0.05)]<-round(as.numeric(p_values_wilcox_test_one_sided_starred[which(p_values_wilcox_test_one_sided_starred>0.05)]),2)
      p_values_t_test_one_sided_starred[which(p_values_t_test_one_sided_starred<=0.05)]<-paste0(round(p_values_t_test_one_sided_starred[which(p_values_t_test_one_sided_starred<=0.05)],2),"*")
      p_values_wilcox_test_one_sided_starred[which(p_values_wilcox_test_one_sided_starred<=0.05)]<-paste0(round(p_values_wilcox_test_one_sided_starred[which(p_values_wilcox_test_one_sided_starred<=0.05)],2),"*")
    }
    ## Intializing the full pvalue vectors 
    p_values_t_test<-rep(NA,nrow(p_val_mat))
    p_values_wilcox_test<-rep(NA,nrow(p_val_mat))
    ## Storing the p values starred
    p_values_t_test[test_pass_ids]<-p_values_t_test_one_sided_starred
    p_values_wilcox_test[test_pass_ids]<-p_values_wilcox_test_one_sided_starred
    return(list(p_values_t_test,p_values_wilcox_test))
  }
  
  ####################################################
  p_val_up<-p_val_starred_generator(p_val_mat = polluter_test_matrix[,c(7,8,9)],file_path = file_path)
  p_val_down<-p_val_starred_generator(p_val_mat = polluter_test_matrix[,c(10,11,12)],file_path = file_path)
  p_val_updown<-p_val_starred_generator(p_val_mat = polluter_test_matrix[,c(13,14,15)],file_path = file_path)
  ## Getting the final dataframe df polluter test
  #df_polluter_test<-data.frame(round(df_polluter_processed[test_pass_ids,c("lon","lat")],2),df_polluter_processed[test_pass_ids,"date"],round(polluter_test_matrix[test_pass_ids,c(1:4)],2),polluter_test_matrix[test_pass_ids,c(5,6)],p_values_t_test_one_sided_starred,p_values_wilcox_test_one_sided_starred)
  df_polluter_test_mean<-data.frame(water_body,county,as.character(df_polluter_processed[,"date"]),round(polluter_test_matrix[,c(1,2)],2),polluter_test_matrix[,c(5,6)],p_val_up[[1]],p_val_down[[1]],p_val_updown[[1]],stringsAsFactors = F)
  names(df_polluter_test_mean)<-c("Affected Water Body","County","Date","Upstream Mean","Downstream Mean","No. of Observations Upstream","No. of Observations Downstream","t test (up) p values","t test (down) p values","t test (updown) p values")
  
  df_polluter_test_median<-data.frame(water_body,county,as.character(df_polluter_processed[,"date"]),round(polluter_test_matrix[,c(3,4)],2),polluter_test_matrix[,c(5,6)],p_val_up[[2]],p_val_down[[2]],p_val_updown[[2]],stringsAsFactors = F)
  names(df_polluter_test_median)<-c("Affected Water Body","County","Date","Upstream Median","Downstream Median","No. of Observations Upstream","No. of Observations Downstream","Wilcox (up) p values","Wilcox (down) p values","Wilcox (updown) p values")
  return(list(df_polluter_test_mean,df_polluter_test_median))
}

########################################################################################################
########################################################################################################
