require(compiler)
enableJIT(3)
source("Br_projected_clust.R")
args = commandArgs()
chunk_ID = as.numeric(args[which(args=="--args")+1])
system.time(wrapper_projected(chunk_ID))

