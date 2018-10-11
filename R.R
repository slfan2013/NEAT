pacman::p_load(data.table,nnet,jsonlite,RCurl)

############### initialization ############### 

in_node_init = c(1,2,3,4)
num_output = 2
out_node_init = max(in_node_init) + 1:num_output
weight_init = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1)

expressed_init = rep(TRUE, length(in_node_init)*length(out_node_init))
innovation_number_init = 1:(length(in_node_init)*length(out_node_init))
innovation_num = length(in_node_init)*length(out_node_init)
############### utility functions ############### 
generate_innovation = function(innovation_num){
  innovation_num <<- innovation_num+1
  return(innovation_num+1)
}

############### fitness score ############### 



############### hyperparameters ############### 
population_size = 150
c1 = 1
c2 = 1
c3 = 0.4
delta = 3
connection_weight_mutate_rate = 0.8
connection_weight_mutate_perturbed_rate = 0.9
connection_weight_mutate_new_value = 1-connection_weight_mutate_perturbed_rate
disabled_rage = 0.75
mutation_rate = 0.25
# crossover_rate = 0.75
interspecies_mating_rate = 0.001
add_new_node_rate = 0.03
link_mutation_rate = 0.05
add_new_link_rate = 0.3
transfer_function = function(x){
  return(1/(1+exp(-4.9*x)))
}


############### node initialization ###############
nodes = data.table(id = c(in_node_init,out_node_init),type = c(rep("INPUT",length(in_node_init)),rep("OUTPUT",length(out_node_init)))) # type can be INPUT, OUTPUT, HIDDEN.
############### connection_node initialization ############### 
connection_node = data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)), weight = weight_init, expressed = expressed_init, innovation_number = innovation_number_init)

############### cynodes  ###############
cylabel = "Compound_Name"

change_cynode_position_by_suid = function(SUID,x,y){
  RCurl::getURL('http://localhost:1234/v1/networks/52/views/156/nodes',customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(
    list(list(
      SUID=SUID,
      view = list(
        list(
          visualProperty= "NODE_X_LOCATION",
          value=x
        ),
        list(
          visualProperty= "NODE_Y_LOCATION",
          value=y
        )
      )
    ))
    ,auto_unbox = T, force = T))
}
change_cynode_position_by_suid(75,-11.738943,-35.5222373)




## get node suid, x, y, h, w (label, label font)
cynodes = data.table(suid = fromJSON("http://localhost:1234/v1/networks/52/nodes"), fromJSON("http://localhost:1234/v1/networks/52/views/156")$elements$nodes$position,label = fromJSON("http://localhost:1234/v1/networks/52/views/156")$elements$nodes$data[[cylabel]])
cynodes = merge(cynodes, data.table(suid = fromJSON("http://localhost:1234/v1/networks/52/views/156/nodes")$SUID,labelfont = as.numeric(sapply(fromJSON("http://localhost:1234/v1/networks/52/views/156/nodes")$view,function(x){x$value[x$visualProperty == "NODE_LABEL_FONT_SIZE"]}))),by = 'suid',sort = FALSE)
cynodes = merge(cynodes, data.table(suid = fromJSON("http://localhost:1234/v1/networks/52/views/156/nodes")$SUID,h = as.numeric(sapply(fromJSON("http://localhost:1234/v1/networks/52/views/156/nodes")$view,function(x){x$value[x$visualProperty == "NODE_HEIGHT"]}))),by = 'suid',sort = FALSE)

cynodes$w = strwidth(cynodes$label, font = cynodes$labelfont, units = 'in')/strheight(cynodes$label, font = cynodes$labelfont, units = 'in') * cynodes$labelfont
cynodes$w[is.nan(cynodes$w)] = cynodes$h[is.nan(cynodes$w)]
# cynodes = cynodes[,c('suid','x','y','h','w')]
cynodes$rx = cynodes$x + 1/2*cynodes$w
cynodes$lx = cynodes$x - 1/2*cynodes$w
cynodes$by = cynodes$y + 1/2*cynodes$h
cynodes$ty = cynodes$y - 1/2*cynodes$h
##
## get node iteration sequence
# 1 calculate the center of the cluster.
center_xy = colMeans(cynodes[,c('x','y')])
# 2 the distance between each node to the center
dists_to_center = sqrt(colSums((t(cynodes[,c('x','y')]) - center_xy)^2))
# 3 get the order
iteration_sequence = order(dists_to_center, decreasing = FALSE)
##
## find the overlapping for each direction (The input of the neural network)
# apply(data.matrix(cynodes),1,function(x){
#   x$rx > cynodes$lx & (x$ty > cynodes$by &  x$by < cynodes$ty)
# })
# # to determine overlapping nodes.
# cynodes[,3]$rx > cynodes$lx & cynodes[,3]$ty < cynodes$by &  cynodes[,3]$by > cynodes$ty
# sapply(l, function(x) sapply(l, function(y) foo(x,y)))
# 

# check overlapping
check_overlapping = function(cynode1, cynode2){
  # If one rectangle is on left side of other
  if (cynode1$lx > cynode2$rx || cynode2$lx > cynode1$rx)
    {return(FALSE)}
  
  # If one rectangle is above other
  if (cynode1$ty > cynode2$by || cynode2$ty > cynode1$by)
    {return(FALSE)}
  
  return(TRUE)
}

# for(i in 1:nrow(cynodes)){
#   print(cynodes$label[i])
#   print(check_overlapping(cynodes[12,],cynodes[i,]))
# }


indices = cbind(combn(1:ncol(cynodes),2,simplify = T),combn(ncol(cynodes):1,2,simplify = T))
indices = indices[,order(indices[1,],decreasing = FALSE)]
overlapping_index = apply(indices,2,function(ind){
  return(check_overlapping(cynode1 = cynodes[ind[1],], cynodes[ind[2],]))
})
overlapping_index_for_each_cynode = by(t(rbind(indices[2,],overlapping_index)),indices[1,],function(x){
  # for(i in 1:nrow(cynodes)){
    # x = t(rbind(indices[2,],overlapping_index))[indices[1,]==1,]
    # print(x[1,][x[2,]==1])
  # }
  return(x[,1][x[,2]==1])
})
## generate input for neural network.
# the inputs are the most left x distance overlapped by others, most right x distance, most top y distance, and most bottom y distance. (by is larger than ty)
generate_nn_inputs_for_each_cynode = function(cynodes){
  return(mapply(function(x,ind){
    if(length(x)>0){
      return(abs(c(
        max(-cynodes[ind,]$w,min(cynodes[x,]$lx - cynodes[ind,]$rx)),
        min(cynodes[ind,]$w,max(cynodes[x,]$rx - cynodes[ind,]$lx)),
        max(-cynodes[ind,]$h,min(cynodes[x,]$ty - cynodes[ind,]$by)),
        min(cynodes[ind,]$h,max(cynodes[x,]$by - cynodes[ind,]$ty))
      )))
    }else{
      return(c(0,0,0,0))
    }
    
  },overlapping_index_for_each_cynode,1:length(overlapping_index_for_each_cynode)))
}

# calculate score (fitness). #overlapping distance + moving distance
moving_dist = 0
generate_overlapping_dist = function(cynodes,overlapping_index_for_each_cynode){
  return(sum(mapply(function(x,ind){
    if(length(x)>0){
      return(
        sum(abs(c(cynodes[x,]$lx - cynodes[ind,]$rx,cynodes[x,]$rx - cynodes[ind,]$lx,cynodes[x,]$ty - cynodes[ind,]$by,cynodes[x,]$by - cynodes[ind,]$ty)))
        )
    }else{
      return(0)
    }
    
  },overlapping_index_for_each_cynode,1:length(overlapping_index_for_each_cynode))))
}
generate_overlapping_dist(cynodes,overlapping_index_for_each_cynode)


# evaluate current status
evaluate_current_status = function(cynodes){
  # check overlapping
  indices = cbind(combn(1:ncol(cynodes),2,simplify = T),
                  combn(ncol(cynodes):1,2,simplify = T))
  indices = indices[,order(indices[1,],decreasing = FALSE)]
  overlapping_index = apply(indices,2,function(ind){
    return(check_overlapping(cynode1 = cynodes[ind[1],], cynodes[ind[2],]))
  })
  overlapping_index_for_each_cynode = by(t(rbind(indices[2,],overlapping_index)),indices[1,],function(x){return(x[,1][x[,2]==1])})
  return(generate_overlapping_dist(cynodes,overlapping_index_for_each_cynode) + moving_dist)
  
  
}



change_cynode_lables_to_cynodes_index = function(){
  
  SUIDs = fromJSON("http://localhost:1234/v1/networks/52/views/156/nodes")$SUID
  for(SUID in SUIDs ){
    RCurl::getURL('http://localhost:1234/v1/networks/52/views/156/nodes',customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(
      list(list(
        SUID=SUID,
        view = list(
          list(
            visualProperty= "NODE_LABEL",
            value=which(cynodes$suid == SUID)
          )
        )
      ))
      ,auto_unbox = T, force = T))
  }
  
  
}



############### initial population ############### 
population_tracker = list() # remember the generations of the population
genomes = list() # each genome in a generation
for(i in 1:population_size){
  genomes[[i]] =  data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)),weight = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1),in_layer = 1, out_layer = 2, in_seq = rep(1:length(in_node_init), length(out_node_init)),out_seq = rep(1:length(out_node_init), each = length(in_node_init)), expressed = expressed_init, innovation_number = innovation_number_init)
}



genome = data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)),weight = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1),in_layer = 1, out_layer = 2, in_seq = rep(1:length(in_node_init), length(out_node_init)),out_seq = rep(1:length(out_node_init), each = length(in_node_init)), expressed = expressed_init, innovation_number = innovation_number_init)
genome$out_layer = genome$out_layer+1
genome = rbind(genome, data.table(in_node = 1,out_node = max(genome$in_node,genome$out_node),weight = runif(1,min=-1,max=1),in_layer = genome$in_layer[genome$in_node==1][1],out_layer=genome$in_layer[genome$in_node==1][1]+1,in_seq=1,out_seq=1,expressed=TRUE,innovation_number = 9))

eval_nodes = function(genome,iteration_sequence){
  # get the input from current status (cynode)
  generate_nn_inputs_for_each_cynode(cynodes)
  #!!!
  iteration_sequence
  sapply(genomes, function(genome){
    
    ## apply nn
    # generate matrix.
    x = 1:4 # (input values)
    for(i in 2:max(genome$out_layer)){
      W = matrix(0,nrow = sum(genome$out_layer==i), ncol = max(genome$in_node[genome$in_layer==i-1]))
      W[genome[genome$out_layer==i,]$in_seq,genome[genome$out_layer==i,]$out_seq] = genome[genome$out_layer==i,]$weight
      
      
      wx = W%*%x
      
    }
    
    
    
  })
  
  
}































############### mutation ############### 
