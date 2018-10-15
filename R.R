##
# to work from office, change 4321 to 4321
#
##

pacman::p_load(data.table,nnet,jsonlite,RCurl, igraph)

############### initialization ############### 

in_node_init = c(1,2,3,4)
num_output = 2
out_node_init = max(in_node_init) + 1:num_output
weight_init = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1)
expressed_init = rep(TRUE, length(in_node_init)*length(out_node_init))
innovation_number_init = 1:(length(in_node_init)*length(out_node_init))
innovation_num = length(in_node_init)*length(out_node_init)

number_of_iteration = 5

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
crossover_rate = 0.75
interspecies_mating_rate = 0.001
add_new_node_rate = 0.03
link_mutation_rate = 0.05
activate_function = function(x){
  return(
    x
  )
}


############### node initialization ###############
nodes = data.table(id = c(in_node_init,out_node_init),type = c(rep("INPUT",length(in_node_init)),rep("OUTPUT",length(out_node_init)))) # type can be INPUT, OUTPUT, HIDDEN.
############### connection_node initialization ############### 
connection_node = data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)), weight = weight_init, expressed = expressed_init, innovation_number = innovation_number_init)

############### cynodes  ###############
cylabel = "Compound_Name"

change_cynode_position_by_suid = function(SUID,x,y){
  RCurl::getURL('http://localhost:4321/v1/networks/52/views/156/nodes',customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(
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
change_cynode_position_by_suid(SUID = 221,x = -6.416754 ,y = -123.30798)




## get node suid, x, y, h, w (label, label font)
cynodes = data.table(suid = fromJSON("http://localhost:4321/v1/networks/52/nodes"), fromJSON("http://localhost:4321/v1/networks/52/views/156")$elements$nodes$position,label = fromJSON("http://localhost:4321/v1/networks/52/views/156")$elements$nodes$data[[cylabel]])
cynodes = merge(cynodes, data.table(suid = fromJSON("http://localhost:4321/v1/networks/52/views/156/nodes")$SUID,labelfont = as.numeric(sapply(fromJSON("http://localhost:4321/v1/networks/52/views/156/nodes")$view,function(x){x$value[x$visualProperty == "NODE_LABEL_FONT_SIZE"]}))),by = 'suid',sort = FALSE)
cynodes = merge(cynodes, data.table(suid = fromJSON("http://localhost:4321/v1/networks/52/views/156/nodes")$SUID,h = as.numeric(sapply(fromJSON("http://localhost:4321/v1/networks/52/views/156/nodes")$view,function(x){x$value[x$visualProperty == "NODE_HEIGHT"]}))),by = 'suid',sort = FALSE)

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
#   print(check_overlapping(cynodes[8,],cynodes[i,]))
# }


indices = cbind(combn(1:nrow(cynodes),2,simplify = T),combn(nrow(cynodes):1,2,simplify = T))
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
generate_nn_inputs_for_each_cynode(cynodes)
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
  indices = cbind(combn(1:nrow(cynodes),2,simplify = T), combn(nrow(cynodes):1,2,simplify = T))
  indices = indices[,order(indices[1,],decreasing = FALSE)]
  overlapping_index = apply(indices,2,function(ind){
    return(check_overlapping(cynode1 = cynodes[ind[1],], cynodes[ind[2],]))
  })
  overlapping_index_for_each_cynode = by(t(rbind(indices[2,],overlapping_index)),indices[1,],function(x){return(x[,1][x[,2]==1])})
  return(generate_overlapping_dist(cynodes,overlapping_index_for_each_cynode) + moving_dist)
  
  
}
evaluate_current_status(cynodes)


change_cynode_lables_to_cynodes_index = function(){
  
  SUIDs = fromJSON("http://localhost:4321/v1/networks/52/views/156/nodes")$SUID
  for(SUID in SUIDs ){
    RCurl::getURL('http://localhost:4321/v1/networks/52/views/156/nodes',customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(
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
# change_cynode_lables_to_cynodes_index()


############### initial population ############### 
population_tracker = list() # remember the generations of the population
genomes = list() # each genome in a generation
for(i in 1:population_size){
  # genomes[[i]] =  data.table(in_node = c(rep(in_node_init, length(out_node_init)), out_node_init), out_node = c(rep(out_node_init, each = length(in_node_init)),out_node_init),weight = c(runif(length(in_node_init)*length(out_node_init), min = -1, max = 1), rep(0,length(out_node_init))),in_layer = 1, out_layer = 2, in_seq = rep(1:length(in_node_init), length(out_node_init)),out_seq = rep(1:length(out_node_init), each = length(in_node_init)), expressed = expressed_init, innovation_number = innovation_number_init)
  genomes[[i]] =  data.table(in_node = c(rep(in_node_init, length(out_node_init)), out_node_init), out_node = c(rep(out_node_init, each = length(in_node_init)),out_node_init),weight = c(runif(length(in_node_init)*length(out_node_init), min = -1, max = 1), rep(0,length(out_node_init))),in_layer = 1, out_layer = 2, expressed = TRUE, innovation_number = c(innovation_number_init,rep(0,length(out_node_init))))
  genomes[[i]][in_node == out_node,in_layer := out_layer]
}


# genome = data.table(in_node = c(1,2,5,3,6,1,1,7,4), out_node = c(4,5,4,6,4,5,7,4,4), in_layer = c(1,1,2,1,2,1,1,2,3), out_layer = c(3,2,3,2,3,2,2,3,3), expressed = TRUE)
# 
# set.seed(1)
# genome$weight = c(runif(nrow(genome)-1),0)
# genome$innovation_num = c(1:(nrow(genome)-1),0)
# genome$value = c(1,2,0,3,0,1,1,0,0)
# genome$expressed[1] = FALSE

eval_genome = function(genome){
  result = genome
  result$in_node[1] = 1 # https://stackoverflow.com/questions/10225098/understanding-exactly-when-a-data-table-is-a-reference-to-vs-a-copy-of-another
  seq_out = 2:max(result$out_layer)
  
  for(i in seq_out){
    # print(i)
    genome_sub = result[out_layer==i,]
    
    out_node_values = result[out_layer==i,out_node]
    for(out_node_value in unique(out_node_values)){
      # print(out_node_value)
      genome_sub2 = genome_sub[out_node==out_node_value & expressed]
      
      gene_output_value = activate_function(sum(genome_sub2[,value] * genome_sub2[,weight]))
      
      result[in_node==out_node_value,value:=gene_output_value]
    }
  }
  return(result)
  
  
}
# result_genome = eval_genome(genome)
cynodes_origin = cynodes


cynodes_inputs = generate_nn_inputs_for_each_cynode(cynodes)



for(genome_index in 1:length(genomes)){
  print(genome_index)
  genome = genomes[[genome_index]]
  for(num_iter in 1:number_of_iteration){
    print(num_iter)
    for(input_index in 1:ncol(cynodes_inputs)){
      input = cynodes_inputs[,input_index]
      genome[!in_node == out_node,value:=input]
      genome[in_node == out_node,value:=0]
      evaluated_genome = eval_genome(genome)
      # change position
      suid_index = input_index
      suid = cynodes$suid[[suid_index]]
      
      dx = evaluated_genome[in_node==out_node,value][1]
      newx = cynodes$x[cynodes$suid==suid] + dx
      dy = evaluated_genome[in_node==out_node,value][2]
      newy = cynodes$y[cynodes$suid==suid] + dy
      if(suid == "73"){
        print(dx)
        print(dy)
      }
      change_cynode_position_by_suid(SUID = suid,x = newx ,y = newy)
      # update cynode
      cynodes$x[cynodes$suid==suid] = newx
      cynodes$y[cynodes$suid==suid] = newy
      moving_dist = moving_dist + abs(dx)
      moving_dist = moving_dist + abs(dy)
    }
  }
  # evaluate the current status (after iterations)
  evaluate_current_status(cynodes)
  print(moving_dist)
}


reset = function(cynodes){
  for(i in 1:nrow(cynodes)){
    change_cynode_position_by_suid(SUID = cynodes$suid[i],x = cynodes$x[i] ,y = cynodes$y[i])
  }
}
reset(cynodes_origin)









genome = data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)),weight = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1),in_layer = 1, out_layer = 2, in_seq = rep(1:length(in_node_init), length(out_node_init)),out_seq = rep(1:length(out_node_init), each = length(in_node_init)), expressed = expressed_init, innovation_number = innovation_number_init)



genome$out_layer = genome$out_layer+1



# adding a new connection.
genome = rbind(genome, data.table(in_node = 1,out_node = max(genome$in_node,genome$out_node)+1,weight = runif(1,min=-1,max=1),in_layer = genome$in_layer[genome$in_node==1][1],out_layer=genome$in_layer[genome$in_node==1][1]+1,in_seq=1,out_seq=1,expressed=TRUE,innovation_number = 9))
genome = rbind(genome, data.table(in_node = 7,out_node = 6,weight = runif(1,min=-1,max=1),in_layer = 2,out_layer=3,in_seq=1,out_seq=2,expressed=TRUE,innovation_number = 10))
# adding a new connection.
genome = rbind(genome, data.table(in_node = 4,out_node = 8,weight = runif(1,min=-1,max=1),in_layer = 1,out_layer=2,in_seq=4,out_seq=2,expressed=TRUE,innovation_number = 11))
genome = rbind(genome, data.table(in_node = 8,out_node = 5,weight = runif(1,min=-1,max=1),in_layer = 2,out_layer=3,in_seq=2,out_seq=1,expressed=TRUE,innovation_number = 12))

visual_net = function(genome){
  nodes = data.frame(id = 1:max(genome$in_node,genome$out_node),label = 1:max(genome$in_node,genome$out_node))
  links = data.frame(from = genome$in_node, to = genome$out_node, weight = genome$weight)
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
  plot(net)
}
visual_net(genome)


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
      
      
      W = matrix(0,nrow = length(unique(genome$out_node[genome$out_layer==i])), ncol = length(unique(genome$in_node[genome$in_layer<i])))
      W[genome[genome$out_layer==i,]$out_seq,genome[genome$out_layer==i,]$in_seq + cumsum(genome[genome$out_layer==i,]$in_layer)] = genome[genome$out_layer==i,]$weight
      
      W = matrix(0,nrow = sum(genome$out_layer==i), ncol = max(genome$in_node[genome$in_layer<i]))
      W[genome[genome$out_layer==i,]$in_seq,genome[genome$out_layer==i,]$out_seq] = genome[genome$out_layer==i,]$weight
      
      
      wx = W%*%x
      
    }
    
    
    
  })
  
  
}































############### mutation ############### 
