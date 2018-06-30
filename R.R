pacman::p_load(data.table,nnet,jsonlite)

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
crossover_rate = 0.75
interspecies_mating_rate = 0.001
add_new_node_rate = 0.03
link_mutation_rate = 0.05
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
  RCurl::getURL('http://localhost:4321/v1/networks/52/views/140/nodes',customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(
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
change_cynode_position_by_suid(221,-6.416754 ,-123.30798)

## get node suid, x, y, h, w (label, label font)
cynodes = data.table(suid = fromJSON("http://localhost:4321/v1/networks/52/nodes"), fromJSON("http://localhost:4321/v1/networks/52/views/140")$elements$nodes$position,label = fromJSON("http://localhost:4321/v1/networks/52/views/140")$elements$nodes$data[[cylabel]])
cynodes = merge(cynodes, data.table(suid = fromJSON("http://localhost:4321/v1/networks/52/views/140/nodes")$SUID,labelfont = as.numeric(sapply(fromJSON("http://localhost:4321/v1/networks/52/views/140/nodes")$view,function(x){x$value[x$visualProperty == "NODE_LABEL_FONT_SIZE"]}))),by = 'suid',sort = FALSE)
cynodes = merge(cynodes, data.table(suid = fromJSON("http://localhost:4321/v1/networks/52/views/140/nodes")$SUID,h = as.numeric(sapply(fromJSON("http://localhost:4321/v1/networks/52/views/140/nodes")$view,function(x){x$value[x$visualProperty == "NODE_HEIGHT"]}))),by = 'suid',sort = FALSE)

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
apply(cynodes,1,function(x){
  x$rx > cynodes$lx & (x$ty > cynodes$by | x$by < cynodes$ty)
})


############### initial population ############### 
population_tracker = list() # remember the generations of the population
genomes = list() # each genome in a generation
for(i in 1:population_size){
  genomes[[i]] =  data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)), weight = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1), expressed = expressed_init, innovation_number = innovation_number_init)
}



genome = data.table(in_node = rep(in_node_init, length(out_node_init)), out_node = rep(out_node_init, each = length(in_node_init)), weight = runif(length(in_node_init)*length(out_node_init), min = -1, max = 1), expressed = expressed_init, innovation_number = innovation_number_init)

eval_nodes = function(genome){
  
}































############### mutation ############### 
