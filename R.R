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
distance_criterion = 3
connection_weight_mutate_rate = 0.8
connection_weight_mutate_perturbed_rate = 0.9
connection_weight_mutate_new_value = 1-connection_weight_mutate_perturbed_rate
disabled_rage = 0.75
mutation_rate = 0.25
# crossover_rate = 0.75
interspecies_mating_rate = 0.001
add_new_node_rate = 0.03
link_mutation_rate = 0.05

activate_function = function(x){
  return(
    x
  )
}
add_new_link_rate = 0.3

minimum_species_mumber = 5
mutation_ratio = 0.25
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
cynodes= cynodes[iteration_sequence]
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



## generate input for neural network.
# the inputs are the most left x distance overlapped by others, most right x distance, most top y distance, and most bottom y distance. (by is larger than ty)
generate_nn_inputs_for_each_cynode = function(cynodes){
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


reset = function(cynodes){
  for(i in 1:nrow(cynodes)){
    change_cynode_position_by_suid(SUID = cynodes$suid[i],x = cynodes$x[i] ,y = cynodes$y[i])
  }
}

reset(cynodes_origin)
genome_scores = c()
species_index = rep(c("S1","S2"),length(genomes)/2)
for(genome_index in 1:length(genomes)){
  print(genome_index)
  reset(cynodes_origin)
  cynodes = cynodes_origin
  moving_dist = 0 # reset to zero for next genome
  genome = genomes[[genome_index]]
  for(num_iter in 1:number_of_iteration){
    print(num_iter)
    cynodes_inputs = generate_nn_inputs_for_each_cynode(cynodes)
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
      cynodes$rx[cynodes$suid==suid]= cynodes$rx[cynodes$suid==suid] + dx
      cynodes$lx[cynodes$suid==suid] = cynodes$lx[cynodes$suid==suid] + dx
      dy = evaluated_genome[in_node==out_node,value][2]
      newy = cynodes$y[cynodes$suid==suid] + dy
      cynodes$by[cynodes$suid==suid] = cynodes$by[cynodes$suid==suid] + dy
      cynodes$ty[cynodes$suid==suid] = cynodes$ty[cynodes$suid==suid] + dy
      change_cynode_position_by_suid(SUID = suid,x = newx ,y = newy)
      # update cynode
      cynodes$x[cynodes$suid==suid] = newx
      cynodes$y[cynodes$suid==suid] = newy
      moving_dist = moving_dist + abs(dx)
      moving_dist = moving_dist + abs(dy)
    }
  }
  # evaluate the current status (after iterations)
  genome_scores[genome_index] = evaluate_current_status(cynodes)
}
########### now it is time to generate the next generation
## speciation
# randomly select a species_representative from each species
species_representative_index = by(1:length(species_index),species_index,function(x){
  return(sample(1:length(unique(x)),1))
})
species_representative_index_temp = paste0(names(species_representative_index),"_",species_representative_index)
DT = data.table(species_index)
DT[, id := seq_len(.N), by = species_index]
species_index_temp = paste0(species_index,"_",DT$id)
species_representative_index = which(species_index_temp%in%species_representative_index_temp)
names(species_representative_index) = names(table(species_index))

distance_measure = function(genome, representative = genomes[[5]]){
  N = max(length(unique(genome[,in_node],genome[,out_node])),length(unique(representative[,in_node],representative[,out_node])))
  gene_match_genome = genome[innovation_number%in%representative$innovation_number]
  gene_match_representative = representative[innovation_number%in%genome$innovation_number]
  E = max(c(genome$innovation_number,representative$innovation_number)) - (min(max(genome$innovation_number, representative$innovation_number)))
  D = max(c(genome$innovation_number,representative$innovation_number))-(max(intersect(gene_match_genome$innovation_number,gene_match_representative$innovation_number)))
  W_bar = sum(abs(gene_match_genome$weight - gene_match_representative$weight))
  
  distance = c1*E/N + c2*D/N + c3*W_bar
  return(distance)
}
# use distance measure to update the species
for(genome_index in 1:length(genomes)){
  genome = genomes[[genome_index]]
  species_representative_distances = c()
  for(representative_index in 1:length(species_representative_index)){
    species_representative_distances[representative_index] = distance_measure(genome, genomes[[species_representative_index[representative_index]]])
  }
  if(any(species_representative_distances<distance_criterion)){ # if the distance is small, it means this genome belongs to a existing species
    species_index[genome_index] = names(species_representative_index)[which.min(species_representative_distances)]
  }else{ # otherwise, the genome is a new species.
    species_index[genome_index] = paste0("S",max(as.numeric(gsub("S", "", unique(species_index))))+1)
  }
}
# The champion of each species with more than five networks was copied into the next generation unchanged.
champion_index = by(genome_scores, species_index,which.min)
champion_index_temp = paste0(names(champion_index),"_",champion_index)
DT = data.table(species_index)
DT[, id := seq_len(.N), by = species_index]
species_index_temp = paste0(species_index,"_",DT$id)
champion_index = which(species_index_temp%in%champion_index_temp)
names(champion_index) = names(table(species_index))
remaining_champion_index = champion_index[table(species_index)>minimum_species_mumber]

next_genomes = list()
next_genomes = genomes[remaining_champion_index]

perturb = function(x,sd=0.2){
  return(x + rnorm(1,mean=0,sd))
}
## generate new genome
# mutation
genome_mutation_index = sample(1:length(genomes), mutation_ratio * length(genomes))

for(genome_index in 1:length(genomes)){
  genome = genomes[[genome_index]]
  
  # mutate 1
  # There was an 80% chance of a genome having its connection weights mutated
  if(runif(1)<connection_weight_mutate_rate){
    # in which case each weight had a 90% chance of being uniformly perturbed
    connection_weight_mutate_perturbed_rate_index = runif(genome$weight)<connection_weight_mutate_perturbed_rate
    for(mutate_index in 1:length(connection_weight_mutate_perturbed_rate_index)){
      if(!genome[mutate_index,in_node == out_node]){ # just to make sure that the currently gene is not a output node link.
        if(connection_weight_mutate_perturbed_rate_index[mutate_index]){
          genome[mutate_index,weight:=perturb(weight)]
        }else{ # a 10% chance of being assigned a new random value.
          genome[mutate_index,weight:=runif(1, min = -1, max = 1)]
        }
      }
    }
  }else{ # mutate 2
    
  }
  
  
  
  
}


























