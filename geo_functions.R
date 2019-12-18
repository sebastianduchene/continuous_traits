get.geo.data.matrix <- function(tr, sim_trait, tip_labels){
    geo.data.matrix <- matrix(NA, nrow(tr$edge), 8)
    colnames(geo.data.matrix) <- c('node_from', 'node_to', 'x1', 'x2', 'y1',
                                   'y2', 'edge_length', 'distance')
    for(i in 1:nrow(geo.data.matrix)){
        if(tr$edge[i, 1] == length(tr$tip.label)+1){ # If 'from' is the root set to 0
            from1 <- 0
            from2 <- 0
        }else{
            from1 <- sim_trait[grep(paste0(tr$edge[i, 1], '_1'), names(sim_trait))]
            from2 <- sim_trait[grep(paste0(tr$edge[i, 1], '_2'), names(sim_trait))]
        }
        if(tr$edge[i, 2] <= length(tr$tip.label)){
            to1 <- sim_trait[grep(paste0(tip_labels[tr$edge[i, 2]], '_1'), names(sim_trait))]
            to2 <- sim_trait[grep(paste0(tip_labels[tr$edge[i, 2]], '_2'), names(sim_trait))]
        }else{
            to1 <- sim_trait[grep(paste0(tr$edge[i, 2], '_1'), names(sim_trait))]
            to2 <- sim_trait[grep(paste0(tr$edge[i, 2], '_2'), names(sim_trait))]
        }
        distance <- sqrt(((from1 - to1)^2) + ((from2 - to2)^2))
        geo.data.matrix[i, ] <- c(tr$edge[i, ], from1, to1, from2, to2, tr$edge.length[i], distance)
    }
    return(geo.data.matrix)
}
