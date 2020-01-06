# 2 traits, 100 tips
# Store Branch lengths in BM space too

library(ape,NELSI)
library(TreeSim)
library(phylotools)
library(phytools)
library(MASS)
library(phangorn)

setwd("~/Documents/GitHub/continuous_traits/test_sequences/2traits_100tips/")
tr <- sim.bdsky.stt(n = 100, lambdasky = 2, deathsky = 1, sampprobsky = 0.1,
timesky = 0)[[1]]
write.tree(tr, file = "ReferenceTree_100tips_labmda2_death1_samprob0.1.nwk")
tr <- read.tree("ReferenceTree_100tips_labmda2_death1_samprob0.1.nwk")
tiplabs <- tr$tip.label

# Generating MVN Simulations

# 100 tips => 2N-2  - 198 (exclude root) 
# But, we assume the root node is fixed at zero so it is not an RV => 198 measurements for each trait

mu <- matrix(c(rep(0,198),rep(0,198)), ncol = 1)
R <- matrix(c(1, 0, 0, 1), 2)
C <- vcvPhylo(tr, anc.nodes = T)
V <- kronecker(R, C)
Ntraits <- 2
Ntips <- length(tr$tip.label)
Nnodes <- 198
delta_trait <- list()
delta_time <- matrix(NA)
traits <- list()

sim_traits <- mvrnorm(n = 5, mu = mu, Sigma = V)

# Break into Ntraits = 2 trait matrices of p obs x n nodes. Store in lists. Each element of the list is the matrix for a given trait

for (i in 1:Ntraits) {
  traits[[i]] <- sim_traits[,((i-1)*Nnodes+1):(Nnodes*i)]
}

# Adding in column for root node
Ntips <- length(tr$tip.label)
Nnodes <- 2*Ntips-2

for (i in 1:length(traits)){
  traits[[i]] <- cbind(traits[[i]][,1:Ntips], rep(0, length(sim_traits[,1])), traits[[i]][,(Ntips+1):Nnodes])
}

# Getting delta trait per branch using edge objct
v1 <- tr$edge[,1]
v2 <- tr$edge[,2]

for (t in 1:length(traits)){
  delta_trait[[t]] <- matrix(NA, nrow = nrow(traits[[1]]), ncol = ncol(traits[[1]]))
}

for (t in 1:length(traits)){
  
  delta_trait[[t]][,(Ntips+1)] <- rep(0, length(traits[[1]][,1]))
  
}

for (t in 1:length(traits)){
  for (i in 1:length(v1)){
    
    delta_trait[[t]][,v2[i]] <- as.numeric(traits[[t]][, v2[i]] - traits[[t]][, v1[i]])
    
  }
}

# NB: Number of columns is now 199 instead of 198 (2N-1 instead of 2N-2)

# Euclidean Distance
distance <- sqrt((delta_trait[[1]])^2 + (delta_trait[[2]])^2)

# Diffusivity (D)
Path <- (sum(tr$edge.length))
D <-  rowSums(distance)/Path

# Generating test sequences and saving with delta traits

DELTA_TRAITS <- matrix()

for (t in 1:(length(traits)-1)){
  DELTA_TRAITS <- cbind(delta_trait[[t]], delta_trait[[t+1]])
}

colnames(DELTA_TRAITS) <- c(paste0("t1_node", 1:199), paste0("t2_node", 1:199))
write.csv(cbind(paste0("N", 1:5), DELTA_TRAITS), row.names = F, file = "branch_delta_100tips_2traits.csv")

tr$tip.label <- paste0(tiplabs, '_', round(traits[[1]][N,1:Ntips],2), '_', round(traits[[2]][N,1:Ntips],2), '_', 
                       round(allnode.times(tr, tipsonly = T, reverse = F), 2))

for (N in 1:length(sim_traits[,1])){

  phylogram <- tr

  
  phylogram$edge.length <- phylogram$edge.length * 1E-2
  aln <- as.DNAbin(simSeq(phylogram, l = 5000))
  
  corr <- R[1,2]/(sqrt(R[1,1]*R[2,2]))
  
  write.dna(aln, file = paste0("BMtest2traits_sigma", R[1,1], '_', R[2,2],"corr", corr, '_', corr, 'D', round(D[N],2),'_100tips','N', N, '.fasta'), 
            format = 'fasta', nbcol = -1, colsep = '')
}
