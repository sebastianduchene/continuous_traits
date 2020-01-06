# 1 trait, 100 tips
# Store Branch lengths in BM space too

library(ape,NELSI)
library(TreeSim)
library(phylotools)
library(phytools)
library(MASS)
library(phangorn)

setwd("~/Documents/GitHub/continuous_traits/test_sequences/1trait_100tips/")
#tr <- sim.bdsky.stt(n = 100, lambdasky = 2, deathsky = 1, sampprobsky = 0.1, timesky = 0)[[1]]
#write.tree(tr, file = "ReferenceTree_100tips_labmda2_death1_samprob0.1.nwk")
tr <- read.tree("ReferenceTree_100tips_labmda2_death1_samprob0.1.nwk")
tiplabs <- tr$tip.label

# Generating MVN Simulations

# 100 tips => 2N-2  - 198 (exclude root) 
# But, we assume the root node is fixed at zero so it is not an RV => 198 measurements for each trait

mu <- rep(0,198)
R <- 1
C <- vcvPhylo(tr, anc.nodes = T)
V <- kronecker(R, C)
Ntraits <- 1
Ntips <- length(tr$tip.label)
Nnodes <- 198
delta_trait <- list()
delta_time <- matrix(NA)
traits <- list()

sim_traits <- mvrnorm(n = 5, mu = mu, Sigma = V)


for (i in 1:Ntraits) {
  traits[[i]] <- sim_traits[,((i-1)*Nnodes+1):(Nnodes*i)]
}

# Adding in column for root node
Ntips <- length(tr$tip.label)
Nnodes <- 2*Ntips-2

for (i in 1:length(traits)){
  traits[[i]] <- cbind(traits[[i]][,1:Ntips], rep(0, length(sim_traits[,1])), traits[[i]][,(Ntips+1):Nnodes])
}

# Saving Node states
node_states <- cbind(c(1:5), traits[[1]])
colnames(node_states) <- c("replicate", paste0("t1_node", 1:199))
write.csv(node_states, row.names = F, file = "1trait_100tips_nodestates.csv")

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

# Diffusivity (D)
Path <- (sum(tr$edge.length))
D <-  rowSums(abs(delta_trait[[1]]))/Path

# Generating test sequences and saving with delta traits
DELTA_TRAITS <- cbind(paste0("N", 1:5), delta_trait[[1]])
colnames(DELTA_TRAITS) <- c("replicate", paste0("t1_node", 1:199))

write.csv(DELTA_TRAITS, row.names = F, file = "branch_delta_100tips_1trait.csv")

tr$tip.label <- paste0(tiplabs, '_', round(traits[[1]][N,1:Ntips],2), '_', 
                       round(allnode.times(tr, tipsonly = T, reverse = F), 2))

for (N in 1:length(sim_traits[,1])){
  
  phylogram <- tr
  
  phylogram$edge.length <- phylogram$edge.length * 1E-2
  aln <- as.DNAbin(simSeq(phylogram, l = 5000))
  
  write.dna(aln, file = paste0("BMtest1trait_sigma", R, 'D', round(D[N],2),'_100tips','N', N, '.fasta'), 
            format = 'fasta', nbcol = -1, colsep = '')
}
