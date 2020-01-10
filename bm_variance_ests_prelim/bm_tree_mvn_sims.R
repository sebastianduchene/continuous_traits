# .libPath("~/R/lib")
library(TreeSim)
library(phylotools)
library(phytools)
library(MASS)
library(motmot)
library(phangorn)
tr <- list()
sim_traits <- list()

# Data storage
bms <- data.frame() # Is where esdtimates will be stored
nodes <- data.frame() # storing ancestral nodes

# BM parameters
mu <- matrix(c(rep(0,198),rep(0,198)), ncol = 1)
R <- matrix(c(1, 0, 0, 1), 2)

for (i in 1:10){
	tr[[i]] <- sim.bdsky.stt(n = 100, lambdasky = 2, deathsky = 1, sampprobsky = 0.1, timesky = 0)[[1]]
	C <- vcvPhylo(tr[[i]], anc.nodes = T)
	V <- kronecker(R, C)
	sim_traits[[i]] <- mvrnorm(n = 10, mu = mu, Sigma = V)
	write.tree(tr[[i]], file = paste0("2traits_100tips_1var_0corr_rep", i, ".nwk"))
}


for (i in 1:length(sim_traits)){
	nodes <- rbind(nodes, sim_traits[[1]])
}

colnames(nodes) <- c(paste0("t1_n", c(1:100, 102:199)), paste0("t2_n", c(1:100, 102:199)))
write.csv(nodes, row.names = F,  file = "nodes_states_2traits_1var_0corr.csv")

for (i in 1:10){
	for (j in 1:length(sim_traits[[i]][,1])){

	tip_states <- cbind(sim_traits[[i]][j,1:100], sim_traits[[i]][j,200:299])
	row.names(tip_states) <- tr[[i]]$tip.label

	bms <- rbind(bms, c(i,j,as.vector(transformPhylo.ML(phy = tr[[i]], y = tip_states, model = "bm")$brownianVariance),
			as.vector(transformPhylo.ML(phy = tr[[i]], y = tip_states, model = "bm")$root.state),
			as.vector(transformPhylo.ML(phy = tr[[i]], y = tip_states, model = "lambda")$Lambda),
			as.vector(transformPhylo.ML(phy = tr[[i]], y = tip_states, model = "lambda")$brownianVariance)))
	}
}

colnames(bms) <- c("tree", "mvn_iteration", "var_t1", "covar_t12", "covar_t12", "var_t2", "t1_root_est", "t2_root_est", "pagels_lambda", "pagels_lambda_lowerCI", 'pagels_lambda_upperCI', "pagelvar_t1", "pagelcovar_t12", "pagelcovar_t12", "pagelvar_t2")
write.csv(bms, row.names = F, file  = "bm_2traits_100tips_10trees_1var_0corr_ests.csv")
