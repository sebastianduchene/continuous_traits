library(NELSI)
library(phytools)
library(phangorn)
library(geiger)

tr <- read.tree('ReferenceTree_20tips_labmda2_death1_samprob0.1.nwk')

aln <- read.dna('BMtest_sigma1_1corr0_0D1.24N9.fasta', format = 'fasta')

x <- as.numeric(sapply(rownames(aln), function(x) strsplit(x, '_')[[1]][2]))
y <- as.numeric(sapply(rownames(aln), function(x) strsplit(x, '_')[[1]][3]))

sim_coords <- cbind(x, y)
rownames(sim_coords) <- tr$tip.label
fit1 <- fitContinuous(tr, sim_coords, model = 'lambda')


beast_estimates <- matrix(c(7.546, -4.398, -4.398, 0.105), 2, 2) / 100

#fit1 <- ace(sim_coords, tr, method = 'ML')

#ace(sim_single[1, ], tr, method = 'ML')
