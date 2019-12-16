library(MASS)
library(NELSI)
library(phytools)
library(phangorn)
library(geiger)

tr <- read.tree(text = '((x1:1.5,x2:1.5):1.5, (x3:0.25,x4:0.25):2.75);')
plot(tr)
edgelabels(text = tr$edge.length)
is.rooted(tr)
is.ultrametric(tr)
vcvPhylo(tr, anc.nodes = T)

# Single trait
C = vcv(tr)
R = 0.1 # variance of trait
V = kronecker(R, C)
a = rep(0, 4) # expected mean of trait
sim_single <- mvrnorm(n = 100, mu = a, Sigma = V)
par(mfrow = c(3, 1), mar = c(2, 2, 1, 1))
plot(tr)
nodelabels()
plot(sim_single[, 1:2])
plot(sim_single[, 3:4])

nodelabels()
plot(sim_single[, 1:2])
plot(sim_single[, 3:4])

# Two traits with covariance
C <- vcv.phylo(tr)
R <- matrix(c(1, -0.99, -0.99, 2), 2)
a <- rep(0, 8) # two traits for 4 species, with means of 0
V = kronecker(R, C)
colnames(V) <- c(paste0('x', 1:4,'_1'), paste0('x', 1:4,'_2'))
rownames(V) <- colnames(V)
sim_multi <- mvrnorm(n = 100, mu = rep(0, 8), Sigma = V)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(tr)
nodelabels()
plot(sim_multi[, 'x1_1'], sim_multi[, 'x1_2'])
plot(sim_multi[, 'x3_1'], sim_multi[, 'x4_1'])
plot(sim_multi[, 'x1_1'], sim_multi[, 'x2_1'])
#plot(sim_multi[, 1],  sim_1[, 2])
##
#         1     2    3    4     5    6      7      8
#       x1_1  x2_1 x3_1 x4_1  x1_2  x2_2   x3_2   x4_2
# x1_1 3.000 1.500 0.00 0.00 2.970 1.485 0.0000 0.0000
# x2_1 1.500 3.000 0.00 0.00 1.485 2.970 0.0000 0.0000

# Second method
e <- eigen(R, EISPACK = T)
f <- eigen(C, EISPACK = T)
A <- f$vectors %*% diag(sqrt(f$values))
B <- e$vectors %*% diag(sqrt(e$values))
sim_multi_m <- matrix(NA, 100, 8)
colnames(sim_multi_m) <- colnames(sim_multi)
for(i in 1:nrow(sim_multi_m)){
    X <- matrix(rnorm(8), 4)
    #t(X %*% chol(R)) %*% chol(C)
    sim_multi_m[i, ] <- as.vector(A %*% t(B %*% t(X)))
}
plot(tr)
plot(sim_multi_m[, 'x1_1'], sim_multi_m[, 'x1_2'])
plot(sim_multi_m[, 'x3_1'], sim_multi_m[, 'x4_1'])
plot(sim_multi_m[, 'x1_1'], sim_multi_m[, 'x2_1'])

par(mfrow = c(2, 2))
hist(sim_multi_m[, 'x1_1'])
hist(sim_multi_m[, 'x3_1'])
hist(sim_multi_m[, 'x1_2'])
hist(sim_multi_m[, 'x3_2'])

var(sim_multi_m[, 'x1_1'])
var(sim_multi_m[, 'x3_1'])
var(sim_multi_m[, 'x1_2'])
var(sim_multi_m[, 'x3_2'])



## Two simulations on larger trees with different variances.

tr <- rcoal(100)
tr$edge.length <- tr$edge.length * 10
max(allnode.times(tr))

par(mfrow = c(1, 1))
plot(tr)
nts <- round(intnode.times(tr), 2)
nodelabels(nts)

#
C <- vcvPhylo(tr, anc.nodes = T)
R <- matrix(c(0.5, 0.01, 0.01, 1), 2)
a <- rep(0, 2*(length(tr$tip.label) + tr$Nnode - 1))
V <- kronecker(R, C)
colnames(V) <- c(paste0(colnames(C), '_1'), paste0(colnames(C), '_2'))
rownames(V) <- colnames(V)

sim_trait <- round(mvrnorm(n = 1, mu = a, Sigma = V), 3)

tr$tip.label <- paste0(tr$tip.label, '_', sim_trait[1:length(tr$tip.label)], '_',
       sim_trait[(length(tr$tip.label)+1):(2*length(tr$tip.label))])

# Can we also get the values at internal nodes?
plot(tr, show.tip.label = F)
nodelabels()
tiplabels()


tip_labels <- gsub('_.+', '', tr$tip.label)
# Get start and end points, branch length, distance
                                        # requires tr, sim_trait

geo.data.matrix <- get.geo.data.matrix(tr, sim_trait, tip_labels)

plot(geo.data.matrix[, 1:2], type = 'n', xlim = range(geo.data.matrix[, c('x1', 'x2')]),
     ylim = range(geo.data.matrix[, c('y1', 'y2')]))
for(i in 1:nrow(geo.data.matrix)){
    lines(geo.data.matrix[i, c('x1', 'x2')], geo.data.matrix[i, c('y1', 'y2')], col = i, lwd = 2)
    points(geo.data.matrix[i, c('x1', 'x2')], geo.data.matrix[i, c('y1', 'y2')], pch = 20, col = i)
}




#############################################
##### Simulate a data set with some parameters
library(TreeSim)
tr <- sim.bdsky.stt(n = 100, lambdasky = 2, deathsky = 1, sampprobsky = 0.1,
                    timesky = 0)[[1]]
par(mfrow = c(1, 1))
ntimes <- round(allnode.times(tr, tipsonly = T), 2)
tr$tip.label <- paste0(tr$tip.label, '@', ntimes, '@')
plot(tr)
max(allnode.times(tr))
#
par(mfrow = c(1, 1))
plot(tr, show.tip.label = F)
nts <- round(intnode.times(tr), 2)
nodelabels(nts)
#
C <- vcvPhylo(tr, anc.nodes = T)
R <- matrix(c(1, 0.9, 0.9, 1), 2)
a <- rep(0, 2*(length(tr$tip.label) + tr$Nnode - 1))
V <- kronecker(R, C)
colnames(V) <- c(paste0(colnames(C), '_1'), paste0(colnames(C), '_2'))
rownames(V) <- colnames(V)
#
sim_trait <- round(mvrnorm(n = 1, mu = a, Sigma = V), 3)
tr$tip.label <- paste0(tr$tip.label, '_', sim_trait[1:length(tr$tip.label)], '_',
       sim_trait[(length(tr$tip.label)+1):(2*length(tr$tip.label))])
# sanity check plot
tip_labels <- gsub('_.+', '', tr$tip.label)
#
geo.data.matrix <- get.geo.data.matrix(tr, sim_trait, tip_labels)
#

par(mfrow = c(1, 3))
plot(geo.data.matrix[, 1:2], type = 'n', xlim = range(geo.data.matrix[, c('x1', 'x2')]),
     ylim = range(geo.data.matrix[, c('y1', 'y2')]))
for(i in 1:nrow(geo.data.matrix)){
    lines(geo.data.matrix[i, c('x1', 'x2')], geo.data.matrix[i, c('y1', 'y2')], col = i, lwd = 2)
    points(geo.data.matrix[i, c('x1', 'x2')], geo.data.matrix[i, c('y1', 'y2')], pch = 20, col = i)
}
hist(geo.data.matrix[, 'distance'] / geo.data.matrix[, 'edge_length'])
plot(tr, show.tip.label = F, edge.col = 'darkgrey', lwd = 2)

#
diffusion_rate <- sum(geo.data.matrix[, 'distance']) / sum(geo.data.matrix[, 'edge_length'])
#
phylogram <- tr
phylogram$edge.length <- phylogram$edge.length * 1E-2
aln <- as.DNAbin(simSeq(phylogram, l = 5000))
length(seg.sites(aln))
write.dna(aln, file = 'test_phylogeo_3.fasta', format = 'fasta', nbcol = -1, colsep = '')

#hist(geo.data.matrix[, 'distance'] / geo.data.matrix[, 'edge_length'])

median(geo.data.matrix[, 'distance'] / geo.data.matrix[, 'edge_length'])


log_cauchy <- tail(read.table('geo_testing/test_phylogeo_4_cauchy.log', head = T), 1000)

hist(log_cauchy$LOCATION.diffusionRate * log_cauchy$treeLength)


diffusion_rate
max(allnode.times(tr))
sum(tr$edge.length)
write.tree(tr, file = 'geo_testing/true_tree.tre')



## Been testing up to here....




median(geo.data.matrix[, 'distance'] / geo.data.matrix[, 'edge_length'])


no_root <- geo.data.matrix[!geo.data.matrix[, 'node_from'] == length(tr$tip.label)+1, ]

mean(no_root[, 'distance'] / no_root[, 'edge_length'])
sum(no_root[, 'distance']) / sum(no_root[, 'edge_length'])



### Map BEAST results onto true tree and compare branches

# Make tree with distance in branch lenghts
dist_tree <- tr
dist_tree$edge.length <- geo.data.matrix[, 'distance']
par(mfrow = c(1, 1))
plot(dist_tree, show.tip.label = F)

write.nexus(dist_tree, file = 'testing_phylogeo_true.tree')



## Estimate values using ML

BMlk <- function(C, inv.C, sigmasq, root.state, data) {
    N <- length(data); # the number of tips
    EX <- rep(root.state, N) # creates a vector of the expected trait value - which under BM is the root state
    V <- C * sigmasq; # multiply the entries in C by the BM rate
    inv.V <- inv.C * sigmasq ^-1; # do the same for the inverted matrix using the inverse of the rate
    lnlNum<- -0.5*(data - EX) %*% inv.V %*% (data - EX)
    lnlDen<- log(sqrt((2*pi)^N*det(V)))
    L<-lnlNum-lnlDen
    return(L);
}

trait1 <- sapply(tr$tip.label, function(x) as.numeric(strsplit(x, split = '_')[[1]][2]))
trait2 <- sapply(tr$tip.label, function(x) as.numeric(strsplit(x, split = '_')[[1]][3]))

c <- vcv.phylo(tr)
inv.c <- solve(c)

BMlk(c, inv.c, sigmasq = 1, root.state = 0, data = trait1)
BMlk(c, inv.c, sigmasq = 3, root.state = 0, data = trait1)

ml_fit_trait1 <- fitContinuous(tr, trait1)
ml_fit_trait2 <- fitContinuous(tr, trait2)

ml_fit_trait1$opt$sigsq * length(tr$tip.label) / (length(tr$tip.label) - 1)
ml_fit_trait2$opt$sigsq * length(tr$tip.label) / (length(tr$tip.label) - 1)


ace(



## Estimate values using BEAST


## Recover correct values
precision <- sum(tr$edge.length)*matrix(c(6.7E-2, -9.6E-4,  -9.6E-4, 0.137), nrow = 2)
#precision <- matrix(c(15.54, 0.119, 0.119, 4.46), nrow = 2)
variance <- solve(precision)
correlation <- variance[1,2] / (sqrt(variance[1,1] * variance[2,2]))

variance
correlation

variance[1,1] / variance[2,2]

variance[2,2] / (sqrt(variance[1,1] * variance[2,2]))

sum(tr$edge.length) * (1 / matrix(c(15.54, 0.119, 0.119, 4.46), nrow = 2))

precision * sum(tr$edge.length)

# Fit a BM for tree
