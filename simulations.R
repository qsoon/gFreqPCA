# packages may be redundant..
library(gasper)
library(GNAR)
library(ggplot2)
library(cccd)
library(igraph)
library(igraphdata)
library(MASS)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)
library(expm)
library(Hmisc)

source("utils.R")

# window filterbank mother function for WGF-based estimation
g <- function(lambda, sigma.sq, m, tau){
  return(exp(-(lambda-m*tau)^2 / sigma.sq))
}

#############################
### 01. graph preparation ###
#############################

###########################
### karate club network ###
###########################
data(karate)
V(karate)$name
E(karate)
edge.attributes(karate)
degree(karate)


N.karate <- length(V(karate)$name)
edge.wt <- igraph::as_data_frame(karate, what="edges")
for(i in 1:nrow(edge.wt)){
  edge.wt[i,1] <- which(V(karate)$name == edge.wt[i,1])
  edge.wt[i,2] <- which(V(karate)$name == edge.wt[i,2])
}
edge.wt <- sapply(edge.wt, as.numeric)
wmat.karate <- matrix(0, nrow=gorder(karate), ncol=gorder(karate))
# colnames(wmat.karate) <- V(karate)$name
# rownames(wmat.karate) <- V(karate)$name

colnames(wmat.karate) <- 1:N.karate
rownames(wmat.karate) <- 1:N.karate

for(i in 1:nrow(edge.wt)){
  wmat.karate[edge.wt[i,1], edge.wt[i,2]] <- edge.wt[i,3]
  wmat.karate[edge.wt[i,2], edge.wt[i,1]] <- edge.wt[i,3]
}

sp.wmat.karate <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat.karate <- rbind(sp.wmat.karate, c(edge.wt[i,1], edge.wt[i,2], 
                                            wmat.karate[edge.wt[i,1], edge.wt[i,2]]))
}

zachary.karate <- list()
set.seed(1)
layout.karate <- layout.fruchterman.reingold(karate)
layout.karate <- layout.norm(layout.karate, -1, 1, -1, 1)
zachary.karate$xy <- data.frame(x=layout.karate[,1],
                                y=layout.karate[,2])

zachary.karate$A <- wmat.karate
zachary.karate$sA <- sp.wmat.karate

plot_graph(zachary.karate)

L.karate <- gasper::laplacian_mat(wmat.karate)
val <- eigensort(L.karate)
lmax.karate <- max(val$evalues)
evalues.karate <- val$evalues
evectors.karate <- val$evectors

eigenres.karate <- eigen(L.karate)

######################
## US sensor network##
######################
data_ustemp <- readMat("Data/stationary/US_Hourly_2010_August_1st.mat") # hourly data on 2010.08.01
UShourlytemp <- list()
UShourlytemp$xy <- cbind(data_ustemp$x, data_ustemp$y)
rownames(UShourlytemp$xy) <- 1:nrow(UShourlytemp$xy) 
UShourlytemp$f1 <- data_ustemp$signal[,15]
UShourlytemp$f2 <- data_ustemp$signal[,3]

distmat.htemp <- distm(UShourlytemp$xy, fun = distHaversine) / 1000
A.htemp <- c()
for(i in 1:(nrow(distmat.htemp)-1)){
  for(j in (i+1):ncol(distmat.htemp)){
    val <- distmat.htemp[i,j]
    A.htemp <- rbind(A.htemp, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.htemp, k=7), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.htemp[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=length(UShourlytemp$f1), ncol=length(UShourlytemp$f1))

colnames(wmat) <- 1:length(UShourlytemp$f1)  
rownames(wmat) <- 1:length(UShourlytemp$f1)

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}

# weight matrix
UShourlytemp$A <- wmat

# sparse weight matrix
UShourlytemp$sA <- sp.wmat

UShourlytemp$dist <- distmat.htemp
UShourlytemp$sdist <- A.htemp

# visualize
plot_graph(UShourlytemp)

L.ush <- gasper::laplacian_mat(UShourlytemp$A)
N.ush <- nrow(UShourlytemp$xy)
val1 <- eigensort(L.ush)
evalues.ush <- val1$evalues
evectors.ush <- val1$evectors
# largest eigenvalue
lmax.ush <- max(evalues.ush)



##############
## gFreqPCA ##
##############

#########################################
### 02. detection of common frequency ###
#########################################

# Karate network
X.karate <- list()


coef10.karate <- c(1, 2.5, 3.5,
                   0, 0, 0, 
                   2.1, 1.4, 2.5,
                   0, 0, 0)

coef20.karate <- c(0,0,0,
                   2, 1.7, 3.2,
                   0.9, 2, 2.2,
                   0,0,0)

set.seed(1)
for(i in 1:3){
  X.karate[[i]] <- coef10.karate[i]*evectors.karate[,10] + coef20.karate[i]*evectors.karate[,20] + rnorm(N.karate, 0, 0.5)
}


set.seed(100)
for(i in 4:6){
  X.karate[[i]] <- coef10.karate[i]*evectors.karate[,10] + coef20.karate[i]*evectors.karate[,20] + rnorm(N.karate, 0, 0.5)
}

set.seed(1000)
for(i in 7:9){
  X.karate[[i]] <- coef10.karate[i]*evectors.karate[,10] + coef20.karate[i]*evectors.karate[,20] + rnorm(N.karate, 0, 0.5)
}

set.seed(10000)
for(i in 10:12){
  X.karate[[i]] <- coef10.karate[i]*evectors.karate[,10] + coef20.karate[i]*evectors.karate[,20] + rnorm(N.karate, 0, 0.5)
}


res.karate.random_window <- GFPCA(X=X.karate, S=L.karate, M=50, sigma=0.5, method="random")

# scree plot?
plot(colSums(res.karate.random_window$tau.hat) / sum(res.karate.random_window$tau.hat), type="o")

round(cumsum(colSums(res.karate.random_window$tau.hat) / sum(res.karate.random_window$tau.hat)),2)


res.karate.random_window <- GFPCA(X=X.karate, S=L.karate, M=50, sigma=0.5, q=2, method="random")

# spectral envelope
par(mfrow=c(1,1))
plot(res.karate.random_window$tau.hat[,1], type="l")

par(mfcol=c(3,4), mar=c(5,4.5,4,2)+0.1)
for(i in 1:12){
  plot(psd.graph(X.karate[[i]], S=L.karate, M=50, sigma=0.5, method="random"), type="l",
       xlab="Graph frequency index", ylab="GPSD", cex.lab=1.5)
}

# scaling
par(mfrow=c(1,3), mar=c(5,5,4,2)+0.1, mgp=c(3,1,0))
# graph spectral envelope
plot(res.karate.random_window$tau.hat[,1], type="l", main="Graph spectral envelope",
     xlab="Graph frequency index", ylab=expression(tau[1]), cex.lab=2, cex.main=2, cex.axis=1.2, lwd=1.5)
abline(v=c(10,20), col="red", lty=2)
par(mar=c(2,5,4,2)+0.1)
# amplitude estimation
barplot(rbind(coef10.karate / norm(coef10.karate, type="2"), -res.karate.random_window$H.hat[10,1,]),
        beside = TRUE, col = c("black", "gray"), main=expression(lambda[10]^"KN"), cex.main=2,
        ylim = c(-0.1, 1))
legend("topright", legend = expression("True", "Estimate"),
       fill = c("black", "gray"),
       cex = 1.7,          # Adjust text size
       box.lwd = 1,       # Adjust the line width of the legend box
       box.col = "black")
title(xlab=list("Dimension index", cex=2), mgp=c(0.1,1,0))
title(ylab=list("Normalized amplitude", cex=2), mgp=c(3,1,0))

par(mar=c(2,5,4,2)+0.1)
barplot(rbind(coef20.karate / norm(coef20.karate, type="2"), res.karate.random_window$H.hat[20,1,]),
        beside = TRUE, col = c("black", "gray"), main=expression(lambda[20]^"KN"), cex.main=2,
        ylim = c(-0.1, 1))
legend("topright", legend = expression("True", "Estimate"),
       fill = c("black", "gray"),
       cex = 1.7,          # Adjust text size
       box.lwd = 1,       # Adjust the line width of the legend box
       box.col = "black")
title(xlab=list("Dimension index", cex=2), mgp=c(0.1,1,0))
title(ylab=list("Normalized amplitude", cex=2), mgp=c(3,1,0))

# reconstruction error
lapply(mapply("-", res.karate.random_window$X, res.karate.random_window$X.hat, SIMPLIFY = FALSE), mean)
plot(res.karate.random_window$X[[1]], type="l")
lines(res.karate.random_window$X.hat[[1]], col="red")


g.data.karate.list <- list()
g.data.karate.list[[1]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[1]], value="Value", ratio=0.6,
                                              min=min(X.karate[[1]]), max=max(X.karate[[1]]), mg=c(4,4,4,6), title=expression(italic(X)[1]), main.title.size = 20)
g.data.karate.list[[5]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[2]], value="Value", ratio=0.6,
                                              min=min(X.karate[[2]]), max=max(X.karate[[2]]), mg=c(4,4,4,6), title=expression(italic(X)[2]), main.title.size = 20)
g.data.karate.list[[9]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[3]], value="Value", ratio=0.6,
                                              min=min(X.karate[[3]]), max=max(X.karate[[3]]), mg=c(4,4,4,6), title=expression(italic(X)[3]), main.title.size = 20)
g.data.karate.list[[2]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[4]], value="Value", ratio=0.6,
                                              min=min(X.karate[[4]]), max=max(X.karate[[4]]), mg=c(4,4,4,6), title=expression(italic(X)[4]), main.title.size = 20)
g.data.karate.list[[6]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[5]], value="Value", ratio=0.6,
                                              min=min(X.karate[[5]]), max=max(X.karate[[5]]), mg=c(4,4,4,6), title=expression(italic(X)[5]), main.title.size = 20)
g.data.karate.list[[10]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[6]], value="Value", ratio=0.6,
                                              min=min(X.karate[[6]]), max=max(X.karate[[6]]), mg=c(4,4,4,6), title=expression(italic(X)[6]), main.title.size = 20)
g.data.karate.list[[3]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[7]], value="Value", ratio=0.6,
                                              min=min(X.karate[[7]]), max=max(X.karate[[7]]), mg=c(4,4,4,6), title=expression(italic(X)[7]), main.title.size = 20)
g.data.karate.list[[7]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[8]], value="Value", ratio=0.6,
                                              min=min(X.karate[[8]]), max=max(X.karate[[8]]), mg=c(4,4,4,6), title=expression(italic(X)[8]), main.title.size = 20)
g.data.karate.list[[11]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[9]], value="Value", ratio=0.6,
                                              min=min(X.karate[[9]]), max=max(X.karate[[9]]), mg=c(4,4,4,6), title=expression(italic(X)[9]), main.title.size = 20)
g.data.karate.list[[4]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[10]], value="Value", ratio=0.6,
                                              min=min(X.karate[[10]]), max=max(X.karate[[10]]), mg=c(4,4,4,6), title=expression(italic(X)[10]), main.title.size = 20)
g.data.karate.list[[8]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[11]], value="Value", ratio=0.6,
                                              min=min(X.karate[[11]]), max=max(X.karate[[11]]), mg=c(4,4,4,6), title=expression(italic(X)[11]), main.title.size = 20)
g.data.karate.list[[12]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[12]], value="Value", ratio=0.6,
                                              min=min(X.karate[[12]]), max=max(X.karate[[12]]), mg=c(4,4,4,6), title=expression(italic(X)[12]), main.title.size = 20)

# for(i in 1:12){
#   g.data.karate.list[[i]] <- plot_graph_custom4(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[i]], value="Passengers", ratio=0.6,
#                                                min=min(X.karate[[i]]), max=max(X.karate[[i]]), mg=c(4,4,4,4), title=expression(italic(X)[i]), main.title.size = 20)
# }
grid.arrange(grobs=g.data.karate.list, nrow=3)

g.res.karate.list <- list()
g.res.karate.list[[1]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[1]]-res.karate.random_window$X.hat[[1]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[1] - hat(italic(X))[1]), main.title.size = 20)
g.res.karate.list[[5]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[2]]-res.karate.random_window$X.hat[[2]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[2] - hat(italic(X))[2]), main.title.size = 20)
g.res.karate.list[[9]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[3]]-res.karate.random_window$X.hat[[3]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[3] - hat(italic(X))[3]), main.title.size = 20)
g.res.karate.list[[2]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[4]]-res.karate.random_window$X.hat[[4]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[4] - hat(italic(X))[4]), main.title.size = 20)
g.res.karate.list[[6]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[5]]-res.karate.random_window$X.hat[[5]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[5] - hat(italic(X))[5]), main.title.size = 20)
g.res.karate.list[[10]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[6]]-res.karate.random_window$X.hat[[6]], value="Value", ratio=0.6,
                                               min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[6] - hat(italic(X))[6]), main.title.size = 20)
g.res.karate.list[[3]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[7]]-res.karate.random_window$X.hat[[7]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[7] - hat(italic(X))[7]), main.title.size = 20)
g.res.karate.list[[7]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[8]]-res.karate.random_window$X.hat[[8]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[8] - hat(italic(X))[8]), main.title.size = 20)
g.res.karate.list[[11]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[9]]-res.karate.random_window$X.hat[[9]], value="Value", ratio=0.6,
                                               min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[9] - hat(italic(X))[9]), main.title.size = 20)
g.res.karate.list[[4]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[10]]-res.karate.random_window$X.hat[[10]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[10] - hat(italic(X))[10]), main.title.size = 20)
g.res.karate.list[[8]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[11]]-res.karate.random_window$X.hat[[11]], value="Value", ratio=0.6,
                                              min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[11] - hat(italic(X))[11]), main.title.size = 20)
g.res.karate.list[[12]] <- plot_graph_custom5(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[12]]-res.karate.random_window$X.hat[[12]], value="Value", ratio=0.6,
                                               min=-0.072, max=0.069, mg=c(4,4,4,6), title=expression(italic(X)[12] - hat(italic(X))[12]), main.title.size = 20)

# for(i in 1:12){
#   g.data.karate.list[[i]] <- plot_graph_custom4(zachary.karate, e.size=1.3, v.size=6, vertex_color = X.karate[[i]], value="Passengers", ratio=0.6,
#                                                min=min(X.karate[[i]]), max=max(X.karate[[i]]), mg=c(4,4,4,4), title=expression(italic(X)[i]), main.title.size = 20)
# }
grid.arrange(grobs=g.res.karate.list, nrow=3)

for(i in 1:12){
  print(range(X.karate[[i]] - res.karate.random_window$X.hat[[i]]))
}
range(X.karate[[1]])


# US sensor network
X.ush <- list()
coef50.ush <- c(3, 1.5, 2, 
                0, 0, 0, 
                0, 0, 0, 
                2, 0, 0)

coef100.ush <- c(0, 0, 0,
                 2, 4, 3, 
                 0, 0, 0, 
                 4, 3, 0)

coef150.ush <- c(0, 0, 0,
                 0, 0, 0, 
                 5, 2, 1.5, 
                 0, 2.5, 0)


set.seed(1)
for(i in 1:3){
  X.ush[[i]] <- coef50.ush[i]*evectors.ush[,50] + coef100.ush[i]*evectors.ush[,100] + 
    coef150.ush[i]*evectors.ush[,150] + rnorm(N.ush, 0, 0.5)
}

set.seed(100)
for(i in 4:6){
  X.ush[[i]] <- coef50.ush[i]*evectors.ush[,50] + coef100.ush[i]*evectors.ush[,100] + 
    coef150.ush[i]*evectors.ush[,150] + rnorm(N.ush, 0, 0.5)
}

set.seed(1000)
for(i in 7:9){
  X.ush[[i]] <- coef50.ush[i]*evectors.ush[,50] + coef100.ush[i]*evectors.ush[,100] + 
    coef150.ush[i]*evectors.ush[,150] + rnorm(N.ush, 0, 0.5)
}

set.seed(10000)
for(i in 10:11){
  X.ush[[i]] <- coef50.ush[i]*evectors.ush[,50] + coef100.ush[i]*evectors.ush[,100] + 
    coef150.ush[i]*evectors.ush[,150] + rnorm(N.ush, 0, 0.5)
}

set.seed(100000)
for(i in 12:12){
  X.ush[[i]] <- 2*rnorm(N.ush, 0, 0.5)
}

res.ush.random_window <- GFPCA(X=X.ush, S=L.ush, M=50, sigma=0.5, method="random")

# scree plot?
plot(colSums(res.ush.random_window$tau.hat) / sum(res.ush.random_window$tau.hat), type="o")
round(cumsum(colSums(res.ush.random_window$tau.hat) / sum(res.ush.random_window$tau.hat)),2)
res.ush.random_window <- GFPCA(X=X.ush, S=L.ush, M=50, sigma=0.5, q=4, method="random")

# spectral envelope
plot(res.ush.random_window$tau.hat[,1], type="l")

par(mfcol=c(3,4), mar=c(5,4.5,4,2)+0.1)
for(i in 1:12){
  plot(psd.graph(X.ush[[i]], S=L.ush, M=50, sigma=0.5, method="random"), type="l",
       xlab="Graph frequency index", ylab="GPSD", cex.lab=1.5)
}


# scaling
par(mfrow=c(2,2), mar=c(5,5,4,2)+0.1, mgp=c(3,1,0))
# graph spectral envelope
plot(res.ush.random_window$tau.hat[,1], type="l", main="Graph spectral envelope",
     xlab="Graph frequency index", ylab=expression(tau[1]), cex.lab=2, cex.main=2, cex.axis=1.2, lwd=1.5)
abline(v=c(50,100,150), col="red", lty=2)
par(mar=c(3,5,4,2)+0.1)
barplot(rbind(coef50.ush / norm(coef50.ush, type="2"), -res.ush.random_window$H.hat[50,1,]),
        beside = TRUE, col = c("black", "gray"), main=expression(lambda[50]^"US"), cex.main=2,
        ylim = c(-0.1, 1))
legend("topright", legend = expression("True", "Estimate"),
       fill = c("black", "gray"),
       cex = 1.4,          # Adjust text size
       box.lwd = 1,       # Adjust the line width of the legend box
       box.col = "black")
title(xlab=list("Dimension index", cex=2), mgp=c(1,1,0))
title(ylab=list("Normalized amplitude", cex=2), mgp=c(3,1,0))

par(mar=c(3,5,4,2)+0.1)
# amplitude estimation
barplot(rbind(coef100.ush / norm(coef100.ush, type="2"), -res.ush.random_window$H.hat[100,1,]),
        beside = TRUE, col = c("black", "gray"), main=expression(lambda[100]^"US"), cex.main=2,
        ylim = c(-0.1, 1))
legend("topright", legend = expression("True", "Estimate"),
       fill = c("black", "gray"),
       cex = 1.4,          # Adjust text size
       box.lwd = 1,       # Adjust the line width of the legend box
       box.col = "black")
title(xlab=list("Dimension index", cex=2), mgp=c(1,1,0))
title(ylab=list("Normalized amplitude", cex=2), mgp=c(3,1,0))

par(mar=c(3,5,4,2)+0.1)
barplot(rbind(coef150.ush / norm(coef150.ush, type="2"), res.ush.random_window$H.hat[150,1,]),
        beside = TRUE, col = c("black", "gray"), main=expression(lambda[150]^"US"), cex.main=2,
        ylim = c(-0.1, 1))
legend("topright", legend = expression("True", "Estimate"),
       fill = c("black", "gray"),
       cex = 1.4,          # Adjust text size
       box.lwd = 1,       # Adjust the line width of the legend box
       box.col = "black")
title(xlab=list("Dimension index", cex=2), mgp=c(1,1,0))
title(ylab=list("Normalized amplitude", cex=2), mgp=c(3,1,0))

# reconstruction error
lapply(mapply("-", res.ush.random_window$X, res.ush.random_window$X.hat, SIMPLIFY = FALSE), mean)
plot(res.ush.random_window$X[[1]], type="l")
lines(res.ush.random_window$X.hat[[1]], col="red")



# two plots together
g1 <- plot_graph_custom3(zachary.karate, e.size=1.3, v.size=3, vertex_color = "black", 
                         value="value", ratio=0.7, signal=FALSE)
g2 <- plot_graph_custom3(UShourlytemp, e.size=1.3, v.size=3, vertex_color = "black", 
                         value="value", ratio=0.7, signal=FALSE)
grid.arrange(g1,g2, nrow=1)


# two scree plots
par(mfrow=c(1,2), mar=c(5,4,4,2)+0.1)
plot(colSums(res.karate.random_window$tau.hat) / sum(res.karate.random_window$tau.hat), type="o", xlab="q", ylab="", main="Karate club network", cex.lab=1.3)
plot(colSums(res.ush.random_window$tau.hat) / sum(res.ush.random_window$tau.hat), type="o", xlab="q", ylab="", main="US sensor network", cex.lab=1.3)

