# may be redundant 

library(RColorBrewer)
library(colourvalues)
library(grDevices)
# library(SDMTools)
library(network)
library(ggraph)
library(tidygraph)
library(dplyr)
library(gasper)
library(readxl)
library(forecast)
library(ggfortify)
library(Metrics)
library(GNAR)
library(DTWBI)
# library(vars)
library(geosphere)
# library(xlsx)
library(scales)
library(igraph)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)
library(expm)
library(Hmisc)
library(cccd)
library(zoo)

source("utils.R")


trade2023 <- read.csv("BACI_HS92_V202501/BACI_HS92_Y2023_V202501.csv", header=TRUE)
trade2023$q <- as.numeric(trade2023$q)
trade2023 <- na.omit(trade2023)

trade2023$value <- trade2023$v * trade2023$q


identical(sort(unique(trade2023$i)), sort(unique(trade2023$j)))
head(trade2023)


countrycode <- read.csv("BACI_HS92_V202501/country_codes_V202501b.csv", header=TRUE)
countrycode_from_ind <- read.csv("economic/Indicators/agriculture.csv", skip=4)

countrycode3digits <- sort(intersect(countrycode$country_iso3, countrycode_from_ind$Country.Code))
countrylist_from_trade <- countrycode[countrycode$country_iso3 %in% countrycode3digits,]$country_name
countrylist_from_ind <- countrycode_from_ind[countrycode_from_ind$Country.Code %in% countrycode3digits,]$Country.Name

country_name_code_tbl <-cbind(countrylist_from_ind, countrycode3digits)
colnames(country_name_code_tbl) <- c("name", "iso3")
country_name_code_tbl <- as.data.frame(country_name_code_tbl)
country_name_code_tbl$lat <- 0
country_name_code_tbl$lon <- 0

coords <- read.csv("country_coord.csv", header=TRUE)
for(i in 1:nrow(country_name_code_tbl)){
  if(country_name_code_tbl$iso3[i]=="CUW"){
    country_name_code_tbl[i,"lat"] <- 12.2135
    country_name_code_tbl[i,"lon"] <- -68.9496
  } else if(country_name_code_tbl$iso3[i]=="SXM"){
    country_name_code_tbl[i,"lat"] <- 18.0347
    country_name_code_tbl[i,"lon"] <- -63.0681
  } else{
    country_name_code_tbl[i,"lat"] <-  coords[which(coords$Alpha.3.code==country_name_code_tbl$iso3[i]),5]
    country_name_code_tbl[i,"lon"] <-  coords[which(coords$Alpha.3.code==country_name_code_tbl$iso3[i]),6]
  }
}

countrycode <- countrycode[-c(18,81,203),] # remove Belgium-Luxembourg (...1998) / Fed. Rep. of Germany (...1990) / Sudan (...2011)
countrycode_used <- countrycode[countrycode$country_iso3 %in% country_name_code_tbl[,2],]


trade2023_used <- trade2023[(trade2023$i %in% countrycode_used$country_code) & (trade2023$j %in% countrycode_used$country_code),]



trade2023_final <- trade2023_used %>%
  # Create new variables i_new and j_new to treat (i, j) and (j, i) as the same pair
  mutate(i_new = pmin(i, j), j_new = pmax(i, j)) %>%
  # Group by these new variables
  group_by(i_new, j_new) %>%
  # Summarize by summing the total
  summarise(
    total = sum(v) / ifelse(sum(i != i_new | j != j_new) > 0, 2, 1),  # Divide by 2 if both directions exist
    .groups = 'drop'
  ) %>%
  # Rename columns back to original names if desired
  rename(i = i_new, j = j_new)


trade2023_final <- trade2023_final %>% arrange(desc(total))
trade2023_final <- trade2023_final[1:round(0.03*nrow(trade2023_final)),]

trade2023_final$total <- trade2023_final$total/10^6
trade2023_final$total_log <- log(trade2023_final$total)

countrycode_final <- countrycode_used[countrycode_used$country_code %in% union(trade2023_final$i, trade2023_final$j),]
countrycode_final_sort <- countrycode_final %>% arrange(country_iso3)

ITW <- list()
# location
country_name_code_tbl_final <- country_name_code_tbl[(country_name_code_tbl$iso3 %in% countrycode_final$country_iso3),]
ITW$xy <- country_name_code_tbl_final[,c("lon","lat")]
rownames(ITW$xy) <- country_name_code_tbl_final$iso3

# weight matrix
e.weight <- matrix(0,nrow=nrow(country_name_code_tbl_final), ncol=nrow(country_name_code_tbl_final))
colnames(e.weight) <- country_name_code_tbl_final$iso3
rownames(e.weight) <- country_name_code_tbl_final$iso3

e.sp.weight <- NULL
for(k in 1:nrow(trade2023_final)){
  i <- which(countrycode_final_sort$country_code==as.numeric(trade2023_final[k,1]))
  j <- which(countrycode_final_sort$country_code==as.numeric(trade2023_final[k,2]))
  e.weight[i,j] <- as.numeric(trade2023_final[k,4])
  e.weight[j,i] <- as.numeric(trade2023_final[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(i,j,as.numeric(trade2023_final[k,4])))
  e.sp.weight <- rbind(e.sp.weight, c(j,i,as.numeric(trade2023_final[k,4])))
}


ITW$A <- e.weight

# sparse weight matrix
ITW$sA <- e.sp.weight[nrow(e.sp.weight):1,]

L.ITW <- laplacian_mat(ITW$A) # laplacian matrix
val1 <- eigensort(L.ITW)
evalues.ITW <- val1$evalues
evectors.ITW <- val1$evectors
# largest eigenvalue
lmax.ITW <- max(evalues.ITW)

N.ITW <- nrow(L.ITW)


plot_graph(ITW)



# load economic data
files <- list.files(path = "./economic/Indicators", pattern = ".csv")
p.ITW <- length(files)
R.ITW <- 34 # 1990~2023
X.ITW <- array(0, c(p.ITW,N.ITW,R.ITW))

tmp2 <- c()
for(i in 1:p.ITW){
  print(paste("i=",i))
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  print((as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                              paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  tmp2 <- c(tmp2, (as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                                        paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # for(l in 1:N.ITW){
  #   
  # }
  # X.ITW[i,,] <- as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
  #     paste("X", 1990:2023, sep="")]) na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
}

variable_remove <- c(4,5,7,17,19,26,27,28,33,34,35,36,37) 

tmp3 <- c()
for(i in 1:p.ITW){
  if(i %in% variable_remove){
    next
  }
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  print((as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                              paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  tmp3 <- c(tmp3, (as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                                        paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # for(l in 1:N.ITW){
  #   
  # }
  # X.ITW[i,,] <- as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
  #     paste("X", 1990:2023, sep="")]) na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
}

country_name_code_tbl_final[unique(tmp3),]
table(tmp3)

country_remove <- unique(tmp3)

ITW2 <- list()
ITW2$xy <- ITW$xy[-country_remove,]

countrycode_final_sort2 <- countrycode_final_sort[-country_remove, ]
rownames(countrycode_final_sort2) <- 1:N.ITW2
trade2023_final2 <- trade2023_final[(trade2023_final$i %in% countrycode_final_sort2$country_code),]
trade2023_final2 <- trade2023_final2[(trade2023_final2$j %in% countrycode_final_sort2$country_code),]
ITW2$A <- ITW$A[-country_remove, -country_remove] 
e.sp.weight <- NULL
for(k in 1:nrow(trade2023_final2)){
  i <- which(countrycode_final_sort2$country_code==as.numeric(trade2023_final2[k,1]))
  j <- which(countrycode_final_sort2$country_code==as.numeric(trade2023_final2[k,2]))
  # e.weight[i,j] <- as.numeric(trade2023_final2[k,4])
  # e.weight[j,i] <- as.numeric(trade2023_final2[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(i,j,as.numeric(trade2023_final2[k,4])))
  e.sp.weight <- rbind(e.sp.weight, c(j,i,as.numeric(trade2023_final2[k,4])))
}

ITW2$sA <- e.sp.weight[nrow(e.sp.weight):1,]



G20_idx <- c(2,8,9,12,16,23,24,31,32,36,37,41,48,64,65,73,75,77)

countrycode_final_sort_G20 <- countrycode_final_sort2[G20_idx,]


ITW_G20 <- list()
ITW_G20$xy <- ITW2$xy[G20_idx,]

rownames(countrycode_final_sort_G20) <- 1:length(G20_idx)


ITW_G20_2023 <- ITW_G20
ITW_G20_2023$A <- matrix(0,0, nrow=N.ITW_G20, ncol=N.ITW_G20)
colnames(ITW_G20_2023$A) <- colnames(ITW_G20$A)
rownames(ITW_G20_2023$A) <- rownames(ITW_G20$A)
e.sp.weight <- NULL
for(k in 1:nrow(trade2023_final)){
  i <- which(countrycode_final_sort_G20$country_code==as.numeric(trade2023_final[k,1]))
  j <- which(countrycode_final_sort_G20$country_code==as.numeric(trade2023_final[k,2]))
  # e.weight[i,j] <- as.numeric(trade2023_final2[k,4])
  # e.weight[j,i] <- as.numeric(trade2023_final2[k,4])
  ITW_G20_2023$A[i,j] <- as.numeric(trade2023_final[k,4])
  ITW_G20_2023$A[j,i] <- as.numeric(trade2023_final[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(min(i,j),max(i,j),as.numeric(trade2023_final[k,4])))
}

ITW_G20_2023$sA <- e.sp.weight[nrow(e.sp.weight):1,]

L.ITW_G20_2023 <- laplacian_mat(ITW_G20_2023$A) # laplacian matrix
val1 <- eigensort(L.ITW_G20_2023)
evalues.ITW_G20_2023 <- val1$evalues
evectors.ITW_G20_2023 <- val1$evectors
# largest eigenvalue
lmax.ITW_G20_2023 <- max(evalues.ITW_G20_2023)

N.ITW_G20 <- nrow(L.ITW_G20_2023)


p.ITW2 <- length(files) - length(variable_remove)
R.ITW <- 34 # 1990~2023
X.ITW2 <- array(0, c(p.ITW2,N.ITW_G20,R.ITW))


for(i in 1:p.ITW){
  if(i %in% variable_remove){
    next
  }
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  print((as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% rownames(ITW$xy[-country_remove,]), 
                                              paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # tmp4 <- c(tmp4, (as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
  #                                                       paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # for(l in 1:N.ITW){
  #   X.ITW[i,,] <- as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3,
  #                               paste("X", 1990:2023, sep="")]) na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
  # }
}


# tmp4 <- c()
j <- 0
for(i in 1:p.ITW){
  if(i %in% variable_remove){
    next
  }
  j <- j+1
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  
  for(l in 1:N.ITW2){
    X.ITW2[j,l,] <- as.matrix(tmp[tmp$Country.Code %in% rownames(ITW$xy[-country_remove,]),
                                  paste("X", 1990:2023, sep="")])[l,] %>%  na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
  }
}



variables_selected <- c(6,7,8,11,2,4,15,5,14,17,22)

var_eng <- c("GDP growth", "GDP per capita growth", "GDP per capita", "GNI per capita",
             "Current account balance", "Exports of goods and services", "Imports of goods and services",
             "Foreign direct investment (net inflows)", "Gross capital formation",
             "Inflation", "Price level ratio")

var_kor <- c("국내총생산 성장률 (연간 %)", "1인당 GDP 성장률 (연간 %)", "1인당 GDP (현재 미국 달러 기준)", "1인당 총국민소득 (구매력 평가 기준, 국제 달러)",
             "경상수지 (국제수지 기준, 현재 미국 달러 기준)", "상품 및 서비스 수출액 (GDP 대비 %)", "상품 및 서비스 수입액 (GDP 대비 %)",
             "해외직접투자 순유입 (현재 미국 달러 기준)", "총자본형성 (GDP 대비 %)",
             "물가상승률 (GDP 디플레이터 기준)", "구매력 기준 환율과 시장 환율의 비율")

p.ITW_G20.std.selected <- length(variables_selected)

X.ITW_G20_2023 <- X.ITW2[variables_selected,G20_idx,34]
for (i in 1:nrow(X.ITW_G20_2023)) {
  v <- X.ITW_G20_2023[i,]
  X.ITW_G20_2023[i,] <- (v - mean(v)) / sd(v)
}

X.ITW_G20_2023.list <- list()
for(i in 1:p.ITW_G20.std.selected){
  X.ITW_G20_2023.list[[i]] <- X.ITW_G20_2023[i,]
}


res.ITW_G20.2023.random_window <- GFPCA(X=X.ITW_G20_2023.list, S=L.ITW_G20_2023, M=50, sigma=0.5, method="random")


# scree plot?
plot(colSums(res.ITW_G20.2023.random_window$tau.hat) / sum(res.ITW_G20.2023.random_window$tau.hat), type="o")

res.ITW_G20.2023.random_window <- GFPCA(X=X.ITW_G20_2023.list, S=L.ITW_G20_2023, q=1, M=50, sigma=0.5, method="random")


plot(res.ITW_G20.2023.random_window$tau.hat[,1], type="l", main="Graph spectral envelope",
     xlab="Graph frequency index", ylab=expression(tau[1]), cex.lab=2, cex.main=2, cex.axis=1.2)
abline(v=c(4,17), col="red", lty=2)

p1 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2023[,4], value="Value", ratio=0.6,
                         min=min(evectors.ITW_G20_2023[,4]), max=max(evectors.ITW_G20_2023[,4]), mg=c(4,4,4,4), title=expression(v[4]^"TN"), main.title.size = 20)

p2 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2023[,17], value="Value", ratio=0.6,
                         min=min(evectors.ITW_G20_2023[,17]), max=max(evectors.ITW_G20_2023[,17]), mg=c(4,4,4,4), title=expression(v[17]^"TN"), main.title.size = 20)

grid.arrange(p1,p2, nrow=1)

p3 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = res.ITW_G20.2023.random_window$Y[[1]], value="Value", ratio=0.6,
                         min=min(res.ITW_G20.2023.random_window$Y[[1]]), max=max(res.ITW_G20.2023.random_window$Y[[1]]), mg=c(4,4,4,4), title="PC 1", main.title.size = 20)

grid.arrange(p1,p2,p3, nrow=1)

par(mfrow=c(1,2))
plot(res.ITW_G20.2023.random_window$H.hat[4,1,], type="l", xlab="Dimension index", ylab="Normalized amplitude", main=expression(lambda[4]^"SM"), cex.main=2, cex.lab=2, cex.axis=1.3)
points(1:p.ITW_G20.std.selected, res.ITW_G20.2023.random_window$H.hat[4,1,], pch=1, cex=1.7)
abline(v=c(1,2,3,10,11), col="red", lty=2)
files[-variable_remove][variables_selected][c(1,2,3,10,11)]
var_eng[c(1,2,3,10,11)] ; var_kor[c(1,2,3,10,11)]

plot(res.ITW_G20.2023.random_window$H.hat[17,1,], type="l", xlab="Dimension index", ylab="Normalized amplitude", main=expression(lambda[17]^"SM"), cex.main=2, cex.lab=2, cex.axis=1.3)
points(1:p.ITW_G20.std.selected, res.ITW_G20.2023.random_window$H.hat[17,1,], pch=1, cex=1.7)
abline(v=c(5,8), col="red", lty=2)

files[-variable_remove][variables_selected][c(5,8)]
var_eng[c(5,8)] ; var_kor[c(5,8)]


## comparison
ITW.glpca.res <- c()
alpha.vec <- c()

# tmp.X.original <- tmp.X
# tmp.X <- tmp.X - rowMeans(tmp.X)
tmp.eigen <- eigen(t(X.ITW_G20_2023)%*%X.ITW_G20_2023)


for(beta.tmp in c(0, 0.0001, 0.001, 0.01, 0.1, 0.2,0.3,0.5,0.7,0.9,0.95)){
  lambda.n <- tmp.eigen$values[1]
  xi.n <- evalues.ITW_G20_2023[N.ITW_G20]
  alpha <- beta.tmp / (1-beta.tmp) * lambda.n / xi.n
  G.alpha <- -t(X.ITW_G20_2023)%*%X.ITW_G20_2023 + alpha*L.ITW_G20_2023
  eigen.G <- eigen(G.alpha)
  
  k <- 1
  Q.star <- as.matrix(eigen.G$vectors[,N.ITW_G20:(N.ITW_G20-k+1)])
  U.star <- as.matrix(X.ITW_G20_2023 %*% Q.star)
  
  ITW.glpca.res <- c(ITW.glpca.res, norm(X.ITW_G20_2023 - U.star%*%t(Q.star), type="F")^2 / norm(X.ITW_G20_2023, type="F")^2)
}

plot(c(0, 0.0001, 0.001, 0.01, 0.1, 0.2,0.3,0.5,0.7,0.9,0.95),
     ITW.glpca.res, type="l")

# > karate.glpca.res[1]
# [1] 0.461812

classic.eigen <- eigen(X.ITW_G20_2023%*%t(X.ITW_G20_2023))

norm(t(X.ITW_G20_2023)-t(X.ITW_G20_2023)%*%as.matrix(classic.eigen$vectors[,1:1])%*%t(classic.eigen$vectors[,1:1]), type="F")^2 / norm(X.ITW_G20_2023, type="F")^2
plot(classic.eigen$vectors[,1], type="l")

tmp.X.hat <- matrix(0, nrow=p.ITW_G20.std.selected, ncol=N.ITW_G20)

for(i in 1:p.ITW_G20.std.selected){
  tmp.X.hat[i,] <- res.ITW_G20.2023.random_window$X.hat[[i]]
}


1-norm(X.ITW_G20_2023 - tmp.X.hat, type="F") / norm(X.ITW_G20_2023, type="F") # 96.8%

norm(X.ITW_G20_2023, type="F")
