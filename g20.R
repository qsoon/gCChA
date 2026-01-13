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
library(expm)

source("utils.R")

# load trading network in 2023
load("tradingnw.RData")
plot_graph(ITW_G20_2023)


# make trading network in 2019
ITW_G20_2019 <- ITW_G20_2023

trade2019 <- read.csv("BACI_HS92_V202501/BACI_HS92_Y2019_V202501.csv", header=TRUE)
trade2019$q <- as.numeric(trade2019$q)
trade2019 <- na.omit(trade2019)

trade2019$value <- trade2019$v * trade2019$q


identical(sort(unique(trade2019$i)), sort(unique(trade2019$j)))
head(trade2019)

countrycode <- read.csv("BACI_HS92_V202501/country_codes_V202501.csv", header=TRUE)
countrycode_used <- countrycode[countrycode$country_iso3 %in% rownames(ITW_G20_2023$xy),]
countrycode_used <- countrycode_used[countrycode_used$country_code!=280,]

trade2019_used <- trade2019[(trade2019$i %in% countrycode_used$country_code) & (trade2019$j %in% countrycode_used$country_code),]


trade2019_final <- trade2019_used %>%
  # Create new variables i_new and j_new to treat (i, j) and (j, i) as the same pair
  mutate(i_new = pmin(i, j), j_new = pmax(i, j)) %>%
  # Group by these new variables
  group_by(i_new, j_new) %>%
  # Summarize by summing the total
  summarise(
    total = sum(value) / ifelse(sum(i != i_new | j != j_new) > 0, 2, 1),  # Divide by 2 if both directions exist
    .groups = 'drop'
  ) %>%
  # Rename columns back to original names if desired
  rename(i = i_new, j = j_new)


trade2019_final <- trade2019_final %>% arrange(desc(total))
# trade2019_final <- trade2019_final[1:round(0.03*nrow(trade2019_final)),]

trade2019_final$total <- trade2019_final$total/10^9
trade2019_final$total_log <- log(trade2019_final$total)

rm(trade2019, trade2019_used)



ITW_G20_2019$A <- matrix(0,0, nrow=N.ITW_G20, ncol=N.ITW_G20)
colnames(ITW_G20_2019$A) <- colnames(ITW_G20_2023$A)
rownames(ITW_G20_2019$A) <- rownames(ITW_G20_2023$A)
e.sp.weight <- NULL
for(k in 1:nrow(trade2019_final)){
  i <- which(countrycode_used$country_code==as.numeric(trade2019_final[k,1]))
  j <- which(countrycode_used$country_code==as.numeric(trade2019_final[k,2]))
  # e.weight[i,j] <- as.numeric(trade2022_final2[k,4])
  # e.weight[j,i] <- as.numeric(trade2022_final2[k,4])
  ITW_G20_2019$A[i,j] <- as.numeric(trade2019_final[k,4])
  ITW_G20_2019$A[j,i] <- as.numeric(trade2019_final[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(min(i,j),max(i,j), as.numeric(trade2019_final[k,4])))
}

ITW_G20_2019$sA <- e.sp.weight[nrow(e.sp.weight):1,]

plot_graph(ITW_G20_2019)


L.ITW_G20_2019 <- laplacian_mat(ITW_G20_2019$A) # laplacian matrix
val1 <- eigensort(L.ITW_G20_2019)
evalues.ITW_G20_2019 <- val1$evalues
evectors.ITW_G20_2019 <- val1$evectors
# largest eigenvalue
lmax.ITW_G20_2019 <- max(evalues.ITW_G20_2019)

N.ITW_G20 <- nrow(L.ITW_G20_2019)


# load economic data in 2019
# load("X.ITW_G20_2023.RData")
X.ITW_G20_2019 <- NULL
# p.ITW_G20 <- 11
p.ITW_G20 <- 5

files <- list.files(path = "./economic/Indicators", pattern = ".csv")
# [1] "balance.csv"            "export.csv"             "foreign_invest.csv"    
# [4] "GDP_growth.csv"         "GDP_per_cap_growth.csv" "GDP_per_cap.csv"       
# [7] "GNI_per_cap.csv"        "gross_capital.csv"      "import.csv"            
# [10] "inflation.csv"          "PLR.csv"   
file_idx <- c(2,3,4,9,10)
for(i in 1:p.ITW_G20){
  tmp <- read.csv(paste("economic/Indicators",files[file_idx[i]], sep="/"), skip=4)
  X.ITW_G20_2019 <- rbind(X.ITW_G20_2019, tmp[tmp$Country.Code %in% rownames(ITW_G20_2023$xy), "X2019"])
}

for (i in 1:nrow(X.ITW_G20_2019)) {
  v <- X.ITW_G20_2019[i,]
  X.ITW_G20_2019[i,] <- (v - mean(v)) / sd(v)
}


# load energy data in 2019
Y.ITW_G20_2019 <- NULL
q.ITW_G20 <- 5

files <- list.files(path = "./energy/Indicators", pattern = ".csv")
# [1] "elec_power_consump.csv" "energy_import.csv"      "energy_use.csv"        
# [4] "fossil_fuel.csv"        "fuel_export.csv"        "GDP_per_energy.csv"    
# [7] "renewable.csv"  
file_idx <- c(2,3,4,6,7)
# file_idx <- c(2,3,4,7)
for(i in 1:q.ITW_G20){
  tmp <- read.csv(paste("energy/Indicators",files[file_idx[i]], sep="/"), skip=4)
  Y.ITW_G20_2019 <- rbind(Y.ITW_G20_2019, tmp[tmp$Country.Code %in% rownames(ITW_G20_2023$xy), "X2019"])
}

for (i in 1:nrow(Y.ITW_G20_2019)) {
  v <- Y.ITW_G20_2019[i,]
  Y.ITW_G20_2019[i,] <- (v - mean(v)) / sd(v)
}


X.ITW_G20_2019.list <- list()
for(i in 1:p.ITW_G20){
  X.ITW_G20_2019.list[[i]] <- X.ITW_G20_2019[i,]
}

Y.ITW_G20_2019.list <- list()
for(i in 1:q.ITW_G20){
  Y.ITW_G20_2019.list[[i]] <- Y.ITW_G20_2019[i,]
}


gcca_res_window <- gCChA1(X=X.ITW_G20_2019.list, Y = Y.ITW_G20_2019.list, evalues = evalues.ITW_G20_2019, evectors = evectors.ITW_G20_2019,
                        S=L.ITW_G20_2019, M=50, sigma=0.5, method="random", cal="sym")


r.ITW_G20 <- min(p.ITW_G20, q.ITW_G20)
Z.ITW_G20_2019 <- matrix(0, nrow=r.ITW_G20, ncol=N.ITW_G20)
for(i in 1:r.ITW_G20){
  Z.ITW_G20_2019[i,] <- as.vector(gcca_res_window$Z[[i]])
}

W.ITW_G20_2019 <- matrix(0, nrow=r.ITW_G20, ncol=N.ITW_G20)
for(i in 1:r.ITW_G20){
  W.ITW_G20_2019[i,] <- as.vector(gcca_res_window$W[[i]])
}


# canonical graph loadings
loadings.mat_Z <- array(NA, dim=c(N.ITW_G20, r.ITW_G20, p.ITW_G20))
loadings.mat_Z_sign <- array(NA, dim=c(N.ITW_G20, r.ITW_G20, p.ITW_G20))
# loadings.mat_Z2 <- array(NA, dim=c(N.ITW_G20, r.ITW_G20, p.ITW_G20))

for(l in 1:N.ITW_G20){
  for(j in 1:r.ITW_G20){
    for(k in 1:p.ITW_G20){
      # seed <- i*j+k+j*k
      # loadings.mat_Z[,j,k] <- coherence.graph(x=Z.ITW_G20_2019[j,], y=X.ITW_G20_2019.list[[k]], S=L.ITW_G20_2019, M=50, sigma=0.5, method="random", seed=seed)
      loadings.mat_Z[l,j,k] <- as.numeric((t(gcca_res_window$H.hat[l,j,]) %*% gcca_res_window$Px[l,,k])^2 / 
                                            (t(gcca_res_window$H.hat[l,j,]) %*% gcca_res_window$Px[l,,] %*% gcca_res_window$H.hat[l,j,]) /
                                            gcca_res_window$Px[l,k,k])
      loadings.mat_Z_sign[l,j,k] <- as.numeric(as.numeric(t(gcca_res_window$H.hat[l,j,]) %*% gcca_res_window$Px[l,,k]) >= 0)
    }
  }
}




loadings.mat_W <- array(NA, dim=c(N.ITW_G20, r.ITW_G20, q.ITW_G20))
loadings.mat_W_sign <- array(NA, dim=c(N.ITW_G20, r.ITW_G20, q.ITW_G20))

for(l in 1:N.ITW_G20){
  for(j in 1:r.ITW_G20){
    for(k in 1:q.ITW_G20){
      seed <- i+j*k+i*k
      # loadings.mat_W[,j,k] <- coherence.graph(x=W.ITW_G20_2019[j,], y=Y.ITW_G20_2019.list[[k]], S=L.ITW_G20_2019, M=50, sigma=0.5, method="random", seed=seed)
      loadings.mat_W[l,j,k] <- as.numeric((t(gcca_res_window$F.hat[l,j,]) %*% gcca_res_window$Py[l,,k])^2 / 
                                            (t(gcca_res_window$F.hat[l,j,]) %*% gcca_res_window$Py[l,,] %*% gcca_res_window$F.hat[l,j,]) /
                                            gcca_res_window$Py[l,k,k])
      loadings.mat_W_sign[l,j,k] <- as.numeric(as.numeric(t(gcca_res_window$F.hat[l,j,]) %*% gcca_res_window$Py[l,,k]) >= 0)
    }
  }
}

canonical_coherence_sign <- array(NA, dim=c(N.ITW_G20, r.ITW_G20))
for(l in 1:N.ITW_G20){
  for(j in 1:r.ITW_G20){
      # seed <- i*j+k+j*k
      # loadings.mat_Z[,j,k] <- coherence.graph(x=Z.ITW_G20_2019[j,], y=X.ITW_G20_2019.list[[k]], S=L.ITW_G20_2019, M=50, sigma=0.5, method="random", seed=seed)
    canonical_coherence_sign[l,j] <- as.numeric(as.numeric(t(gcca_res_window$H.hat[l,j,]) %*% gcca_res_window$Pxy[l,,] %*% gcca_res_window$F.hat[l,j,])>=0)
    }
  }



############################################################################
############################################################################
par(mfrow=c(2,4), mar=c(5,5,4,2)+0.1)
par(xpd=FALSE)
tol=0.2
plot(loadings.mat_Z[1,1,], type="b", ylim=c(0,1), main=expression(hat(Z)[1] * ", " * lambda[1]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("X[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_Z[1,1,], col=ifelse(loadings.mat_Z_sign[1,1,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_W[1,1,], type="b", ylim=c(0,1), main=expression(hat(W)[1] * ", " * lambda[1]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("Y[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_W[1,1,], col=ifelse(loadings.mat_W_sign[1,1,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_Z[1,2,], type="b", ylim=c(0,1), main=expression(hat(Z)[2] * ", " * lambda[1]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("X[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_Z[1,2,], col=ifelse(loadings.mat_Z_sign[1,2,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_W[1,2,], type="b", ylim=c(0,1), main=expression(hat(W)[2] * ", " * lambda[1]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("Y[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_W[1,2,], col=ifelse(loadings.mat_W_sign[1,2,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_Z[7,1,], type="b", ylim=c(0,1), main=expression(hat(Z)[1] * ", " * lambda[7]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("X[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_Z[7,1,], col=ifelse(loadings.mat_Z_sign[7,1,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_W[7,1,], type="b", ylim=c(0,1), main=expression(hat(W)[1] * ", " * lambda[7]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("Y[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_W[7,1,], col=ifelse(loadings.mat_W_sign[7,1,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_Z[7,2,], type="b", ylim=c(0,1), main=expression(hat(Z)[2] * ", " * lambda[7]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("X[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_Z[7,2,], col=ifelse(loadings.mat_Z_sign[7,2,] == 1, "red", "blue"), pch=19, cex=2)

plot(loadings.mat_W[7,2,], type="b", ylim=c(0,1), main=expression(hat(W)[2] * ", " * lambda[7]), ylab="Canonical graph loadings",xaxt='n', xlab="", cex.lab=1.4, cex.main=1.6)
axis(1, at=1:5, labels=parse(text=paste0("Y[", 1:5, "]")),cex.axis=1.4)
abline(h=tol, lty=2)
points(loadings.mat_W[7,2,], col=ifelse(loadings.mat_W_sign[7,2,] == 1, "red", "blue"), pch=19, cex=2)


p1 <- plot_graph_custom4(ITW_G20_2019, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2019[,1], value="Value", ratio=0.6,
                         min=min(evectors.ITW_G20_2019[,7]), max=max(evectors.ITW_G20_2019[,7]), mg=c(4,4,4,4), title=expression(v[1]^"TN"), main.title.size = 20)


p2 <- plot_graph_custom4(ITW_G20_2019, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2019[,7], value="Value", ratio=0.6,
                   min=min(evectors.ITW_G20_2019[,7]), max=max(evectors.ITW_G20_2019[,7]), mg=c(4,4,4,4), title=expression(v[7]^"TN"), main.title.size = 20)


grid.arrange(p1,p2, nrow=1)


plot_graph_custom3(ITW_G20_2019, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2019[,4], value="Value", ratio=0.6,
                   min=min(ITW_G20_2019$sA[,3], ITW_G20_2019$sA[,3]), max=max(ITW_G20_2019$sA[,3], ITW_G20_2019$sA[,3]), mg=c(4,4,4,4), title="", main.title.size = 20, signal=FALSE)


# graph canonical coherence
par(mfrow=c(1,1), mar=c(5,5,4,8)+0.1)
par(xpd=TRUE)  
plot(gcca_res_window$gamma.hat[,5], type="l", xlab = expression(Graph~frequency~index~"\u2113"),
     ylab = expression(Graph~canonical~coherence~~hat(gamma)[i](lambda["\u2113"])), 
     cex.lab=1.2, col="cyan", ylim=c(0,1))
points(gcca_res_window$gamma.hat[,5], pch=25, cex=1.5, col="cyan", bg="cyan")

lines(gcca_res_window$gamma.hat[,4], col="magenta")
points(gcca_res_window$gamma.hat[,4], pch=23, cex=1.5, col="magenta", bg="magenta")

lines(gcca_res_window$gamma.hat[,3], col="blue")
points(gcca_res_window$gamma.hat[,3], pch=22, cex=1.5, col="blue", bg="blue")

lines(gcca_res_window$gamma.hat[,2], col="red")
points(gcca_res_window$gamma.hat[,2], pch=24, cex=1.5, col="red", bg="red")

lines(gcca_res_window$gamma.hat[,1], col="black")
points(gcca_res_window$gamma.hat[,1], pch=19, cex=1.5, col="black")

legend("topright",
       inset=c(-0.35,0),  # inset을 음수로 주면 plot 밖으로 밀림
       legend=c(expression((hat(Z)[1]*", "*hat(W)[1])), expression((hat(Z)[2]*", "*hat(W)[2])),
                expression((hat(Z)[3]*", "*hat(W)[3])), expression((hat(Z)[4]*", "*hat(W)[4])),
                expression((hat(Z)[5]*", "*hat(W)[5]))),
       col=c("black","red","blue","magenta","cyan"), 
       pt.bg=c("black","red","blue","magenta","cyan"), lty=1, 
       pch=c(19,24,22,23,25), bty="n", cex=1.1)

plot(gcca_res_window$gamma.hat[,1], type="l")
lines(gcca_res_window$gamma.hat[,2], col="red")

