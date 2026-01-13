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



#######################################
# USPS dataset
#######################################
install.packages("IMIFA")
library(IMIFA)
data(USPSdigits)

USPSdigits$train[,-1]
USPSdigits$train[,1]
data(usps)

table(usps$label)
dim(usps$data) # 11000 x 256 (each image = 16 x 16)



# visualize
opar <- par(no.readonly=TRUE, mfrow=c(1,3), pty="s")
image(t(matrix(usps$data[digit_indices[7],],nrow=16)[16:1,])) # last of digit 4
image(t(matrix(usps$data[7502,],nrow=16)[16:1,])) # last of digit 9
image(t(matrix(usps$data[6600,],nrow=16)[16:1,])) # last of digit 6
par(opar)

digit_indices <- sapply(0:9, function(d){
  head(which(USPSdigits$train[,1] == d), 1)
})

# plot 10 digits in 2-by-5 grid
op <- par(no.readonly=TRUE)
# par(mfrow=c(2,5), pty="s", mar=c(1,1,1,1))
par(mfrow=c(2,5), pty="s", mar=c(0.3,0.3,0.8,0.3), oma=c(0,0,0,0))
# par(mfrow=c(2,5), pty="s", mai=c(0.05, 0.05, 0.05, 0.05))
for (i in 1:10) {
  idx <- digit_indices[i]
  img <- matrix(as.numeric(USPSdigits$train[idx,-1]), nrow=16)
  img <- img[,16:1]
  # img <- t(img[1:16,])
  # img <- t(img[nrow(img):1, ])# rotate to correct USPS orientation
  # Convert to 0–1 scale for plotting
  mat <- (img+1) / 2   # -1→0, 1→1
  
  # Plot (black background, white foreground)
  image(mat, axes = FALSE, col = gray(seq(0,1,len=256)))
}




digit_indices <- sapply(0:9, function(d){
  tail(which(usps$label == d), 1)
})

# plot 10 digits in 2-by-5 grid
op <- par(no.readonly=TRUE)
par(mfrow=c(2,5), pty="s", mar=c(1,1,1,1))
par(mfrow=c(2,5), pty="s", mar=c(0.3,0.3,0.8,0.3), oma=c(0,0,0,0))

for (i in 1:10) {
  idx <- digit_indices[i]
  img <- matrix(usps$data[idx,], nrow=16)
  img <- t(img[16:1,])   # rotate to correct USPS orientation
  image(img, col=gray(0:256/256), axes=FALSE)
  # title(paste("Digit", i-1), cex.main=1.0)
}



avg_accuracy_list_gmcca <- NULL
avg_accuracy_list_gccha <- NULL



for(prop in c(1/4,3/8,1/2,5/8,3/4)){
  cat("############ Proportion:", prop, "############", end="\n")
  accuracy_list_gmcca <- c()
  accuracy_list_gccha <- c()
  reduced_dim <- 10
  for(iter_ind in 1:50){
    set.seed(100*iter_ind)
    cat("### iter:", iter_ind, "###", end="\n")
    # 각 클래스별 샘플 수 계산
    n_classes <- 10
    n_train_per_class <- 400 / n_classes  # 40
    
    train_idx <- unlist(lapply(0:(n_classes-1), function(lbl){
      idx <- which(USPSdigits$train[,1] == lbl)
      sample(idx, n_train_per_class)
    }))
    
    # train
    train_small_images <- as.matrix(USPSdigits$train[train_idx, -1])
    train_small_labels <- USPSdigits$train[train_idx,1]
    
    n_train.usps <- nrow(train_small_images)
    p.usps <- 256*prop # first (16*prop) rows
    q.usps <- ncol(train_small_images) - p.usps
    
    X_train <- train_small_images[,1:p.usps]
    Y_train <- train_small_images[,(p.usps+1):ncol(train_small_images)]
    
    # 결과 edge 리스트 저장용
    edge_list <- list()
    weight_list <- numeric()
    
    # 클래스별로 처리
    unique_classes <- sort(unique(train_small_labels))
    
    for (cls in unique_classes) {
      
      cat("Processing class", cls, "\n")
      
      # 해당 클래스 인덱스
      idx <- which(train_small_labels == cls)
      
      # 클래스 내 데이터 (n_train.usps × 256)
      Xc <- train_small_images[idx, ]
      
      # cosine similarity matrix (n_train.usps × n_train.usps)
      # crossprod for speed: cos(A,B) = A·B / (||A|| ||B||)
      
      # 정규화
      Xc_norm <- Xc / sqrt(rowSums(Xc^2))
      
      # similarity matrix
      S <- Xc_norm %*% t(Xc_norm)
      
      # self similarity 제거
      diag(S) <- 0
      
      # 모든 쌍을 edge로 추가 (39-NN)
      n_class <- length(idx)
      for (i in 1:n_class) {
        # i → 모든 j
        js <- setdiff(1:n_class, i)
        
        # global index 변환
        from <- rep(idx[i], length(js))
        to <- idx[js]
        w <- S[i, js]
        
        edge_list[[length(edge_list) + 1]] <- cbind(from, to)
        weight_list <- c(weight_list, w)
      }
    }
    
    # edge list 하나로 합치기
    edges <- do.call(rbind, edge_list)
    
    # igraph 생성
    g <- graph_from_data_frame(
      data.frame(
        from = edges[,1],
        to = edges[,2],
        weight = weight_list
      ),
      directed = FALSE
    )
    
    # 1. 가중치 adjacency matrix 만들기
    A <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
    
    # 2. Degree matrix
    D <- diag(rowSums(A))
    
    # 3. Unnormalized graph Laplacian
    L_usps <- D - A
    
    X_total_train <- list(t(X_train), t(Y_train))
    # gmcca_res <- gmcca(X_total_train, L_usps, gamm=0.1, M=2, d=30, eps = 1e-4) 
    # 
    # dim(gmcca_res$U[[1]]) ; dim(gmcca_res$U[[2]]) # 64 x 30 / 192 x 30
    
    
    gamma_candidate <- c(0.001,0.01,0.1,1,10,100,1000)
    accuracy <- NULL
    for(gamma in gamma_candidate){
      print(paste("gamma = ",gamma))
      gmcca_res <- gmcca(X_total_train, L_usps, gamm=gamma, M=2, d=reduced_dim, eps = 1e-4) 
      
      reduced_train <- rbind(t(gmcca_res$U[[1]]) %*% t(X_train),
                             t(gmcca_res$U[[2]]) %*% t(Y_train))
      
      k <- 10
      
      # 10-NN search: 각 test sample마다 가장 가까운 k train samples index
      # nn_res <- get.knnx(t(reduced_train), t(reduced_train), k = k)
      nn_res_raw <- get.knnx(t(reduced_train), t(reduced_train), k = k+1)
      
      predicted_labels <- sapply(1:ncol(reduced_train), function(i){
        # idx <- nn_res$nn.index[i, ]   # test sample i의 10 nearest train indices
        # dist <- nn_res$nn.dist[i, ]   # 해당 거리
        # neighbor_labels <- train_small_labels[idx]
        
        idx_raw <- nn_res_raw$nn.index[i, ]
        dist_raw <- nn_res_raw$nn.dist[i, ]
        
        # 자기 자신 제거: index가 i인 경우 필터링
        keep <- idx_raw != i
        idx <- idx_raw[keep][1:k]
        dist <- dist_raw[keep][1:k]
        neighbor_labels <- train_small_labels[idx]
        
        # 거리 순으로 정렬 (가장 가까운 것 먼저)
        ord <- order(dist)
        neighbor_labels <- neighbor_labels[ord]
        
        # 테이블 생성
        tab <- table(neighbor_labels)
        max_count <- max(tab)
        max_labels <- names(tab)[tab == max_count]
        
        # 동률이면 가장 가까운 neighbor label로 결정
        predicted_label <- if(length(max_labels) == 1){
          max_labels
        } else {
          # max_labels 중 첫 번째로 가장 가까운 neighbor label 선택
          for(lbl in neighbor_labels){
            if(lbl %in% max_labels){
              predicted_label <- lbl
              break
            }
          }
          predicted_label
        }
        as.numeric(predicted_label)
      })
      
      # 정확도 계산 (실제 test label이 있을 경우)
      accuracy <- c(accuracy, mean(predicted_labels == train_small_labels))
    }
    
    accuracy_list_gmcca[iter_ind] <- max(accuracy)
    
    
    # gCChA
    val_usps <- eigensort(L_usps)
    evalues_usps <- val_usps$evalues
    evectors_usps <- val_usps$evectors
    
    gcca_res_window_usps <- gCChA2(X=X_train, Y = Y_train, 
                                        S=L_usps, evalues=evalues_usps, evectors=evectors_usps,
                                        M=100, sigma=0.05, method="random", cal="asym")
    
    
    r <- reduced_dim
    reduced_train_gcca_window <- rbind(t(gcca_res_window_usps$Z[,1:r]),
                                       t(gcca_res_window_usps$W[,1:r]))
    
    k <- 10
    
    # 10-NN search: 각 test sample마다 가장 가까운 k train samples index
    # nn_res <- get.knnx(t(reduced_train), t(reduced_train), k = k)
    nn_res_raw_gcca_window <- get.knnx(t(reduced_train_gcca_window), t(reduced_train_gcca_window), k = k+1)
    
    predicted_labels_gcca_window_train <- sapply(1:ncol(reduced_train_gcca_window), function(i){
      # idx <- nn_res$nn.index[i, ]   # test sample i의 10 nearest train indices
      # dist <- nn_res$nn.dist[i, ]   # 해당 거리
      # neighbor_labels <- train_small_labels[idx]
      
      idx_raw <- nn_res_raw_gcca_window$nn.index[i, ]
      dist_raw <- nn_res_raw_gcca_window$nn.dist[i, ]
      
      # 자기 자신 제거: index가 i인 경우 필터링
      keep <- idx_raw != i
      idx <- idx_raw[keep][1:k]
      dist <- dist_raw[keep][1:k]
      neighbor_labels <- train_small_labels[idx]
      
      # 거리 순으로 정렬 (가장 가까운 것 먼저)
      ord <- order(dist)
      neighbor_labels <- neighbor_labels[ord]
      
      # 테이블 생성
      tab <- table(neighbor_labels)
      max_count <- max(tab)
      max_labels <- names(tab)[tab == max_count]
      
      # 동률이면 가장 가까운 neighbor label로 결정
      predicted_label <- if(length(max_labels) == 1){
        max_labels
      } else {
        # max_labels 중 첫 번째로 가장 가까운 neighbor label 선택
        for(lbl in neighbor_labels){
          if(lbl %in% max_labels){
            predicted_label <- lbl
            break
          }
        }
        predicted_label
      }
      as.numeric(predicted_label)
    })
    
    # 정확도 계산 (실제 test label이 있을 경우)
    accuracy_gcca_window_train <- mean(predicted_labels_gcca_window_train == train_small_labels)
    
    accuracy_list_gccha[iter_ind] <- accuracy_gcca_window_train
  }
  
  avg_accuracy_list_gmcca <- rbind(avg_accuracy_list_gmcca, c(prop, mean(accuracy_list_gmcca), sd(accuracy_list_gmcca)))
  avg_accuracy_list_gccha <- rbind(avg_accuracy_list_gccha, c(prop, mean(accuracy_list_gccha), sd(accuracy_list_gccha)))
}


cbind(avg_accuracy_list_gmcca, avg_accuracy_list_gccha)
