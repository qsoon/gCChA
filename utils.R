# localization operator
localization.op <- function(evalues, evectors, g, i=NULL){ # support of T_i g is centered at node i
  if(is.null(i)){
    res <- evectors %*% (t(Conj(evectors))*g(evalues)) # ith column: T_i g
  } else{
    res <- as.vector(evectors %*% (Conj(evectors)[i,]*g(evalues))) # T_i g vector when i is specified
  }
  return(res)
}

# graph windowed Fourier transform
graph.window.FT <- function(x, S, g, M){
  # eigenres <- eigen(S)
  # evalues <- eigenres$values
  # evectors <- eigenres$vectors
  # lmax <- max(evalues)
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  lmax <- max(evalues)
  
  C <- NULL
  normconst <- c()
  for(m in 1:M){
    gm <- function(lambda){
      return(g(lambda, sigma.sq=lmax*(M+1)/M^2, m, tau=lmax*(M+1)/M^2))
    }
    Tgm <- localization.op(evalues, evectors, gm)
    C <- cbind(C, t(Conj(Tgm))%*%x) # where C[i,m] = <x, Ti gm>
    normconst <- c(normconst, norm(Tgm, type="F")^2)
  }
  return(list(C=C, normconst=normconst))
}

# cpsd estimate
cpsd.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # x, y should be R realizations (N x R matrix)
    # eigenres <- eigen(S)
    # x.tilde <- t(Conj(eigenres$vectors))%*%x
    # y.tilde <- t(Conj(eigenres$vectors))%*%y
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    x.tilde <- t(Conj(evectors))%*%x
    y.tilde <- t(Conj(evectors))%*%y
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  } else if(method=="window"){
    # x, y should be one realization (N x 1 vector)
    C1 <- graph.window.FT(x, S, g, M)
    C2 <- graph.window.FT(y, S, g, M)
    cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
  } else if(method=="random"){
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    x.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(x))
    y.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(y))
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  }
  return(cpsd)
}

cpsd.graph.fast <- function(x, y, evalues, evectors, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # x, y should be R realizations (N x R matrix)
    # eigenres <- eigen(S)
    # x.tilde <- t(Conj(eigenres$vectors))%*%x
    # y.tilde <- t(Conj(eigenres$vectors))%*%y
    x.tilde <- t(Conj(evectors))%*%x
    y.tilde <- t(Conj(evectors))%*%y
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  } else if(method=="window"){
    # x, y should be one realization (N x 1 vector)
    C1 <- graph.window.FT(x, S, g, M)
    C2 <- graph.window.FT(y, S, g, M)
    cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
  } else if(method=="random"){
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    x.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(x))
    y.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(y))
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  }
  return(cpsd)
}


# psd estimate
psd.graph <- function(x, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  return(cpsd.graph(x=x, y=x, S=S, g=g, M=M, sigma=sigma, method=method, seed=seed))
}

# coherence estimate
coherence.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, seed=seed)
  } else if(method=="window"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, method=method, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, method=method, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, method=method, seed=seed)
  } else if(method=="random"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, M=M, sigma=sigma, method=method, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, M=M, sigma=sigma, method=method, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, M=M, sigma=sigma, method=method, seed=seed)
  }
  return(cpsd*Conj(cpsd) / psd.x / psd.y)
}

# cross spectrum analysis All in one
cross.spectrum.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  return(list(cpsd=cpsd, psd.x = psd.x, psd.y = psd.y, coherence=cpsd*Conj(cpsd) / psd.x / psd.y))
}


# random window generation
windowbank.random <- function(N, M, V, sigma, seed=1){
  res <- matrix(0, nrow=M, ncol=N)
  set.seed(seed)
  for(i in 1:M){
    W.tilde <- diag(N) + rnorm(N^2, 0, sigma)
    res[i,] <- diag(V %*% W.tilde %*% t(Conj(V)))
  }
  return(res)
}


plot_graph_custom <- function (z, size = 0.75, edge_color, vertex_color=NULL) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$color <- factor(edge_color, levels = unique(edge_color))
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, colour = color), 
                                              linewidth = 2, data = df2) +
    scale_color_manual(values=color.cand, labels = paste("Line", 1:8, sep=""),
                       name = "Line number") + 
    geom_point(aes(fill=vertex_color), size = size, shape=21) + 
    scale_fill_gradient(low="white", high="black", na.value = "yellow", name = "People") +
    theme_void() +
    theme(legend.margin = margin(10,10,10,10))
  print(p1)
}

plot_graph_custom3 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  # if(signal==FALSE){
  #   p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
  #                                                   xend = y1, yend = y2),  linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
  #     geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
  #   print(p1)
  # }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w),  linewidth=e.size, data = df2) + ggtitle(title) +
      scale_color_gradient(low="lightgray", high="black", name="Edge Weight", guide="none") + 
      geom_point(size = v.size, fill="white", colour="black", shape=21, stroke=1.2) + theme_void()+
      geom_label_repel(aes(label = w), size = 3, max.overlaps = Inf, box.padding = 0.3, segment.color = "gray50")+
      theme(plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5),aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}


plot_graph_custom4 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom5 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2, show.legend=FALSE) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(color=guide_colourbar(order=1)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}


gCChA1 <- function(X, Y, S, evalues, evectors, q=NULL, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1, cal="asym"){
  if(is.null(method)){
    # X[[i]] should be R realizations (n x R matrix)
    # p-dimensional graph signal
    p <- length(X)
    q <- length(Y)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    if(nrow(Y[[1]])!=nrow(S)){
      stop("The number of rows of list Y's element matrices should be equal to the number of rows of GSO!")
    }
    if(ncol(X[[1]])!=ncol(Y[[1]])){
      stop("The number of columns of list X's element matrices should be equal to the number of columns of list Y's element matrices!")
    }
    
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    P_array_X <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        P_array_X[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S) 
      }
    }
    
    P_array_Y <- array(NA, dim = c(n,q,q))
    for(i in 1:q){
      for(j in 1:q){
        P_array_Y[,i,j] <- cpsd.graph(x=Y[[i]], y=Y[[j]], S=S) 
      }
    }
    
    P_array_XY <- array(NA, dim = c(n,p,q))
    for(i in 1:p){
      for(j in 1:q){
        P_array_XY[,i,j] <- cpsd.graph(x=X[[i]], y=Y[[j]], S=S) 
      }
    }
    
    r <- min(p,q)
    H.hat <- array(NA, dim = c(n,r,p))
    F.hat <- array(NA, dim = c(n,r,q))
    eta.hat <- array(NA, dim = c(n,q,r))
    gamma.hat <- matrix(NA, nrow=n, ncol=r)
    
    P_array_YX <- array(NA, dim = c(n,q,p))
    for(l in 1:n){
      P_array_YX[l,,] <- Conj(t(P_array_XY[l,,]))
      if(cal=="sym"){
        targetmat <- solve(expm::sqrtm(P_array_Y[l,,])) %*% P_array_YX[l,,] %*% solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,]))
        tmp <- eigensort(targetmat)
        eta.hat[l,,] <- tmp$evectors[,q:(q-r+1)]
        gamma.hat[l,] <- rev(tmp$evalues)[1:r]
        
        # H.hat[l,,] <- t(solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        # H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        H.hat[l,,] <- t(solve(expm::sqrtm(P_array_X[l,,])) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        H.hat[l,,] <- t(solve(expm::sqrtm(P_array_X[l,,])) %*% t(H.hat[l,,]))
        
        F.hat[l,,] <- t(eta.hat[l,,])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
        F.hat[l,,] <- t(solve(expm::sqrtm(P_array_Y[l,,])) %*% t(F.hat[l,,]))
        
        # F.hat[l,,] <- t(solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        # F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
      } else if(cal=="asym"){
        targetmat1 <- solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(P_array_Y[l,,]) %*% P_array_YX[l,,]
        targetmat2 <- solve(P_array_Y[l,,]) %*% P_array_YX[l,,] %*% solve(P_array_X[l,,]) %*% P_array_XY[l,,]
        
        tmp1 <- eigensort(targetmat1)
        gamma.hat[l,] <- rev(tmp1$evalues)[1:r]
        H.hat[l,,] <- t(tmp1$evectors[,q:(q-r+1)])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        
        tmp2 <- eigensort(targetmat2)
        F.hat[l,,] <- t(tmp2$evectors[,q:(q-r+1)])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
      }
    }

    
    Z <- list()
    W <- list()
    # Z[[i]] : canonical graph signal corresponding to X for R realizations (n x R matrix) / length r
    # W[[i]] : canonical graph signal corresponding to Y for R realizations (n x R matrix) / length r
    for(i in 1:r){
      Z[[i]] <- matrix(0, nrow=n, ncol=R)
      W[[i]] <- matrix(0, nrow=n, ncol=R)
      
      for(j in 1:p){
        Z[[i]] <- Z[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
      }
      for(j in 1:q){
        W[[i]] <- W[[i]] + evectors %*% diag(F.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
      }
    }
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$Z <- Z ; res$W <- W
    res$eta.hat <- eta.hat ; res$gamma.hat <- gamma.hat ; res$H.hat <- H.hat ; res$F.hat <- F.hat
    res$Px <- P_array_X ; res$Py <- P_array_Y ; res$Pxy <- P_array_XY ; res$Pyx <- P_array_YX
    return(res)
  } 
  
  else if(method=="random"){
    # X[[i]] should be one realization (n x 1 vector)
    # p-dimensional graph signal
    if(is.vector(X[[1]])){
      X <- lapply(X, as.matrix)
    }
    p <- length(X)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The length of list X's element vector OR the number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    
    if(is.vector(Y[[1]])){
      Y <- lapply(Y, as.matrix)
    }
    q <- length(Y)
    if(nrow(Y[[1]])!=nrow(S)){
      stop("The number of rows of list Y's element matrices should be equal to the number of rows of GSO!")
    }
    
    if(ncol(X[[1]])!=ncol(Y[[1]])){
      stop("The number of columns of list X's element matrices should be equal to the number of columns of list Y's element matrices!")
    }
    
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    # val <- eigensort(S)
    # evalues <- val$evalues
    # evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    print("calculating P_array_X...")
    P_array_X <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        # P_array_X[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S, g=g, M=M, sigma=sigma, method=method, seed=seed) 
        P_array_X[,i,j] <- cpsd.graph.fast(x=X[[i]], y=X[[j]], evalues=evalues, evectors=evectors, g=g, M=M, sigma=sigma, method=method, seed=seed) 
      }
    }
    
    print("calculating P_array_Y...")
    P_array_Y <- array(NA, dim = c(n,q,q))
    for(i in 1:q){
      for(j in 1:q){
        # P_array_Y[,i,j] <- cpsd.graph(x=Y[[i]], y=Y[[j]], S=S, g=g, M=M, sigma=sigma, method=method, seed=seed) 
        P_array_Y[,i,j] <- cpsd.graph.fast(x=Y[[i]], y=Y[[j]], evalues=evalues, evectors=evectors, g=g, M=M, sigma=sigma, method=method, seed=seed) 
      }
    }
    
    print("calculating P_array_XY...")
    P_array_XY <- array(NA, dim = c(n,p,q))
    for(i in 1:p){
      for(j in 1:q){
        # P_array_XY[,i,j] <- cpsd.graph(x=X[[i]], y=Y[[j]], S=S, g=g, M=M, sigma=sigma, method=method, seed=seed) 
        P_array_XY[,i,j] <- cpsd.graph.fast(x=X[[i]], y=Y[[j]], evalues=evalues, evectors=evectors, g=g, M=M, sigma=sigma, method=method, seed=seed) 
      }
    }
    
    r <- min(p,q)
    H.hat <- array(NA, dim = c(n,r,p))
    F.hat <- array(NA, dim = c(n,r,q))
    eta.hat <- array(NA, dim = c(n,q,r))
    gamma.hat <- matrix(NA, nrow=n, ncol=r)
    
    P_array_YX <- array(NA, dim = c(n,q,p))
    for(l in 1:n){
      print(paste("iteration:", l))
      P_array_YX[l,,] <- Conj(t(P_array_XY[l,,]))
      if(cal=="sym"){
        targetmat <- solve(expm::sqrtm(P_array_Y[l,,])) %*% P_array_YX[l,,] %*% solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,]))
        tmp <- eigensort(targetmat)
        eta.hat[l,,] <- tmp$evectors[,q:(q-r+1)]
        gamma.hat[l,] <- rev(tmp$evalues)[1:r]
        
        H.hat[l,,] <- t(solve(expm::sqrtm(P_array_X[l,,])) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        H.hat[l,,] <- t(solve(expm::sqrtm(P_array_X[l,,])) %*% t(H.hat[l,,]))
        
        F.hat[l,,] <- t(eta.hat[l,,])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
        F.hat[l,,] <- t(solve(expm::sqrtm(P_array_Y[l,,])) %*% t(F.hat[l,,]))
      } else if(cal=="asym"){
        targetmat1 <- solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(P_array_Y[l,,]) %*% P_array_YX[l,,]
        targetmat2 <- solve(P_array_Y[l,,]) %*% P_array_YX[l,,] %*% solve(P_array_X[l,,]) %*% P_array_XY[l,,]
        
        tmp1 <- eigensort(targetmat1)
        gamma.hat[l,] <- rev(tmp1$evalues)[1:r]
        H.hat[l,,] <- t(tmp1$evectors[,q:(q-r+1)])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        
        tmp2 <- eigensort(targetmat2)
        F.hat[l,,] <- t(tmp2$evectors[,q:(q-r+1)])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
      }
    }
    
    
    Z <- list()
    W <- list()
    # Z[[i]] : canonical graph signal corresponding to X for R realizations (n x R matrix) / length r
    # W[[i]] : canonical graph signal corresponding to Y for R realizations (n x R matrix) / length r
    for(i in 1:r){
      Z[[i]] <- matrix(0, nrow=n, ncol=R)
      W[[i]] <- matrix(0, nrow=n, ncol=R)
      
      for(j in 1:p){
        Z[[i]] <- Z[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
      }
      for(j in 1:q){
        W[[i]] <- W[[i]] + evectors %*% diag(F.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
      }
    }
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$Z <- Z ; res$W <- W
    res$eta.hat <- eta.hat ; res$gamma.hat <- gamma.hat ; res$H.hat <- H.hat ; res$F.hat <- F.hat
    res$Px <- P_array_X ; res$Py <- P_array_Y ; res$Pxy <- P_array_XY ; res$Pyx <- P_array_YX

    return(res)
  }
}


gCChA2 <- function(X, Y, S, evalues, evectors, q=NULL, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1, cal="asym"){
  if(is.null(method)){
    # X[[i]] should be R realizations (n x R matrix)
    # p-dimensional graph signal
    p <- length(X)
    q <- length(Y)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    if(nrow(Y[[1]])!=nrow(S)){
      stop("The number of rows of list Y's element matrices should be equal to the number of rows of GSO!")
    }
    if(ncol(X[[1]])!=ncol(Y[[1]])){
      stop("The number of columns of list X's element matrices should be equal to the number of columns of list Y's element matrices!")
    }
    
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    P_array_X <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        P_array_X[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S) 
      }
    }
    
    P_array_Y <- array(NA, dim = c(n,q,q))
    for(i in 1:q){
      for(j in 1:q){
        P_array_Y[,i,j] <- cpsd.graph(x=Y[[i]], y=Y[[j]], S=S) 
      }
    }
    
    P_array_XY <- array(NA, dim = c(n,p,q))
    for(i in 1:p){
      for(j in 1:q){
        P_array_XY[,i,j] <- cpsd.graph(x=X[[i]], y=Y[[j]], S=S) 
      }
    }
    
    r <- min(p,q)
    H.hat <- array(NA, dim = c(n,r,p))
    F.hat <- array(NA, dim = c(n,r,q))
    eta.hat <- array(NA, dim = c(n,q,r))
    gamma.hat <- matrix(NA, nrow=n, ncol=r)
    
    P_array_YX <- array(NA, dim = c(n,q,p))
    for(l in 1:n){
      P_array_YX[l,,] <- Conj(t(P_array_XY[l,,]))
      if(cal=="sym"){
        targetmat <- solve(expm::sqrtm(P_array_Y[l,,])) %*% P_array_YX[l,,] %*% solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,]))
        tmp <- eigensort(targetmat)
        eta.hat[l,,] <- tmp$evectors[,q:(q-r+1)]
        gamma.hat[l,] <- rev(tmp$evalues)[1:r]
        
        # H.hat[l,,] <- t(solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        # H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        H.hat[l,,] <- t(solve(expm::sqrtm(P_array_X[l,,])) %*% P_array_XY[l,,] %*% solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        H.hat[l,,] <- t(solve(expm::sqrtm(P_array_X[l,,])) %*% t(H.hat[l,,]))
        
        F.hat[l,,] <- t(eta.hat[l,,])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
        F.hat[l,,] <- t(solve(expm::sqrtm(P_array_Y[l,,])) %*% t(F.hat[l,,]))
        
        # F.hat[l,,] <- t(solve(expm::sqrtm(P_array_Y[l,,])) %*% eta.hat[l,,])
        # F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
      } else if(cal=="asym"){
        targetmat1 <- solve(P_array_X[l,,]) %*% P_array_XY[l,,] %*% solve(P_array_Y[l,,]) %*% P_array_YX[l,,]
        targetmat2 <- solve(P_array_Y[l,,]) %*% P_array_YX[l,,] %*% solve(P_array_X[l,,]) %*% P_array_XY[l,,]
        
        tmp1 <- eigensort(targetmat1)
        gamma.hat[l,] <- rev(tmp1$evalues)[1:r]
        H.hat[l,,] <- t(tmp1$evectors[,q:(q-r+1)])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        
        tmp2 <- eigensort(targetmat2)
        F.hat[l,,] <- t(tmp2$evectors[,q:(q-r+1)])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
      }
    }
    
    
    Z <- matrix(0, nrow=n, ncol=r)
    W <- matrix(0, nrow=n, ncol=r)
    # Z[[i]] : canonical graph signal corresponding to X for R realizations (n x R matrix) / length r
    # W[[i]] : canonical graph signal corresponding to Y for R realizations (n x R matrix) / length r
    for(i in 1:r){
      for(j in 1:p){
        Z[,i] <- Z[,i] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X
      }
      for(j in 1:q){
        W[[i]] <- W[[i]] + evectors %*% diag(F.hat[,i,j]) %*% t(Conj(evectors)) %*% Y
      }
    }
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$Z <- Z ; res$W <- W
    res$eta.hat <- eta.hat ; res$gamma.hat <- gamma.hat ; res$H.hat <- H.hat ; res$F.hat <- F.hat
    res$Px <- P_array_X ; res$Py <- P_array_Y ; res$Pxy <- P_array_XY ; res$Pyx <- P_array_YX
    return(res)
  } 

  else if(method=="random"){
    # X[[i]] should be one realization (n x 1 vector)
    # p-dimensional graph signal
    p <- ncol(X)
    if(nrow(X)!=nrow(S)){
      stop("The length of list X's element vector OR the number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    
    q <- ncol(Y)
    if(nrow(Y)!=nrow(S)){
      stop("The number of rows of list Y's element matrices should be equal to the number of rows of GSO!")
    }
    
    # if(ncol(X)!=ncol(Y)){
    #   stop("The number of columns of list X's element matrices should be equal to the number of columns of list Y's element matrices!")
    # }
    
    n <- nrow(S)
    R <- 1
    
    # val <- eigensort(S)
    # evalues <- val$evalues
    # evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    
    m <- nrow(WB)
    out_X <- array(0, c(nrow(X), ncol(X), m)) # n x p x m
    out_Y <- array(0, c(nrow(Y), ncol(Y), m)) # n x q x m
    
    
    for (i in 1:m) {
      out_X[,,i] <- sweep(X, 1, WB[i, ], `*`)
      out_Y[,,i] <- sweep(Y, 1, WB[i, ], `*`)
    }
    
    
    print("calculating P_array_X...")
    P_array_X <- array(NA, dim = c(n,p,p))
    
    
    # apply evectors to all slices simultaneously
    out_X.tilde <- t(Conj(evectors)) %*% matrix(out_X, nrow=n, ncol=p*m)
    
    # reshape back into n x p x m
    out_X.tilde <- array(out_X.tilde, dim = c(n, p, m))
    
    for (k in 1:n) {
      P_array_X[k,,] <- out_X.tilde[k,,] %*% Conj(t(out_X.tilde[k,,])) / m
    }
    
    
    print("calculating P_array_Y...")
    P_array_Y <- array(NA, dim = c(n,q,q))
    # apply evectors to all slices simultaneously
    out_Y.tilde <- t(Conj(evectors)) %*% matrix(out_Y, nrow=n, ncol=q*m)
    
    # reshape back into n x p x m
    out_Y.tilde <- array(out_Y.tilde, dim = c(n, q, m))
    
    for (k in 1:n) {
      P_array_Y[k,,] <- out_Y.tilde[k,,] %*% Conj(t(out_Y.tilde[k,,])) / m
    }
    
    
    print("calculating P_array_XY, H and F...")
    P_array_XY <- array(NA, dim = c(n,p,q))
    P_array_YX <- array(NA, dim = c(n,q,p))
    # apply evectors to all slices simultaneously

    
    r <- min(p,q)
    H.hat <- array(NA, dim = c(n,r,p))
    F.hat <- array(NA, dim = c(n,r,q))
    eta.hat <- array(NA, dim = c(n,q,r))
    gamma.hat <- matrix(NA, nrow=n, ncol=r)
    
    
    for(l in 1:n){
      print(paste("iteration:", l))
      P_array_XY[l,,] <- out_X.tilde[l,,] %*% Conj(t(out_Y.tilde[l,,])) / m
      P_array_YX[l,,] <- Conj(t(P_array_XY[l,,]))
      if(cal=="sym"){
        targetmat <- ginv(Re(expm::sqrtm(P_array_Y[l,,]))) %*% P_array_YX[l,,] %*% ginv(Re(P_array_X[l,,])) %*% P_array_XY[l,,] %*% ginv(Re(expm::sqrtm(P_array_Y[l,,])))
        tmp <- eigensort(targetmat)
        eta.hat[l,,] <- tmp$evectors[,q:(q-r+1)]
        gamma.hat[l,] <- rev(tmp$evalues)[1:r]

        H.hat[l,,] <- t(ginv(Re(expm::sqrtm(P_array_X[l,,]))) %*% P_array_XY[l,,] %*% ginv(Re(expm::sqrtm(P_array_Y[l,,]))) %*% eta.hat[l,,])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))
        H.hat[l,,] <- t(ginv(Re(expm::sqrtm(P_array_X[l,,]))) %*% t(H.hat[l,,]))

        F.hat[l,,] <- t(eta.hat[l,,])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
        F.hat[l,,] <- t(ginv(Re(expm::sqrtm(P_array_Y[l,,]))) %*% t(F.hat[l,,]))
      } else if(cal=="asym"){
        targetmat1 <- ginv(Re(P_array_X[l,,])) %*% P_array_XY[l,,] %*% ginv(Re(P_array_Y[l,,])) %*% P_array_YX[l,,]
        targetmat2 <- ginv(Re(P_array_Y[l,,])) %*% P_array_YX[l,,] %*% ginv(Re(P_array_X[l,,])) %*% P_array_XY[l,,]

        tmp1 <- eigensort(targetmat1)
        gamma.hat[l,] <- rev(tmp1$evalues)[1:r]
        H.hat[l,,] <- t(tmp1$evectors[,p:(p-r+1)])
        H.hat[l,,] <- H.hat[l,,] / sqrt(rowSums(H.hat[l,,]^2))

        tmp2 <- eigensort(targetmat2)
        F.hat[l,,] <- t(tmp2$evectors[,q:(q-r+1)])
        F.hat[l,,] <- F.hat[l,,] / sqrt(rowSums(F.hat[l,,]^2))
      }
    }


    Z <- matrix(0, nrow=n, ncol=r)
    W <- matrix(0, nrow=n, ncol=r)
    # Z[[i]] : canonical graph signal corresponding to X for R realizations (n x R matrix) / length r
    # W[[i]] : canonical graph signal corresponding to Y for R realizations (n x R matrix) / length r
    VH_X <- t(Conj(evectors)) %*% X
    VH_Y <- t(Conj(evectors)) %*% Y
    
    for(i in 1:r){
      print(paste("i =", i))
      
      H_sub <- H.hat[,i,]   # n x p
      F_sub <- F.hat[,i,]   # n x q

      Z[,i] <- rowSums(evectors %*% (H_sub * VH_X))
      W[,i] <- rowSums(evectors %*% (F_sub * VH_Y))
    }
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$Z <- Z ; res$W <- W
    res$eta.hat <- eta.hat ; res$gamma.hat <- gamma.hat
    res$H.hat <- H.hat ; res$F.hat <- F.hat
    res$Px <- P_array_X ; res$Py <- P_array_Y ; res$Pxy <- P_array_XY ; res$Pyx <- P_array_YX
    
    return(res)
  }
}

gmcca <- function(X, Lg, gamm, M, d, eps = 1e-4) {
  # X: list of length M, X[[i]] = D_i x N matrix
  # Lg: N x N Laplacian
  # M: number of views
  # d: dimension of latent common sources
  # eps: small regularizer
  
  N <- ncol(X[[1]])  # number of samples
  SUMD <- matrix(0, N, N)
  
  # Precompute D_i
  D <- sapply(X, nrow)
  
  # SUMD = sum_i X_i' * (X_i X_i' + eps I)^(-1) * X_i
  for (i in 1:M) {
    temp <- X[[i]]                # D_i x N
    A <- temp %*% t(temp)         # (D_i x D_i)
    A_reg <- A + eps * diag(D[i]) # regularized
    
    # temp' / A_reg %*% temp   == temp' %*% solve(A_reg) %*% temp
    temp2 <- t(temp) %*% solve(A_reg, temp)
    SUMD <- SUMD + temp2
  }
  
  # Ms = SUMD - gamm * Lg
  Ms <- SUMD - gamm * Lg
  
  # eigen decomposition
  eig_res <- eigen(Ms)
  
  # sort eigenvalues descending
  ord <- order(eig_res$values, decreasing = TRUE)
  
  # St = N x d eigenvectors
  St <- eig_res$vectors[, ord[1:d]]
  
  # S = d x N
  S <- t(St)
  
  # Solve U_i
  U <- vector("list", M)
  for (i in 1:M) {
    temp <- X[[i]]
    A <- temp %*% t(temp)
    A_reg <- A + eps * diag(D[i])
    
    # U_i = (temp temp' + eps I)^(-1) temp S'
    U[[i]] <- solve(A_reg, temp %*% t(S))
  }
  
  return(list(U = U, S = S))
}
