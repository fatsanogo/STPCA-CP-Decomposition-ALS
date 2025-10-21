library(fda)
library(ggplot2)
library(ggthemes)
library(sf)
library(dplyr)
library(matrixStats)

data_rang <- function(data, low=1, high=7,
                             use_quantiles=TRUE, q_lo=0.01, q_hi=0.99,
                             eps=1e-12) {
  out <- vector("list", length(data)); names(out) <- names(data)
  for (i in seq_along(data)) {
    M <- data[[i]]
    if (use_quantiles) {
      lo <- rowQuantiles(M, probs=q_lo, na.rm=TRUE)
      hi <- rowQuantiles(M, probs=q_hi, na.rm=TRUE)
      M  <- pmin(pmax(M, lo), hi)  # winsorize per row
      a <- lo; b <- hi
    } else {
      a <- rowMins(M, na.rm=TRUE); b <- rowMaxs(M, na.rm=TRUE)
    }
    rng <- pmax(b - a, eps)
    Z <- sweep(M, 1, a, FUN="-")
    Z <- sweep(Z, 1, rng, FUN="/")
    Z <- low + (high - low) * Z
    Z[(b - a) < eps, ] <- (low + high)/2  # flat rows → center (4)
    rownames(Z) <- rownames(M); colnames(Z) <- colnames(M)
    out[[i]] <- Z
  }
  out
}

matrix_A <- function(Data,K=25)
{
  p <- length(Data)
  T <- dim(Data[[1]])[1]
  time <- (1:T)-0.5
  A <- NULL
  for (i in 1:p)
  {
    dane <- Data[[i]]
    baza <- create.fourier.basis(c(0,T),K, period = 365.2425)
    danefd <- smooth.basis(time, dane, baza)
    X <- danefd$fd
    C <- X$coefs
    C <- t(C)
    A <- cbind(A,C)
  }
  return(A)
}

matrix_M <- function(MS,K=4,param=NULL)
{
  M <- matrix(0,ncol=2,nrow=K)
  zw_S <- eigen(MS) 
  wart <- zw_S$values
  if (is.null(param))
  {
    suma <- sum(wart)
    for (i in 1:K)
    {
      M[i,1] <- wart[i]/suma
      M[i,2] <- sum(wart[1:i])/suma
    }
    rownames(M) <- 1:K
  }
  else if (param=='plus')
  {
    suma <- sum(wart[wart>0])
    for (i in 1:K)
    {
      M[i,1] <- wart[i]/suma
      M[i,2] <- sum(wart[1:i])/suma
    }
    rownames(M) <- 1:K
  }
  else if (param=='minus')
  {
    suma <- sum(wart[wart<0])
    n <- length(wart)
    for (i in n:(n-K+1))
    {
      M[n-i+1,1] <- wart[i]/suma
      M[n-i+1,2] <- sum(wart[n:i])/suma
    }
    rownames(M) <- n:(n-K+1)
  }
  M <- M*100
  colnames(M) <- c("Proportion","Cumulative")
  return(M) 
}

contribution <- function(MS,Var_number=1,K=25)
{
  zw_S <- eigen(MS) 
  gamma <- zw_S$vectors[,Var_number] 
  p <- length(gamma)/K
  M <- matrix(0,ncol=2,nrow=p)
  for (i in 1:p)
  {
    M[i,1] <- sqrt(sum(gamma[((i-1)*K+1):(i*K)]^2))
    M[i,2] <- sum(gamma[((i-1)*K+1):(i*K)]^2)*100
  }
  rownames(M) <- 1:p
  colnames(M) <- c("length","%")
  return(M)
}

matrix_W_dist <- function(MW,type=1,param=4)
{
  if (type==1)
  {
  N <- ncol(MW)
  MOW <- unlist(MW)
  MOW <- as.numeric(MOW)
  MOWS <- sort(MOW)
  MOWS <- unique(MOWS)
  MOWS <- MOWS[-1]
  W <- matrix(0,nrow=N,ncol=N)
  for (i in 1:N)
  {
    MOi <- as.numeric(MW[i,])
    MOiS <- sort(MOi)
    MOiS <- MOiS[-1]
    for (j in 1:N)
    {
      if (MW[i,j]<=MOiS[param]) W[i,j] <- 1
    }
  }
  for (i in 1:N) W[i,i] <- 0
  }
  if (type==2)
  {
    N <- ncol(MW)
    W <- matrix(0,nrow=N,ncol=N)
    for (i in 1:N)
    {
      MOi <- as.numeric(MW[i,])
      MOiS <- sort(MOi)
      MOiS <- MOiS[-1]
      for (j in 1:N)
      {
        if (MW[i,j]<=MOiS[param]) W[i,j] <- 1
      }
    }
    for (j in 1:N)
    {
      MOj <- as.numeric(MW[,j])
      MOjS <- sort(MOj)
      MOjS <- MOjS[-1]
      for (i in 1:N)
      {
        if (MW[i,j]<=MOjS[param]) W[i,j] <- 1
      }
    }
    for (i in 1:N) W[i,i] <- 0
  }
  if (type==3)
  {
    N <- ncol(MW)
    W <- matrix(0,nrow=N,ncol=N)
    MOL <- as.dist(MW)
    d <- mean(MOL)+sd(MOL)
    for (i in 1:N)
      for (j in 1:N)
      {
        if (MW[i,j]<=d) W[i,j] <- 1
      }
    for (i in 1:N) W[i,i] <- 0
  }
  if (type==4)
  {
    N <- ncol(MW)
    W <- matrix(0,nrow=N,ncol=N)
    for (i in 1:N)
      for (j in 1:N)
      {
        W[i,j] <- MW[i,j]^(-param)
      }
    for (i in 1:N) W[i,i] <- 0
  }
  if (type==5)
  {
    N <- ncol(MW)
    W <- matrix(0,nrow=N,ncol=N)
    MO <- matrix(0,nrow=N,ncol=N)
    suma <- sum(MW)/(N*(N-1))
    MO <- MW/suma
    for (i in 1:N)
      for (j in 1:N)
      {
        W[i,j] <- exp(-param*MO[i,j])
      }
    for (i in 1:N) W[i,i] <- 0
  }
  if (type==6)
  {
    N <- ncol(MW)
    W <- matrix(0,nrow=N,ncol=N)
    MOL <- as.dist(MW)
    d <- mean(MOL)+sd(MOL)
    for (i in 1:N)
      for (j in 1:N)
      {
        if (MW[i,j]<=d) W[i,j] <- (1-(MW[i,j]/d)^param)^param
      }
    for (i in 1:N) W[i,i] <- 0
  }
  return(W)
}

matrix_W_bound <- function(MW,type=1)
{
  if (type==1)
  {
    N <- ncol(MW)
    for (i in 1:N) MW[i,i] <- 0
    W <- matrix(0,nrow=N,ncol=N)
    for (i in 1:N)
    {
      for (j in 1:N)
      {
        if (MW[i,j]>0) W[i,j] <- 1
      }
    }
    for (i in 1:N) W[i,i] <- 0
  }
  if (type==2)
  {
    N <- ncol(MW)
    for (i in 1:N) MW[i,i] <- 0
    W <- matrix(0,nrow=N,ncol=N)
    for (i in 1:N)
    {
      for (j in 1:N)
      {
        W[i,j] <- MW[i,j]/sum(MW[i,])
      }
    }
    for (i in 1:N) W[i,i] <- 0
  }
  return(W)
}

matrix_W_dist_bound <- function(MWD,MWB)
{
  N <- ncol(MWD)
  MWP <- matrix(0,nrow=N,ncol=N)
  for (i in 1:N)
    for (j in 1:N)
    {
      MWP[i,j] <- MWD[i,j]^(-1)
    }
  for (i in 1:N) MWP[i,i] <- 0
  for (i in 1:N) MWB[i,i] <- 0
    W <- matrix(0,nrow=N,ncol=N)
    for (i in 1:N)
      for (j in 1:N)
      {
        W[i,j] <- (MWB[i,j]*MWP[i,j])/sum(MWB[i,]*MWP[i,])
      }
  for (i in 1:N) W[i,i] <- 0
  return(W)
}

matrix_Moran <- function(MA,MW)
{
  W <- MW
  W <- (1/2)*(W+t(W))
  W <- W/sum(W)
  E <- t(MA)%*%W%*%MA
  return(E)
}

ggprojection <- function(MA,MS,No1=1,No2=2,OX=1,OY=1)
{
  zw_S <- eigen(MS) 
  m <- MA%*%zw_S$vectors[,c(No1,No2)]
  m <- scale(m,scale=F)
  m <- as.data.frame(m)
  m[,2] <- OX*m[,2]
  m[,1] <- OY*m[,1]
  N <- nrow(m)
  gg <- ggplot(data = m,aes(x=m[,1],y=m[,2],label=1:N)) +
    geom_point(pch=21,size=9) +
    geom_text(size=4) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    labs(title='',x=bquote(U[1]),y=bquote(U[2])) +
    theme_bw() + 
    theme(text=element_text(size=14))
  gg
}

ggchoro <- function(MA,MS,No1,No2,mapa,OX=1,OY=1)
{
  zw_S <- eigen(MS) 
  m <- MA%*%zw_S$vectors[,c(No1,No2)]
  m <- scale(m,scale=F)
  m[,2] <- OX*m[,2]
  m[,1] <- OY*m[,1]
  n <- nrow(MA)
  ct <- vector("numeric",n)
  for (i in 1:n)
  {
    if (m[i,1]>=0 & m[i,2]>=0) 
    {
      ct[i] <- 1
    } else if (m[i,1]>=0 & m[i,2]<0) 
    {
      ct[i] <- 2
    } else if (m[i,1]<0 & m[i,2]>=0) 
    {
      ct[i] <- 3
    }
    else if (m[i,1]<0 & m[i,2]<0) 
    {
      ct[i] <- 4
    }
  }
  ggplot(data = mapa) +
    geom_sf(fill=ct) +
    theme_tufte() +
    theme_void()
}

