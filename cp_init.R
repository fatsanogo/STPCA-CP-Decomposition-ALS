load("eigen_stpca_No3.RData")
#load("eigen_stpca_No4.RData")
#load("eigen_stpca_No6.RData")
load("my_data.RData")
library(fda); library(rTensor); library(R.matlab); library(R.matlab)

stpca_to_cp_init <- function(input_tensor, zw_C, K, R, period=NULL, time_points=NULL,
                             pick=c("abs","plus","minus")) {
  pick <- match.arg(pick)
  dims <- dim(input_tensor); Tlen <- dims[1]; S <- dims[2]; V <- dims[3]
  
  if (is.null(time_points)) time_points <- (1:Tlen)-0.5
  if (is.null(period)) period <- Tlen
  basis <- fda::create.fourier.basis(c(0, Tlen), nbasis=K, period=period)
  Phi   <- fda::eval.basis(time_points, basis)  # T × K
  
  # block starts: every K columns
  block_starts <- seq(1, K*V, by=K)
  
  # choose eigenpairs
  vals <- zw_C$values
  ord <- switch(pick,
                abs   = order(abs(vals), decreasing=TRUE),
                plus  = order(vals, decreasing=TRUE),
                minus = order(vals, decreasing=FALSE))
  idx <- ord[seq_len(R)]
  
  X <- rTensor::as.tensor(input_tensor)
  
  U_time  <- matrix(NA_real_, Tlen, R)
  U_space <- matrix(NA_real_, S,    R)
  U_var   <- matrix(NA_real_, V,    R)
  
  split_blocks <- function(gamma) {
    lapply(seq_len(V), function(v) {
      j0 <- block_starts[v]; gamma[j0:(j0+K-1)]
    })
  }
  
  for (r in seq_len(R)) {
    gamma  <- zw_C$vectors[, idx[r]]       # length K·V
    pieces <- split_blocks(gamma)
    
    # Reconstruct T×V slab from Fourier coeffs
    Tmat <- matrix(0, Tlen, V)
    for (v in 1:V) Tmat[,v] <- as.vector(Phi %*% matrix(pieces[[v]], ncol=1))
    
    # Best rank-1 across (time × variable)
    sv <- svd(Tmat, nu=1, nv=1)
    u_time <- sv$u[,1]
    u_var  <- sv$v[,1]
    
    # Spatial factor by projecting data given (u_time, u_var)
    tmp <- ttl(X, list(matrix(u_time, nrow=1), matrix(u_var, nrow=1)), ms=c(1,3)) # 1×S×1
    u_space <- as.vector(tmp@data)
    
    # unit-norm columns
    U_time[,r]  <- u_time  / sqrt(sum(u_time^2))
    U_var[,r]   <- u_var   / sqrt(sum(u_var^2))
    U_space[,r] <- u_space / sqrt(sum(u_space^2))
  }
  
  list(U_time=U_time, U_space=U_space, U_var=U_var)
}


X <- simplify2array(my_data)
init_res <- stpca_to_cp_init(
  input_tensor = X,
  R    = 2,          # your target CP rank
  zw_C = zw_C,
  K    = 25,
  period = 365.2425, # daily seasonal;
  pick = "abs"      
)
# save for R=3
#save(init_res, file = "stpca_init3.RData")
writeMat(con = "stpca_init3.mat", 
         U_time = init_res$U_time,
         U_space = init_res$U_space,
         U_var = init_res$U_var)
# save for R=2
#save(init_res, file = "stpca_init2.RData")
writeMat(con = "stpca_init2.mat", 
         U_time = init_res$U_time,
         U_space = init_res$U_space,
         U_var = init_res$U_var)
