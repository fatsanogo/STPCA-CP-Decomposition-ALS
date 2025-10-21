source("STPCA.r")
load("my_data.RData")
#load("boundary_PL.RData")
load("new_england_centroid_dist_km.RData")

# The first step: functional PCA
data <- data_rang(my_data)
A <- matrix_A(data,K=25)
A <- scale(A,scale=F)
n <- nrow(A)
S <- ((n-1)/n)*var(A)

# Eigenvalues of the matrix $\pmb{S}=\hat{\pmb{\Sigma}}$
zw_S <- eigen(S) 
#round(zw_S$values,6)

# Proportion of variability: 
matrix_M(S,K=3)

# First functional PCA
contribution(S,Var_number = 1)

# Second functional PCA
contribution(S,Var_number = 2)

# spatio-temporal PCA

# Spatial weight matrix No 3; Radial Distance Weights
  # Eigenvalues of the matrix $\pmb{C}$
W <- matrix_W_dist(D_km,type=3,param=4)
C_stpca <- matrix_Moran(A,W)
zw_C <- eigen(C_stpca)
#round(zw_C$values,6)
save(zw_C, A, file = "eigen_stpca_No3.RData")

# Positive spatio-temporal principal components
# Proportion of variability: 
matrix_M(C_stpca,K=3,param="plus")

# Length of the subvectors:
 # First STPCA
contribution(C_stpca,Var_number = 1)

# Second STPCA
contribution(C_stpca,Var_number = 2)

# Negative spatio-temporal principal components
# Proportion of variability: 
matrix_M(C_stpca,K=3,param="minus")

# Length of the subvectors:
  # First STPCA
contribution(C_stpca,Var_number = 39)

# Second STPCA
contribution(C_stpca,Var_number = 38)

# Spatial weight matrix No 4; Power Distance Weights; $\alpha=2$**
  # Eigenvalues of the matrix $\pmb{C}$
W <- matrix_W_dist(D_km,type=4,param=2)
C <- matrix_Moran(A,W)
zw_C <- eigen(C)
#round(zw_C$values,6)
save(zw_C, A, file = "eigen_stpca_No4.RData")

# Positive spatio-temporal principal components
# Proportion of variability: 
matrix_M(C,K=3,param="plus")

# Length of the subvectors:
  # First STPCA
contribution(C,Var_number = 1)

# Second STPCA
contribution(C,Var_number = 2)

# Negative spatio-temporal principal components
# Proportion of variability: 
matrix_M(C,K=3,param="minus")

# Length of the subvectors:
  # First STPCA
contribution(C,Var_number = 39)

# Second STPCA
contribution(C,Var_number = 38)

# Spatial weight matrix No 6; Double-Power Distance Weights; k=3
  # Eigenvalues of the matrix $\pmb{C}$
W <- matrix_W_dist(D_km,type=6,param=3)
C <- matrix_Moran(A,W)
zw_C <- eigen(C)
#round(zw_C$values,6)
save(zw_C, A, file = "eigen_stpca_No6.RData")

# Positive spatio-temporal principal components
# Proportion of variability: 
matrix_M(C,K=3,param="plus")

# Length of the subvectors:
  # First STPCA
contribution(C,Var_number = 1)

# Second STPCA
contribution(C,Var_number = 2)

# Negative spatio-temporal principal components
# Proportion of variability: 
matrix_M(C,K=3,param="minus")

# Length of the subvectors:
  # First STPCA
contribution(C,Var_number = 39)

# Second STPCA
contribution(C,Var_number = 38)

