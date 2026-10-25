# HM estimation of factor models for HF returns with EIV
#
# The function returns the Z-matrix of higher-moments proposed by
# Dagenais and Dagenais (1997, JE)
#
generate_Z_mat <- function(
    X_mat,		        # n by k matrix of independent variables
    Y_vec = NULL,     # n by 1 vector of dependent variable
    n,                # sample size 
    k,                # nb of predictors
    iv_list = c(1, 4) # list of higher-moments
)
{
  #
  # nb of cols in X_mat must be equal to k
  if(k != dim(as.matrix(X_mat))[2]){
    stop("The number of columns in X_mat must be equal to k.")
  }
  #
  # calculate Y_vec and X_mat in mean deviations form
  y_vec <- scale(Y_vec, scale = FALSE)
  x_mat <- scale(X_mat, scale = FALSE)
  #
  # calculate higher-moment instruments
  #
  # z0
  # calculate z0 = vector of ones
  Z_mat <- matrix(rep(1, n), nrow = n, ncol = 1)
  # label z0
  colnames(Z_mat) <- "z0"
  #
  # column labels index for additional instruments
  col_index <- 1
  #
  # z1
  if(is.element(1, iv_list)){
    # calculate z1
    z1 <- x_mat * x_mat
    # append Z_mat
    Z_mat <- cbind(Z_mat, z1)
    # label z1
    new_cols <- paste0("z1", 1:k)
    colnames(Z_mat)[(col_index + 1):(col_index + k)] <- new_cols
    col_index <- col_index + k
  }
  # z2
  if(is.element(2, iv_list)){
    # calculate z2
    z2 <- matrix(0, nrow = n, ncol = k)
    for(j in 1:k){
      z2[,j] <- x_mat[,j] * y_vec
    }
    # append Z_mat
    Z_mat <- cbind(Z_mat, z2)
    # label z2
    new_cols <- paste0("z2", 1:k)
    colnames(Z_mat)[(col_index + 1):(col_index + k)] <- new_cols
    col_index <- col_index + k
  }
  # z3
  if((is.element(3, iv_list))){
    # calculate z3
    z3 <- y_vec * y_vec
    # append Z_mat
    Z_mat <- cbind(Z_mat, z3)
    # label z3
    colnames(Z_mat)[col_index + 1] <- "z31"
    col_index <- col_index + 1
  }
  # z4
  if((is.element(4, iv_list))){
    # calculate z4
    # z4_1 and z4_2
    z4_1 <- x_mat * x_mat * x_mat
    z4_2 <- 3 * x_mat %*% (((t(x_mat) %*% x_mat) / n) * diag(k))
    #
    z4 <- z4_1 - z4_2
    # append Z_mat
    Z_mat <- cbind(Z_mat, z4)
    # label z4
    new_cols <- paste0("z4", 1:k)
    colnames(Z_mat)[(col_index + 1):(col_index + k)] <- new_cols
    col_index <- col_index + k
  }
  # z5
  if((is.element(5, iv_list))){
    # calculate z5
    # z5_1
    z5_1 <- matrix(0, nrow = n, ncol = k)
    for(j in 1:k){
      z5_1[,j] <- x_mat[,j] * x_mat[,j] * y_vec
    }
    # z5_2
    z5_2 <- 2 * (x_mat %*% 
                   diag(
                     as.vector(
                       (1/n) * (t(x_mat) %*% y_vec)
                     ),
                     nrow = k, ncol = k
                   )
    )
    # z5_3
    # z5_3_1 and z5_3_2
    z5_3_1 <- matrix(1, nrow = k, ncol = 1)
    z5_3_2 <- ((t(x_mat) %*% x_mat) / n) * diag(k)
    # 
    z5_3 <- y_vec %*% (t(z5_3_1) %*% z5_3_2)
    # 
    z5 <- z5_1 - z5_2 - z5_3
    # append Z_mat
    Z_mat <- cbind(Z_mat, z5)
    # label z5
    new_cols <- paste0("z5", 1:k)
    colnames(Z_mat)[(col_index + 1):(col_index + k)] <- new_cols
    col_index <- col_index + k
  }
  # z6
  if((is.element(6, iv_list))){
    # calculate z6
    # z6_1
    z6_1 <- matrix(0, nrow = n, ncol = k)
    for(j in 1:k){
      z6_1[,j] <- x_mat[,j] * y_vec * y_vec
    }
    # z6_2
    scalar <- (1/n) * (t(y_vec) %*% y_vec)
    z6_2 <- x_mat * as.numeric(scalar)
    # z6_3
    z6_3 <- 2 * (
      y_vec %*% ((1/n) * (t(y_vec) %*% x_mat))
    )
    #
    z6 <- z6_1 - z6_2 - z6_3
    # append Z_mat
    Z_mat <- cbind(Z_mat, z6)
    # label z6
    new_cols <- paste0("z6", 1:k)
    colnames(Z_mat)[(col_index + 1):(col_index + k)] <- new_cols
    col_index <- col_index + k
  }
  # z7
  if((is.element(7, iv_list))){
    # calculate z7
    z7 <- matrix(0, nrow = n, ncol = 1)
    #
    z7_1 <- y_vec * y_vec * y_vec
    #
    scalar <- (1/n) * (t(y_vec) %*% y_vec)
    z7_2 <- 3 * y_vec * as.numeric(scalar)
    #
    z7 <- z7_1 - z7_2
    # append Z_mat
    Z_mat <- cbind(Z_mat, z7)
    # label z7
    colnames(Z_mat)[col_index + 1] <- "z71"
    col_index <- col_index + 1
  }
  #
  return(Z_mat)
}


