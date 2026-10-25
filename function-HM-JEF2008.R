# HM estimation of factor models for HF returns with EIV
#
# The function performs the HM estimation procedure described in Carmichael
# and Coen (2008, JEF), mathematical details provided in the appendix
#
fit_HM_JEF <- function(
      Y_vec,       # n x 1 vector of dependent variable
      F_mat,       # n x k matrix of explanatory variables
      Z_mat        # matrix of HM instruments
      )
  {
   #
   Y_vec <- as.matrix(Y_vec)
   F_mat <- as.matrix(F_mat)
   Z_mat <- as.matrix(Z_mat)
   #
   # check the nb of non-empty rows in Y and F
   if(dim(Y_vec)[1] != dim(F_mat)[1]){
      stop("Y and F must have the same number of rows.")
   }
   # nb of observations
   nb_obs <- as.numeric(dim(Y_vec)[1])
   #
   # calculate PZ_mat
   PZ_mat <- Z_mat %*% solve(t(Z_mat) %*% Z_mat) %*% t(Z_mat)
   #
   # calculate MZ_mat
   I_mat  <- diag(nb_obs)
   MZ_mat <- I_mat - Z_mat %*% solve(t(Z_mat) %*% Z_mat) %*% t(Z_mat)
   #
   # W_mat
   W_mat  <- MZ_mat %*% F_mat
   #
   # estimating theta coefficients
   X_mat <- cbind(rep(1, nb_obs), F_mat)
   theta <- solve(t(X_mat) %*% PZ_mat %*% X_mat) %*%
            t(X_mat) %*% PZ_mat %*% Y_vec
   #
   # estimating psi coefficients
   MX_mat<- I_mat - X_mat %*% solve(t(X_mat) %*% X_mat) %*% t(X_mat)
   psi   <- solve(t(F_mat) %*% MZ_mat %*% MX_mat %*% MZ_mat %*% F_mat) %*%
            t(F_mat) %*% MZ_mat %*% MX_mat %*% Y_vec
   #
   # predicted values, residuals, and residual SE
   pred_Y_vec <- X_mat %*% theta + W_mat %*% psi
   resid_HM <- Y_vec - pred_Y_vec
   nb_predictor <- as.numeric(length(theta)) - 1 + as.numeric(length(psi))
   sigma_HM <- sqrt(1 / (nb_obs - nb_predictor) * sum(resid_HM^2))
   #
   # Omega matrix
   Omega_res <- Y_vec - X_mat %*% theta
   Omega <- diag(c(Omega_res^2),
                 nrow = nb_obs,
                 ncol = nb_obs)
   #
   # labels
   labels_theta <- c("theta0",
                     paste("theta",
                           1:as.numeric(dim(F_mat)[2]),
                           sep = ""))
   labels_psi   <- c(paste("psi",
                           1:as.numeric(dim(F_mat)[2]),
                           sep = ""))
   #
   # VCV matrix theta coefs 
   theta_VCV <- solve(t(X_mat) %*% PZ_mat %*% X_mat) %*%
      (t(X_mat) %*% PZ_mat %*% Omega %*% PZ_mat %*% X_mat) %*%
      solve(t(X_mat) %*% PZ_mat %*% X_mat)
   theta_SE  <- diag(theta_VCV)^0.5
   #
   # VCV matrix psi coefs
   psi_VCV <- solve(t(F_mat) %*% MZ_mat %*% MX_mat %*% MZ_mat %*% F_mat) %*%
      (t(F_mat) %*% MZ_mat %*% MX_mat %*% Omega %*% MX_mat %*% MZ_mat %*% F_mat) %*%
      solve(t(F_mat) %*% MZ_mat %*% MX_mat %*% MZ_mat %*% F_mat)
   psi_SE <- diag(psi_VCV)^0.5
   #
   # DWH test
   DWH_stat <- t(psi) %*%
      (t(F_mat) %*% MZ_mat %*% MX_mat %*% MZ_mat %*% F_mat) %*%
      solve(t(F_mat) %*% MZ_mat %*% MX_mat %*% Omega %*%
               MX_mat %*% MZ_mat %*% F_mat) %*%
      (t(F_mat) %*% MZ_mat %*% MX_mat %*% MZ_mat %*% F_mat) %*% psi
   DWH_pval <- 1 - pchisq(DWH_stat, as.numeric(length(psi)))
   #
   # collecting output
   ## theta coefs
   theta_fit <- matrix(theta, ncol = 1)
   theta_SE  <- matrix(theta_SE, ncol = 1)
   theta_t   <- matrix(theta / theta_SE, ncol = 1)
   theta_out <- cbind(theta_fit, theta_SE, theta_t)
   rownames(theta_out) <- c(labels_theta)
   colnames(theta_out) <- c("est", "SE", "t-stat")
   ## psi (e.g. EIV) coeffs
   psi_fit <- matrix(psi, ncol = 1)
   psi_SE  <- matrix(psi_SE, ncol = 1)
   psi_t   <- matrix(psi / psi_SE, ncol = 1)
   psi_out <- cbind(psi_fit, psi_SE, psi_t)
   rownames(psi_out) <- c(labels_psi)
   colnames(psi_out) <- c("est", "SE", "t-stat")
   ## DWH test
   DWH_test <- matrix(c(DWH_stat, DWH_pval), ncol = 1)
   rownames(DWH_test) <- c("statistic", "p-value")
   colnames(DWH_test) <- c("DWH test, H0: No EIV")
   #
   # report output
   return(list(theta    = theta_out,
               psi      = psi_out,
               sigma    = sigma_HM,
               DWH_test = DWH_test)
          )
   }


