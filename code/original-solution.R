###################################################################
########################## SepBART_blk ############################
###################################################################
##### Instead of calculating the inverse of an N by N matrix, we generate multiple 100 by 100 matrices and solve them.
###################################################################
########################### Arguments #############################
###################################################################
##### Y: Outcome
##### X: n by p covariate matrix
##### W: q dimensional exposure matrix
##### Xtest: Testing locations (not be used in data analysis)
##### W1, W0: Two exposure levels under which we want to compare the outcomes
##### n_blocks: The number of blocks which we want to decompose an n by n large matrix into.
##### nMCMC: The total number of MCMC iterations
##### BurnIn_portion: The proportion of iterations to be discarded
##### stepsize: The step size after the burn-in period.
##### Other arguments control the tree parameters (not be used in data analysis)
###################################################################

SepBART_blk <- function(Y,X,W, Xtest = NULL, W1 = NULL, W0 = NULL, W0_quantile = NULL,
                        n_blocks = NULL,
                        Mod_ind = NULL,
                        nMCMC=2000,BurnIn_portion=0.2,stepsize=5){
  require("SoftBart")
  require("ModBart")
  require("mvtnorm")
  start.time.MCMC <- Sys.time()
  
  num_obs <- length(Y) ## Sample size
  p <- ncol(X) ## The number of covariates
  q <- ncol(W) ## The number of exposures
  BurnIn <- nMCMC*BurnIn_portion ## The number of MCMC iterations to be discarded.
  save_index <- seq(BurnIn+1,nMCMC,stepsize) ## A sequence of iteration indices to be stored.
  
  ## If n_blocks is not given,
  ## find the total number of blocks necessary to partition into 100 by 100 sub-matrices.
  if (is.null(n_blocks)) n_blocks <- floor(num_obs/100)
  
  if (is.null(Mod_ind)) Mod_ind <- 1:p
  
  X_names <- colnames(X)
  W_names <- colnames(W)
  Mod_names <- colnames(X)[Mod_ind]
  
  ## We will normalize the covariates and exposures to fall between 0 and 1 using quantiles.
  ## Define the quantile function.
  trank <- function(x) {
    x_unique <- unique(x)
    x_ranks <- rank(x_unique, ties.method = "max")
    tx <- x_ranks[match(x,x_unique)] - 1
    
    tx <- tx / length(unique(tx))
    tx <- tx / max(tx)
    
    return(tx)
  }
  quantile_normalize_bart <- function(X) {
    apply(X = X, MARGIN = 2, trank)
  }
  
  ## We scale the outcome and store the original standard deviation.
  sd_Y <- sd(Y)
  Y_orig <- Y
  Y <- Y/sd_Y
  
  ## Save the quantiles of exposures and covariates
  X_quantiles_orig <- list()
  W_quantiles_orig <- list()
  for (j in 1:p){
    X_quantiles_orig[[j]] <- quantile(X[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
  }
  for (j in 1:q){
    W_quantiles_orig[[j]] <- quantile(W[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
  }
  
  ## default W1 and W0
  if (is.null(W1)) W1 = apply(W, 2, quantile, 0.75, na.rm = TRUE)
  if (is.null(W0) & is.null(W0_quantile)){
    W0 = apply(W, 2, quantile, 0.25, na.rm = TRUE)
  }else if (is.null(W0) & !is.null(W0_quantile)){
    W0 = apply(W, 2, quantile, W0_quantile, na.rm = TRUE)
  }
  
  
  ## Normalize the covariates and exposures.
  if (is.null(Xtest)==FALSE){
    norm_X_Xtest <- quantile_normalize_bart(rbind(X,Xtest))
    norm_W_W0W1 <- quantile_normalize_bart(rbind(W,W0,W1))
    X <- norm_X_Xtest[1:num_obs,]
    Xtest <- norm_X_Xtest[-(1:num_obs),]
    W <- norm_W_W0W1[1:num_obs,]
    W0 <- norm_W_W0W1[num_obs+1,]
    W1 <- norm_W_W0W1[num_obs+2,]
  }else{
    norm_X <- quantile_normalize_bart(X)
    norm_W_W0W1 <- quantile_normalize_bart(rbind(W,W0,W1))
    X <- norm_X[1:num_obs,]
    W <- norm_W_W0W1[1:num_obs,]
    W0 <- norm_W_W0W1[num_obs+1,]
    W1 <- norm_W_W0W1[num_obs+2,]
  }
  
  ## W quantiles after normalizing. This will be used to calculate E(Y(W)-Y(W0))
  X_sample_size <- min(500, num_obs)
  W_quantiles <- list()
  for (j in 1:q){
    W_quantiles[[j]] <- quantile(W[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
  }
  W_quantiles_length <- length(W_quantiles[[j]])
  W1_rep <- t(W1)[rep(1,num_obs),] ## n by q matrix with each row of W1. Use this to calculate h(X,W1)
  W0_rep <- t(W0)[rep(1,num_obs),] ## n by q matrix with each row of W0. Use this to calculate h(X,W0)
  if (is.null(Xtest)==FALSE){
    W1_rep_test <- t(W1)[rep(1,dim(Xtest)[1]),]
    W0_rep_test <- t(W0)[rep(1,dim(Xtest)[1]),]
  }
  
  
  ## Our model is
  ## E(Y|X) = f(X) + G(W) + h(X,W) where h(X,W) = sum{h_j(X_j,W)} where X_j denotes the j-th covariate.
  ## To make the model identifiable, we define X* and W* as the reference vectors
  ## so that f(X*) = G(W*) = 0 and h(X*, . ) = h( . , W*) = 0.
  ## X* and W* can be considered as the baseline of the model.
  ## We define X* and W* to be the sample (column) means of X and W.
  X_star <- apply(X,2,function(x) sort(x)[floor(num_obs/2)])
  W_star <- apply(W,2,function(x) sort(x)[floor(num_obs/2)])
  
  ## The interaction term h(.,.) equals 0 for any rows of three matrices below. 
  ## n by (p+q) matrix with the i-th row of (W_i, X*) where W_i denotes the i-th observation's exposure vector.
  obs_grid_x0 <- 0
  for (i in 1:num_obs){
    obs_grid_x0 <- rbind(obs_grid_x0, c(W[i,],X_star))
  }
  obs_grid_x0 <- obs_grid_x0[-1,]
  ## Similarly define the n by (p+q) matrix with the i-th row of (W*, X_i).
  obs_grid_w0 <- 0
  for (i in 1:num_obs){
    obs_grid_w0 <- rbind(obs_grid_w0, c(W_star,X[i,]))
  }
  obs_grid_w0 <- obs_grid_w0[-1,]
  ## n by (p+q) matrix with each row of (W*,X*)
  obs_grid_00 <- matrix(rep(c(W_star,X_star),num_obs), nrow = num_obs, byrow = T)
  
  block_index <- sample(rep(1:n_blocks, len=num_obs))
  X_sampled <- list()
  W_sampled <- list()
  X_sampled_grid <- list()
  W_sampled_grid <- list()
  W0_grid <- list()
  W_star_sampled <- list()
  X_star_sampled <- list()
  W0_sampled <- list()
  for (block in 1:n_blocks){
    sample_index <- which(block_index==block)
    N <- sum(block_index==block)
    X_sampled[[block]] <- X[sample_index,]
    W_sampled[[block]] <- W[sample_index,]
    X_sampled_grid[[block]] <- X_sampled[[block]][rep(1:N,N),]
    W_sampled_grid[[block]] <- W_sampled[[block]][rep(1:N,each=N),]
    W0_grid[[block]] <- t(W0)[rep(1,N^2),]
    W_star_sampled[[block]] <- t(W_star)[rep(1, N),]
    X_star_sampled[[block]] <- t(X_star)[rep(1, N),]
    W0_sampled[[block]] <- t(W0)[rep(1, N),]
  }
  
  ## Make 3+p matrices to record the function values evaluated at n observations during MCMC scans.
  ## One matrix for each constant, f(X), g(W), and p matrices for h(X,W).
  mu_hat_c <- matrix(NA, nrow=nMCMC, ncol=num_obs)
  mu_hat_X <- matrix(NA, nrow=nMCMC, ncol=num_obs)
  mu_hat_W <- matrix(NA, nrow=nMCMC, ncol=num_obs)
  mu_hat_intact <- list()
  for (j in Mod_ind){
    mu_hat_intact[[j]] <- matrix(NA, nrow=nMCMC, ncol=num_obs)
  }
  
  ## Set the starting values
  mu_hat_c[1,] <- 0
  mu_hat_X[1,] <- rep(0,num_obs)
  mu_hat_W[1,] <- rep(0,num_obs)
  for (j in Mod_ind){
    mu_hat_intact[[j]][1,] <- rep(0,num_obs)
  }
  
  ## We will also calculate
  ## CATE_hat: E(Y(W1)-Y(W0)|X)
  ## psi_j = (Heterogeneity without X_j)/(Total heterogeneity): Variable importance measures
  if (is.null(Xtest)==FALSE){
    Xtest <- as.matrix(Xtest)
    CATE_test_hat <- E_Y_W0_test <- E_Y_W1_test <- matrix(NA, nrow=nMCMC, ncol=dim(Xtest)[1])
  }else{CATE_test_hat <- E_Y_W0_test <- E_Y_W1_test <- NA}
  CATE_hat <- E_Y_W1 <- E_Y_W0 <- matrix(NA, nrow=nMCMC, ncol=num_obs)
  Hetero_total <- rep(NA, nMCMC)
  Hetero_wo <- list()
  psi_hat <- list()
  CATE_h_hat <- list()
  for (j in Mod_ind){
    Hetero_wo[[j]] <- rep(NA, nMCMC)
    psi_hat[[j]] <- rep(NA, nMCMC)
    CATE_h_hat[[j]] <- matrix(NA, nrow=nMCMC, ncol=num_obs)
  }
  
  
  ## Estimate the conditional distribution X_j|X_{-j} using the normality assumption
  ## This will be used to calculate the heterogeneity without the j-th covariate.	
  normal_weight <- list()
  for (j in Mod_ind){
    normal_cd <- lm(X[,j]~X[,-j])
    normal_weight[[j]] <- list()
    for (block in 1:n_blocks){
      N <- sum(block_index==block)
      normal_weight[[j]][[block]] <- matrix(NA, nrow = N, ncol = N)
      for (w in 1:N){
        normal_weight[[j]][[block]][w,] <- dnorm(X_sampled[[block]][,j], mean = c(1,X_sampled[[block]][w,-j]) %*% normal_cd$coefficients, sd=sigma(normal_cd))
        normal_weight[[j]][[block]][w,] <- normal_weight[[j]][[block]][w,]/sum(normal_weight[[j]][[block]][w,])
      }
    }
    
  }
  
  check_list = list()
  check_list$W_count = c()
  check_list$X_count = c()
  check_list$inter_count = c()
  ## Start MCMC
  ## We use the Bayesian backfitting method.
  for (k in 2:nMCMC){
    ## Update the constant term "c"
    R_int_a <- 0 ## The sum of the interactions not updated for the current scan
    R_int_b <- 0 ## The sum of the interactions updated for the current scan
    for (j in Mod_ind) {R_int_a <- R_int_a + mu_hat_intact[[j]][k-1,]} ## All interactions are not updated.
    R <- Y - mu_hat_X[k-1,] - mu_hat_W[k-1,] - R_int_a - R_int_b ## Calculate the residuals w.o. the constant term.
    mu_hat_c[k,] <- mu_hat_c[k-1,] + mean(R) ## Update "c" by taking the mean of residuals.
    
    ## Update g(W)
    R <- Y - mu_hat_c[k] - mu_hat_X[k-1,] - R_int_a - R_int_b ## Residuals w.o. g(W).
    if (k==2) { ## Make a tree when we start MCMC.
      opts_soft_W <- SoftBart::Opts(num_burn = 0, num_save = 1, num_thin = 1, update_sigma = TRUE,
                                    update_sigma_mu = FALSE) ## Options for the tree
      hypers_W <- SoftBart::Hypers(W,Y)
      forest_W <- SoftBart::MakeForest(hypers_W,opts_soft_W)
    }
    mu_hat_W[k,] <- forest_W$do_gibbs(W,R,W,1) ## Update g(W) by do_gibbs(explanatory, response, evaluated at, # of scans)
    mu_hat_W_w0 <- forest_W$do_predict(matrix(W_star,nrow=1)) ## Calculate g(W*)
    
    shared_sigma <- forest_W$get_sigma() ## f, g, h will share the same residual variance.
    
    ## Update f(X)
    R <- Y - mu_hat_c[k] - mu_hat_W[k,] - R_int_a - R_int_b
    if (k==2) {
      opts_soft_X <- SoftBart::Opts(num_burn = 0, num_save = 1, num_thin = 1, update_sigma = TRUE,
                                    update_sigma_mu = FALSE) ## Options for the tree
      hypers_X <- SoftBart::Hypers(X,Y, num_tree = 50, sigma_hat = shared_sigma) ## f will use the residual variance from g.
      forest_X <- SoftBart::MakeForest(hypers_X,opts_soft_X)
    }else{
      forest_X$set_sigma(shared_sigma) ## f will use the residual variance from g.
    }
    mu_hat_X[k,] <- forest_X$do_gibbs(X,R,X,1) ## Update f(X)
    mu_hat_X_x0 <- forest_X$do_predict(matrix(X_star,nrow=1)) ## Calculate f(X*)
    
    check_list$W_count = forest_W$get_counts()
    check_list$X_count = forest_X$get_counts()
    
    ## Update h(X,W) = sum{h_j(X_j,W)}
    if (k==2) {
      opts_mod <- ModBart::Opts(num_burn = 0, num_save = 1, num_thin = 1, update_sigma = FALSE)
      hypers <- ModBart::Hypers(W,Y, sigma_hat = shared_sigma)
      
      hypers$length_scale <- 4 / pi
      mean_ell_sq <- hypers$length_scale^2
      hypers$shape_length_scale <- 1
      hypers$rate_length_scale <- 1 / (mean_ell_sq)
      
      forest <- list()
      for (j in Mod_ind) {forest[[j]] <- ModBart::MakeForest(hypers = hypers, opts = opts_mod)}
    }else {
      for (j in Mod_ind) {forest[[j]]$set_sigma(shared_sigma)}
    }
    rd_order <- sample(Mod_ind,length(Mod_ind)) ## Random order for the update
    
    mu_hat_intact_grid_x0 <- list()
    mu_hat_intact_grid_w0 <- list()
    mu_hat_intact_grid_00 <- list()
    
    temp_count = c()
    for (j in rd_order){
      R_int_a <- R_int_a - mu_hat_intact[[j]][k-1,]
      R <- Y - mu_hat_c[k] - mu_hat_X[k,] - mu_hat_W[k,] - R_int_a - R_int_b
      mu_hat_intact[[j]][k,] <- forest[[j]]$do_gibbs(W, R, X[,j], W, X[,j], 1)
      R_int_b <- R_int_b + mu_hat_intact[[j]][k,]
      
      mu_hat_intact_grid_x0[[j]] <- forest[[j]]$predict(obs_grid_x0[,1:q],obs_grid_x0[,q+j]) ## Calculate h_j(X*_j, W)
      mu_hat_intact_grid_w0[[j]] <- forest[[j]]$predict(obs_grid_w0[,1:q],obs_grid_w0[,q+j]) ## Calculate h_j(X_j, W*)
      mu_hat_intact_grid_00[[j]] <- forest[[j]]$predict(obs_grid_00[,1:q],obs_grid_00[,q+j]) ## Calculate h_j(X*_j, W*)
      temp_count = cbind(temp_count, forest[[j]]$get_counts())
    }
    check_list$inter_count = temp_count
    if(k %%100==0){
      print(check_list)
    }
    
    
    ## For the iterations that will not be discarded, calculate CATE E(Y(W1)-Y(W0)|X), the variable importance measures (VIM), psi_j, for each covariate.
    ## psi_j = (Heterogeneity without X_j)/(Total heterogeneity)
    ## where Total heterogeneity = E_W[Var_X{E(Y(W)-Y(W0)|X)}], the average (over W) variability (over X) of CATE with fixed W0.
    ## We similarly define "Heterogeneity without X_j" but integrating out X_j.
    ## We calculate the VIM for each sub-matrix and take the average of VIMs.
    if (k %in% save_index){
      ## Calculate CATE E(Y(W1)-Y(W0)|X) and CATE_hj E(Y(W1)-Y(W0)|X_j)
      ## In our model E(Y(W1)-Y(W0)|X) = g(W1)-g(W0) + h(X,W1)-h(X,W0)
      ## Also calculate E(Y(W)-Y(W0)|X=EX) = g(W)-g(W0)
      
      h_W0 <-0
      h_W1 <-0
      for (j in Mod_ind) {
        temp1 <- forest[[j]]$predict(W1_rep, X[,j])
        temp0 <- forest[[j]]$predict(W0_rep, X[,j])
        h_W1 <- h_W1 + temp1
        h_W0 <- h_W0 + temp0
        CATE_h_hat[[j]][k,] <- (temp1 - temp0) -
          rep(forest[[j]]$predict(t(W1), X_star[j]) - forest[[j]]$predict(t(W0), X_star[j]), num_obs)
      }
      E_Y_W1[k,] <- mu_hat_c[k,] + forest_X$do_predict(X) +
        rep(forest_W$do_predict(t(W1)),num_obs) + h_W1
      E_Y_W0[k,] <- mu_hat_c[k,] + forest_X$do_predict(X) +
        rep(forest_W$do_predict(t(W0)),num_obs) + h_W0
      
      CATE_hat[k,] <- E_Y_W1[k,] - E_Y_W0[k,]
      cat("ATE = ", mean(CATE_hat[k,])* sd_Y, "E_Y_W0=", mean(E_Y_W0[k,])* sd_Y, "E_Y_W1=", mean(E_Y_W1[k,])* sd_Y, "\r")
      
      if (is.null(Xtest)==FALSE){
        CATE_test_hat_int <- 0
        for (j in 1:p) {
          CATE_test_hat_int <- CATE_test_hat_int + forest[[j]]$predict(W1_rep_test, Xtest[,j]) - forest[[j]]$predict(W0_rep_test, Xtest[,j])
        }
        CATE_test_hat[k,] <- c(forest_W$do_predict(t(W1))-forest_W$do_predict(t(W0))) + CATE_test_hat_int
      }
      
      Hetero_total_blk <- rep(NA, n_blocks)
      Hetero_wo_blk <- matrix(NA, nrow = p, ncol = n_blocks)
      for (block in 1:n_blocks){
        ## First, calculate the total heterogeneity.
        ## E(Y(W)-Y(W0)|X) = g(W)-g(W0) + h(X,W)-h(X,W0)
        ##                 = g_tilde + h_tilde.
        ## Since g and h (also f) are not yet identifiable, we make some adjustments during the calculation.
        N <- sum(block_index==block)
        g_tilde_int1 <- 0
        g_tilde_int2 <- 0
        for (j in Mod_ind) {
          g_tilde_int1 <- g_tilde_int1 + forest[[j]]$predict(W_sampled[[block]], X_star_sampled[[block]][,j])
          g_tilde_int2 <- g_tilde_int2 + forest[[j]]$predict(t(W0), X_star[j])
        }
        g_tilde <-
          (forest_W$do_predict(W_sampled[[block]]) + g_tilde_int1 ) -
          (rep(forest_W$do_predict(t(W0)) + g_tilde_int2, N) )
        h_tilde <- list()
        for (j in Mod_ind){
          h_tilde[[j]] <-
            (forest[[j]]$predict(W_sampled_grid[[block]], X_sampled_grid[[block]][,j]) - 
               rep(forest[[j]]$predict(W_sampled[[block]], X_star_sampled[[block]][,j]), each=N ) ) -
            (forest[[j]]$predict(W0_grid[[block]], X_sampled_grid[[block]][,j]) - 
               rep(forest[[j]]$predict(W0_sampled[[block]], X_star_sampled[[block]][,j]), each=N) )
        }
        
        ## Take the mean of variance of CATE=(g_tilde + h_tilde).
        Hetero_vec_int <- 0
        for (j in Mod_ind){ ## Calculate 
          Hetero_vec_int <- Hetero_vec_int + h_tilde[[j]]
        }
        Hetero_total_vec <- 
          rep(g_tilde, each=N) + Hetero_vec_int # c(tau(X1,W1),tau(X2,W1),...,tau(X1,W2),...,tau(XN,WN))
        Hetero_total_mat <- matrix(Hetero_total_vec, nrow = N)
        Hetero_total_blk[block] <- mean(apply(Hetero_total_mat, 2, var))
        
        ## Second, calculate the heterogeneity in E(Y(W)-Y(W0)|X) without X_j.
        ## We are going to calculate E{E(Y(W)-Y(W0)|X)|X_{-j}}
        ##    where the outer conditional expectation is w.r.t. the conditional distribution of X_j|X_{-j}
        ##    where X_{-j} denotes the (p-1) covariates without X_j.
        for (j in Mod_ind){
          ## Assume X_j|X_{-j}~N(.,.) where the mean and the variance are estimated from the linear regression model.
          ## Integrating out X_j only affects h_tilde[[j]].
          ## Hence, we decompose E{E(Y(W)-Y(W0)|X)|X_{-j}} into two parts.
          
          ## (I) the part that does not involve X_j. -> remains as it is.
          Hetero_wo_j_vec_int = Hetero_vec_int - h_tilde[[j]]
          Hetero_wo_j_vec <-
            rep(g_tilde, each=N) + Hetero_wo_j_vec_int
          Hetero_wo_j_mat_base <- matrix(Hetero_wo_j_vec, nrow = N)
          
          ## (II) the part that does involve X_j. -> integrated out.
          h_j_tilde_mat <- matrix(h_tilde[[j]], nrow = N)
          E_h_j_tilde_mat <- normal_weight[[j]][[block]] %*% h_j_tilde_mat # Integrating out the CATE using the conditional density.
          
          Hetero_wo_j_mat_ind <- Hetero_wo_j_mat_base + E_h_j_tilde_mat # n by n matrix
          Hetero_wo_blk[j,block] <- mean(apply(Hetero_wo_j_mat_ind, 2, var)) # Take the variance for each W value and then take the mean of n variances.
        }
      }
      Hetero_total[k] <- mean(Hetero_total_blk) # Take the mean of the total heterogeneity.
      for (j in Mod_ind){
        Hetero_wo[[j]][k] <- mean(Hetero_wo_blk[j,]) # Take the mean of the heterogeneity w.o X_j.
        psi_hat[[j]][k] <- mean(1-Hetero_wo_blk[j,]/Hetero_total_blk) # Take the ratio.
      }
    }
    
    
    #### shifting ####
    ## This step is needed to make all functions identifiable
    ##    in the sense that f(X*) = G(W*) = 0 and h(X*, . ) = h( . , W*) = 0.
    sum_mu_hat_intact_grid_00 <- 0
    sum_mu_hat_intact_grid_x0 <- 0
    sum_mu_hat_intact_grid_w0 <- 0
    
    for (j in Mod_ind){
      sum_mu_hat_intact_grid_00 <- sum_mu_hat_intact_grid_00 + mu_hat_intact_grid_00[[j]]
      sum_mu_hat_intact_grid_x0 <- sum_mu_hat_intact_grid_x0 + mu_hat_intact_grid_x0[[j]]
      sum_mu_hat_intact_grid_w0 <- sum_mu_hat_intact_grid_w0 + mu_hat_intact_grid_w0[[j]]
    }
    
    mu_hat_c[k,] <- mu_hat_c[k,] + c(mu_hat_X_x0) + c(mu_hat_W_w0) + sum_mu_hat_intact_grid_00
    mu_hat_X[k,] <- mu_hat_X[k,] - c(mu_hat_X_x0) + sum_mu_hat_intact_grid_w0 - sum_mu_hat_intact_grid_00
    mu_hat_W[k,] <- mu_hat_W[k,] - c(mu_hat_W_w0) + sum_mu_hat_intact_grid_x0 - sum_mu_hat_intact_grid_00
    for (j in Mod_ind){
      mu_hat_intact[[j]][k,] <- mu_hat_intact[[j]][k,] - mu_hat_intact_grid_x0[[j]] - mu_hat_intact_grid_w0[[j]] + mu_hat_intact_grid_00[[j]]
    }
    
    if(k%%100==0) cat(k, "/", nMCMC, "Done", "\r")
    
  }
  
  ## Burning and stepping
  mu_hat_c <- mu_hat_c[save_index,] * sd_Y
  mu_hat_X <- mu_hat_X[save_index,] * sd_Y
  mu_hat_W <- mu_hat_W[save_index,] * sd_Y
  X_quantiles <- list()
  for (j in Mod_ind){
    mu_hat_intact[[j]] <- mu_hat_intact[[j]][save_index,] * sd_Y
    
    X_quantiles[[j]] <- quantile(X[,j], seq(0,1,0.01), type=1, na.rm = TRUE)
    Index_quantile <- match(X_quantiles[[j]], X[,j])
    CATE_h_hat[[j]] <- (CATE_h_hat[[j]][save_index, Index_quantile])* sd_Y
  }
  
  CATE_hat <- CATE_hat[save_index,] * sd_Y
  E_Y_W0 <- E_Y_W0[save_index,] * sd_Y
  E_Y_W1 <- E_Y_W1[save_index,] * sd_Y
  if (is.null(Xtest)==FALSE){CATE_test_hat <- CATE_test_hat[seq(BurnIn+1,nMCMC,stepsize),] * sd_Y }
  Hetero_total <- Hetero_total[save_index]
  for (j in Mod_ind){
    Hetero_wo[[j]] <- Hetero_wo[[j]][save_index]
    psi_hat[[j]] <- psi_hat[[j]][save_index]
  }
  
  ## Calculate the predicted outcomes
  sum_mu_hat_intact <- 0
  for (j in Mod_ind){
    sum_mu_hat_intact <- sum_mu_hat_intact + mu_hat_intact[[j]]
  }
  Y.hat <- mu_hat_c + mu_hat_X + mu_hat_W + sum_mu_hat_intact
  
  ## Calculate the posterior means
  CATE_hat.mean <- apply(CATE_hat,2,mean)
  ATE_hat <- apply(CATE_hat,1,mean)
  if (is.null(Xtest)==FALSE){CATE_test_hat.mean <- apply(CATE_test_hat,2,mean)
  }else{CATE_test_hat.mean <- NA}
  Hetero_total.mean <- mean(Hetero_total)
  Y.hat.mean <- apply(Y.hat,2,mean)
  mu.hat.c.mean <- apply(mu_hat_c,2,mean)
  mu.hat.X.mean <- apply(mu_hat_X,2,mean)
  mu.hat.W.mean <- apply(mu_hat_W,2,mean)
  mu.hat.intact.mean <- list()
  Hetero_wo.mean <- list()
  psi_hat.mean <- list()
  CATE_h_hat.mean <- list()
  for (j in Mod_ind){
    mu.hat.intact.mean[[j]] <- apply(mu_hat_intact[[j]],2,mean,na.rm=T)
    Hetero_wo.mean[[j]] <- mean(Hetero_wo[[j]],na.rm=T)
    psi_hat.mean[[j]] <- mean(psi_hat[[j]])
    CATE_h_hat.mean[[j]] <- apply(CATE_h_hat[[j]],2,mean,na.rm=T)
  }
  
  end.time.MCMC <- Sys.time()
  (running.time.MCMC <- end.time.MCMC-start.time.MCMC)
  print(running.time.MCMC)
  
  return(list(
    # Y.hat.mean = Y.hat.mean,
    #    mu.hat.c.mean = mu.hat.c.mean,
    #    mu.hat.X.mean = mu.hat.X.mean,
    #    mu.hat.W.mean = mu.hat.W.mean,
    #    mu.hat.intact.mean = mu.hat.intact.mean,
    #    mu_hat_c = mu_hat_c,
    #    mu_hat_X = mu_hat_X,
    #    mu_hat_W = mu_hat_W,
    #    mu_hat_intact = mu_hat_intact,
    CATE_h_hat = CATE_h_hat,
    CATE_h_hat.mean = CATE_h_hat.mean,
       CATE_hat.mean = CATE_hat.mean,
    #    CATE_test_hat.mean = CATE_test_hat.mean,
    Hetero_total.mean = Hetero_total.mean,
    Hetero_wo.mean = Hetero_wo.mean,
    psi_hat.mean = psi_hat.mean,
    ATE_hat = ATE_hat,
    #    CATE_hat = CATE_hat,
    #    CATE_test_hat = CATE_test_hat,
    Hetero_total = Hetero_total,
    #    Hetero_wo = Hetero_wo,
    psi_hat = psi_hat,
    sd_Y=sd_Y,
    res_sd = shared_sigma,
    X_quantiles_orig = X_quantiles_orig,
    W_quantiles_orig = W_quantiles_orig,
    X_names = X_names,
    W_names = W_names,
    Mod_names = Mod_names,
    # E_Y_W0 = E_Y_W0,
    # E_Y_W1 = E_Y_W1
    W1 = W1,
    W0 = W0
  ))
  
}
