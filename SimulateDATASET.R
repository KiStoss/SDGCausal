generate_data <- function(N = 10) {
#This code generates a data set containing N+1 time series, the time series also have a ExpInd for the case one wants to apply it to invariant causal prediction
  set.seed(123)
  
  phi_list<- lapply(1:N, function(x) runif(4, 0, 0.5))
  phi_y <- runif(4, 0.5, 1)
  phi_z <- runif(4, 0, 1)
  Nobs <- 5000 #Number of observations in each time series 
  max_rows <- 0
  X_stacked <- NULL
  Y_stacked <- NULL
  Y_lagged_stacked <- NULL
  ExpInd <- integer(0)
  burn_in<-2
  for (i in 1:4) {
    X_list <- vector("list", N)
    
    for (j in 1:N) {
      phi_x_i <- phi_list[[j]][i]
      rand_val = runif(1,-100,100)
      X <- as.numeric(rep(0, times = Nobs+burn_in))
      #X = numeric(Nobs+1)
      X <-X+arima.sim(model = list(ar = phi_x_i), n = Nobs + burn_in)
      X_list[[j]] <- X
    }
    
    cause_ind1 <- 1
    cause_ind2 <- 4
    cause_ind3 <- 100
    effect_ind<-50
    rand_val = runif(1,-10,10)
    Y <- as.numeric(rep(0, times = Nobs+burn_in))
    #Y <- numeric(Nobs+burn_in)
    Y <-Y+arima.sim(model = list(ar = phi_y[i]), n = Nobs + burn_in)
    
    lagged_cause_ind1 <- c(NA, X_list[[cause_ind1]][-(length(X_list[[cause_ind1]])-1)])
    
    lagged_cause_ind2 <- c(NA, X_list[[cause_ind2]][-(length(X_list[[cause_ind2]])-1)])
    
    lagged_cause_ind3 <- c(NA, X_list[[cause_ind3]][-(length(X_list[[cause_ind3]])-1)])
    
    Y[(burn_in):length(Y)] <- Y[(burn_in):length(Y)] + 
      0.7 * lagged_cause_ind1[(burn_in):length(Y)] + 
      0.6 * lagged_cause_ind2[(burn_in):length(Y)] +
      0.75* lagged_cause_ind3[(burn_in):length(Y)]
    X_list[[effect_ind]][(1+burn_in):length(Y)]<-X_list[[effect_ind]][(1+burn_in):length(Y)]+0.7*Y[burn_in:(length(Y)-1)]
    #Z <- numeric(Nobs + 1)
    #Z <- arima.sim(model = list(ar = phi_z[i]), n = Nobs + 1)
    #Z[2:length(Z)] <- 0.5 * Y[1:(Nobs)]
    
    X_matrix <- do.call(cbind, X_list)
    colnames(X_matrix) <- paste0("X_", seq_along(X_list))
    
    Y_matrix <- matrix(Y[(1+burn_in):length(Y)], ncol = 1)
    Y_lagged_matrix <- matrix(Y[burn_in:(length(Y) - 1)], ncol = 1)
    
    max_rows <- max(max_rows, nrow(X_matrix))
    
    X_stacked <- rbind(X_stacked, X_matrix[(1+burn_in):length(Y), ])
    Y_stacked <- rbind(Y_stacked, Y_matrix)
    Y_lagged_stacked <- rbind(Y_lagged_stacked, Y_lagged_matrix)
    
    ExpInd <- c(ExpInd, rep(i, Nobs))
  }
  
  return(list(X = X_stacked, Y = Y_stacked, Y_lagged = Y_lagged_stacked,ExpInd = ExpInd))
  
}

data_simul <-generate_data(N=100)