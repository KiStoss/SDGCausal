library(bnlearn)

conditional_independence_test <- function(X, Y, ylagged,effect) {
  num_predictors <- ncol(X)
  p_values <- numeric(num_predictors)
  names <- unique(colnames(X))
  result_df <- data.frame(Effect = character(num_predictors), Predictor = character(num_predictors), P_Value = numeric(num_predictors))
  for (i in 1:num_predictors) {
    data <- cbind(X,Y,ylagged)
    result <- ci.test(x = "Y", y = names[i], z = c(names[-i],"ylagged"),data=data.frame(data),,test="smc-zf")
    
    p_values[i] <- result$p.value

  }
  result_df$Effect <- effect
  result_df$Predictor <- names
  result_df$P_Value <- p_values
  
  
  return(result_df)
}

test_all_condip <- function(data) {
  #This function takes a dataframe containing colums Name, Time,Value,Country where Name specifies the name of the time series it then reshapes it to the correct form and applies a conditional independence test
  names <- unique(data$Name)
  all_results <- list()
  timepoints <- unique(data$Time)
  lag <- 1
  result_df <- NULL  
  
  for (effect in names) {
    icpdat <- stack_gen_mat(data, lag = lag, timepoints = (lag + min(as.numeric(timepoints))):max(as.numeric(timepoints)),effect)
    X <- as.matrix(icpdat$wide_df)
    Y <- as.numeric(icpdat$y)
    ylagged <- as.numeric(icpdat$ylagged)
    ExpInd <- as.numeric(icpdat$ExpInd)
    
    incondres <- conditional_independence_test(X, Y, ylagged,effect)

    all_results[[effect]] <- incondres
    
    if (is.null(result_df)) {
      result_df <- incondres
    } else {
      result_df <- rbind(result_df, incondres)
    }
  }
  
  return(result_df)
}

res <- test_all_condip(df)
#data_simul<-data.frame(data_simul)
#X <- as.matrix(subset(data_simul,select=-c(Y,Y_lagged,ExpInd)))
#Y <- as.matrix(subset(data_simul, select=Y))
#ylagged <-as.matrix(subset(data_simul, select=Y_lagged))
#res <- conditional_independence_test(X,Y,ylagged,"Y")
