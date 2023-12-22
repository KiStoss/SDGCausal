library(tidyr)
library(grpreg)
library(glmnet)
library(lmtest)
library(dplyr)
library(vars)
library(bruceR)
library(countrycode)
library(xts)
library(progress)
N=110
gen_mat <- function(data, lag, timepoint,name){
  current_data <- data[data$Name == name, ]
  current_data <- current_data[current_data$Time ==timepoint, ]
  lagged_data <- data[data$Name != name, ]
  lagged_data <- lagged_data[lagged_data$Time ==timepoint- lag, ]
  wide_df <- pivot_wider(lagged_data, id_cols = Country, names_from = Name, values_from = Value)
  wide_df <- as.matrix(wide_df[ , -1])
  y <- current_data[['Value']]
  ExpInd <- current_data[['numeric_continent']]
  return(list(wide_df = wide_df, y = y, ExpInd = ExpInd))
}

stack_gen_mat <- function(data, lag, timepoints,name) {
  result_list <- list()
  
  for (tp in timepoints) {
    gen_result <- gen_mat(data, lag, tp,name)
    result_list[[as.character(tp)]] <- gen_result
  }
  
  stacked_wide_df <- do.call(rbind, lapply(result_list, function(x) x[[1]]))
  stacked_y <- do.call(c, lapply(result_list, function(x) x[[2]]))
  stacked_ExpInd <- do.call(c, lapply(result_list, function(x) x[[3]]))
  
  return(list(wide_df = stacked_wide_df, y = stacked_y, ExpInd = stacked_ExpInd))
}


generate_multivariate_ts <- function(Nobs, N, coef,lag_order = 1, country = "CountryA") {
  # Initialize data frame to store the reshaped time series data
  ts_data <- data.frame()
  continent <- as.character(countrycode(country, "country.name", "region"))
  continent_mapping <- c("Sub-Saharan Africa" = 1, "Antarctica" = 2, "East Asia & Pacific" = 3, "Europe & Central Asia" = 4, "North America" = 5, "Oceania" = 6, "Latin America & Caribbean" = 7)
  num_cont <- continent_mapping[continent]
  random_ar_coeff <-coef[[num_cont]]$ar
  random_ma_coeff <-coef[[num_cont]]$ma
  # Generate AR processes with random coefficients and reshape the data
  for (i in 1:N) {
    rarcoeff<-random_ar_coeff[i]
    rmacoeff<-random_ma_coeff[i]
    
    
    rarcoeff <- 0.5 * sin(rarcoeff * pi) 
    ts_values <- arima.sim(model = list(ar = c(rarcoeff, 0)), n = Nobs+1)
    ts_df <-data.frame(
      Name = paste0("TS", i),
      Time = 0:Nobs,
      Value = as.numeric(ts_values),
      Country = rep(country, Nobs+1)
    )
    ts_data <- rbind(ts_data, ts_df)
  }
  lol<-lag(ts_data$Value[ts_data$Name == "TS1"], 1)
  ts_data$Value[ts_data$Name == "TS101"] <- ts_data$Value[ts_data$Name == "TS101"] + 0.5 * lag(ts_data$Value[ts_data$Name == "TS1"], 1)
  
  ts_data$Value[ts_data$Name == "TS102"] <- ts_data$Value[ts_data$Name == "TS102"] - 0.5 * lag(ts_data$Value[ts_data$Name == "TS12"], 1) + 0.5 * lag(ts_data$Value[ts_data$Name == "TS22"], 1)
  
  ts_data$Value[ts_data$Name == "TS103"] <- ts_data$Value[ts_data$Name == "TS103"] + 0.5 * lag(ts_data$Value[ts_data$Name=="TS33"], 1)
  
  ts_data$Value[ts_data$Name == "TS104"] <- ts_data$Value[ts_data$Name == "TS104"] + 0.5 * lag(ts_data$Value[ts_data$Name == "TS33"], 1) + 0.5 * lag(ts_data$Value[ts_data$Name == "TS44"], 1)
  
  ts_data$Value[ts_data$Name == "TS105"] <- ts_data$Value[ts_data$Name == "TS105"] + 0.5 * lag(ts_data$Value[ts_data$Name == "TS5"], 1) - 0.5 * lag(ts_data$Value[ts_data$Name == "TS15"], 1)
  
  ts_data$Value[ts_data$Name == "TS106"] <- ts_data$Value[ts_data$Name == "TS106"] - 0.5 * lag(ts_data$Value[ts_data$Name == "TS20"], 1) + 0.5 * lag(ts_data$Value[ts_data$Name == "TS30"], 1)
  
  ts_data$Value[ts_data$Name == "TS107"] <- ts_data$Value[ts_data$Name == "TS107"] + 0.5 * lag(ts_data$Value[ts_data$Name == "TS40"], 1) - 0.5 * lag(ts_data$Value[ts_data$Name == "TS50"], 1)
  
  ts_data$Value[ts_data$Name == "TS108"] <- ts_data$Value[ts_data$Name == "TS108"] + 0.5 * lag(ts_data$Value[ts_data$Name == "TS60"], 1) + 0.5 * lag(ts_data$Value[ts_data$Name == "TS70"], 1)
  
  ts_data$Value[ts_data$Name == "TS109"] <- ts_data$Value[ts_data$Name == "TS109"] - 0.5 * lag(ts_data$Value[ts_data$Name == "TS80"], 1) + 0.5 * lag(ts_data$Value[ts_data$Name == "TS90"], 1)
  
  ts_data$Value[ts_data$Name == "TS110"] <- ts_data$Value[ts_data$Name == "TS110"] + 0.5 * lag(ts_data$Value[ts_data$Name == "TS10"], 1)
  ts_data <- ts_data[ts_data$Time != 0, ]
  
  return(ts_data)
}
generate_multivariate_ts_all_countries <- function(Nobs, N, coef,lag_order = 1, countries = c("CountryA", "CountryB", "CountryC")) {
  
  # Initialize data frame to store the reshaped time series data
  all_ts_data <- data.frame()
  
  # Generate data for each country
  for (country in countries) {
    ts_data <- generate_multivariate_ts(Nobs, N,coef, lag_order, country)
    all_ts_data <- rbind(all_ts_data, ts_data)
  }
  
  return(all_ts_data)
}

generate_random_coefficients <- function(N, ar_coeff_range = c(0.2, 0.9)) {
  random_ar_coefficients <- runif(N, min = ar_coeff_range[1], max = ar_coeff_range[2])
  random_ma_coefficients <- runif(N, min = -0.5, max = 0.5)
  return(list(ar = random_ar_coefficients, ma = random_ma_coefficients))
}

coefficients_list <- generate_random_coefficients(7 * N)
coefficient_map <- lapply(1:7, function(i) {
  start_index <- (i - 1) * N + 1
  end_index <- i * N
  list(ar = coefficients_list$ar[start_index:end_index], ma = coefficients_list$ma[start_index:end_index])
})




granger_high_dim_univariate <- function(data, lag = 1, N) {
  results <-list()
  
  indices <- unique(data$Name)
  
  for (name in indices) {
    BigMat <- stack_gen_mat(data, lag, c((lag + 1):(N - lag)), name)
    X <- BigMat$wide_df
    y <- BigMat$y
    
    lasso_model <- cv.glmnet(X, y, alpha = 1)
    selected_features <- coef(lasso_model, s = "lambda.min", intercept=FALSE)
    nonzero_features <- rownames(selected_features)[selected_features[, 1] != 0]
    intercept_index <- grep("^\\(Intercept\\)$", nonzero_features)
    
    nonzero_features <- nonzero_features[-intercept_index]
    
    for (pred in nonzero_features) {
      x = data[data$Name==pred,]
      y = data[data$Name==name,]
      gran <- granger_test(y$Value~x$Value, lags=1, test.reverse=FALSE, data=data)
      results <- c(results,data.frame(gran))
    }
  }
  
  return(results)
}

granger_high_dim_multi <- function(data, lag = 1, N) {
  results <- data.table()
  
  indices <- unique(data$Name)
  
  # Load the progress package
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(indices))
  }
  
  for (name in indices) {
    if (exists("pb")) pb$tick()  # Increment the progress bar
    
    BigMat <- stack_gen_mat(data, lag, c((lag + 1):(N - lag)), name)
    X <- BigMat$wide_df
    y <- BigMat$y
    
    lasso_model <- cv.glmnet(X, y, alpha = 1)
    selected_features <- coef(lasso_model, s = "lambda.min", intercept=FALSE)
    nonzero_features <- rownames(selected_features)[selected_features[, 1] != 0]
    intercept_index <- grep("^\\(Intercept\\)$", nonzero_features)
    nonzero_features <- nonzero_features[-intercept_index]
    
    # Check if nonzero_features is empty
    if (length(nonzero_features) == 0) {
      # Add "No causes" to result_row and go to the next iteration
      result_row <- data.table(
        Name = name,
        Granger = "No Causes found in preselection model fitted",
        pval = 1
      )
      
      #results <- rbindlist(list(results, result_row), use.names = TRUE, fill = TRUE)
      next
    }
    
    features <- data %>%
      filter(Name %in% c(nonzero_features, name))
    
    time_series_objects <- list()
    for (feature_name in c(nonzero_features,name)) {
      subset_df <- features[features$Name == feature_name, ]
      time_series_objects[[feature_name]] <- xts(subset_df$Value, order.by = as.POSIXct(subset_df$Time))
    }
    
    # Combine time series into a matrix with named columns
    varser_matrix <- do.call(cbind, lapply(names(time_series_objects), function(ts_name) {
      colnames(time_series_objects[[ts_name]]) <- ts_name
      coredata(time_series_objects[[ts_name]])
    }))
    
    # Try fitting the VAR model and handle errors
    tryCatch({
      tsVar <- VAR(varser_matrix, p = 1)
      cause_variable <- names(time_series_objects)[1] 
      causagran <- causality(tsVar, cause = nonzero_features)$Granger
    }, error = function(e) {
      result_row <- data.table(
        Name = name,
        Granger = "No VAR model fitted",
        pval = 1
      )
      #causagran <- "VAR model could not be computed"
    })
    
    # Create a data.table with results
    result_row <- data.table(
      Name = name,
      Granger = causagran[["method"]],
      pval = causagran[["p.value"]]
    )
    
    results <- rbindlist(list(results, result_row), use.names = TRUE, fill = TRUE)
  }
  
  if (exists("pb")) pb$terminate()  # Close the progress bar
  
  return(results)
}
nobs = 100
data <- generate_multivariate_ts_all_countries(nobs,110,coefficient_map,1,countries = "Switzerland")
res <-granger_high_dim_multi(data,N=nobs)
library(HDGCvar)
HDGC_mat <- data  %>% dplyr::select(-Country)%>%
  pivot_wider(names_from = Name, values_from = Value)%>% dplyr::select(-Time)
HDGC_mat <- as.matrix(HDGC_mat)
interest_variables=list("GCto"="TS106","GCfrom"="TS1")
#HDGC_VAR_all_I0(HDGC_mat,p=1,bound=0.2*nrow(data))
#res<-HDGC_VAR(GCpair = interest_variables, data=as.matrix(HDGC_mat),p=2,d=1,bound = 0.1*nrow(HDGC_mat))
library(igraph)
#Plot_GC_all(res, Stat_type="FS_cor",alpha=0.01, multip_corr=list(F),directed=T, layout=layout.circle, main="Network",edge.arrow.size=.2,vertex.size=5, vertex.color=c("lightblue"), vertex.frame.color="blue",vertex.label.size=2,vertex.label.color="black",vertex.label.cex=0.6, vertex.label.dist=1, edge.curved=0,cluster=list(T,5,"black",0.8,1,0)) 


