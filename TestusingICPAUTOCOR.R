

library(glmnet)
library(dplyr)
library(mboost)

normalize_matrix_by_group <- function(matrix, groups) {
  matrix <- as.matrix(matrix)
  unique_groups <- unique(groups)
  centered_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  
  for (group in unique_groups) {
    group_indices <- which(groups == group)
    group_data <- matrix[group_indices, , drop = FALSE]
    
    # Center the group data along each column (subtract the mean)
    centered_group_data <- sweep(group_data, 2, colMeans(group_data), "-")
    
    # Assign the centered group data back to the corresponding rows in the result matrix
    centered_matrix[group_indices, ] <- centered_group_data
  }
  
  return(centered_matrix)
}

getblanketboosting <- function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
  p <- ncol(X)
  
  if(p <=maxNoVariables)
  {
    usevar <- 1:p
  }else
  {
    usevar <- selLmBoost(X = cbind(X,Y), k = p + 1, output = FALSE, pars = list(atLeastThatMuchSelected = 0.00, atMostThatManyNeighbors = maxNoVariables))
    usevar <- which(usevar)
    
  }
  
  testsets <- list()
  if(length(usevar)>0)
  {
    for (ic in ((1:2^length(usevar))-1))
    {
      testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
    }
  }
  testsets <- unique(testsets)
  le <- sapply(testsets,length)
  testsets <- testsets[ keep <- which(le>0 & le <= maxNoVariablesSimult) ]
  testsets <- testsets[order(le[keep])]
  return(testsets)
}

selLmBoost <- function(X,k,pars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),output = FALSE)
{
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
  if(output)
  {
    cat("Performing variable selection for variable", k, ": \n")
  }
  result <- list()
  p <- dim(as.matrix(X))
  if(p[2] > 1)
  {
    selVec <- rep(FALSE, p[2])
    a <- X[,-k]
    #names(a) <- LETTERS[1:(dim(X)[2]-1)]
    b <- X[,k]
    #names(b) <- "b"
    modfitLm <- train_LMboost(X[,-k],X[,k],pars)
    cc <- unique(modfitLm$model$xselect())
    if(output)
    {
      cat("The following variables \n")
      show(cc)
    }
    nstep <- length(modfitLm$model$xselect())
    howOftenSelected <- rep(NA,length(cc))
    for(i in 1:length(cc))
    {
      howOftenSelected[i] <- sum(modfitLm$model$xselect() == cc[i])/nstep
    }
    if(output)
    {
      cat("... have been selected that many times: \n")
      show(howOftenSelected)
    }
    howOftenSelectedSorted <- sort(howOftenSelected, decreasing = TRUE)
    if( sum(howOftenSelected>pars$atLeastThatMuchSelected) > pars$atMostThatManyNeighbors)
    {
      cc <- cc[howOftenSelected>howOftenSelectedSorted[pars$atMostThatManyNeighbors + 1]]
    } else
    {
      cc <- cc[howOftenSelected>pars$atLeastThatMuchSelected]
    }
    if(output)
    {
      cat("We finally choose as possible parents: \n")
      show(cc)
      cat("\n")
    }
    tmp <- rep(FALSE,p[2]-1)
    tmp[cc] <- TRUE
    selVec[-k] <- tmp
  } else
  {
    selVec <- list()
  }
  return(selVec)
}


train_LMboost <- function(X,y,pars = list()) #
{
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
  
  y <- y - rep( mean(y), length(y))
  
  ## begin old version
  #dat <- data.frame(cbind(yy,X))
  #gb <- glmboost(yy ~ .,data=dat)
  # EXPLANATION: surprisingly, it turned out that this cannot be applied to large p (private discussion with T. Hothorn in Sep 2013)
  yy <- as.vector(y)
  options(warn=-1)
  gb <- glmboost(X,yy, center = TRUE)
  options(warn=1)
  ## end old version
  
  ## begin new version
  #dat <- as.data.frame(X)
  #bl <- lapply(dat, bols)
  #gb <- mboost_fit(bl, y)
  ## end new version
  
  result <- list()
  result$Yfit <- gb$fitted()
  result$residuals <- gb$resid()
  result$model <- gb
  return(result)
}

getblanketlasso <-
  function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
    p <- ncol(X)
    
    if(p <=maxNoVariables){
      usevar <- 1:p
    }else{
      lobs <- coef(glmnet(X,Y,family= if(!is.factor(Y)) "gaussian" else "binomial"))[-1,]
      nnz <- apply(lobs!=0,2,sum)
      nnzsel <- 0
      usevarcandidate <- numeric(0)
      usevar <- numeric(0)
      while(length(usevarcandidate)<maxNoVariables & nnzsel<max(nnz)){
        nnzsel <- nnzsel+1
        sel <- which( nnz==nnzsel)
        if(length(sel)>0) usevarcandidate <- sort(unique(c(usevarcandidate, which(apply(lobs[,sel,drop=FALSE]!=0,1,any)))))
        if(length(usevarcandidate)<=maxNoVariables) usevar <- usevarcandidate
      }
      # the following applies if several variables enter at once
    }
    testsets <- list()
    if(length(usevar)>0){
      for (ic in ((1:2^length(usevar))-1)){
        testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
      }
    }
    testsets <- unique(testsets)
    le <- sapply(testsets,length)
    testsets <- testsets[ keep <- which(le>0 & le <= maxNoVariablesSimult) ]
    testsets <- testsets[order(le[keep])]
    return(testsets)
  }
pvalfunc <- function( x, y, test="normal"){
  
  if(is.function( test)){
    pval <- test(x,y)
  }else{
    if(test %in% c("ranks","ks")){
      if( test=="ranks"){
        pval <- 2*min( t.test(x,y)$p.value, var.test(x,y)$p.value)
      }else{
        pval <- ks.test(x,y)$p.value
      }
    }else{
      pval <- 2*min( t.test(x,y)$p.value, var.test(x,y)$p.value)
    }
  }
  return(pval)
}


getpval <-
  function(Y,X,IN,test="exact",maxNoObs=200){
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
    
    if (ncol(X) != 1) {
      for (ki in 1:length(IN)) {
        ARLinm <- glm(Y[IN[[ki]]] ~ X[IN[[ki]], 2, drop = FALSE] - 1, family = "gaussian", control = glm.control(maxit = 10, epsilon = 10^(-6)))
        ARPred <- predict(ARLinm)
        Y[IN[[ki]]] <- Y[IN[[ki]]] - ARPred
      }
    }
    X <- as.matrix(X[,-2])
    linm <- glm( Y ~ X -1, family= if(is.factor(Y)) "binomial" else "gaussian", control=glm.control(maxit=10,epsilon=10^(-6)))
    coefficients <- coef(linm)[-1]
    coefficientsvar <- summary(linm)$coefficients[,2][-1]
    pred <- predict(linm, type="response")
    
    if(is.factor(Y)){
      resid <- ((as.numeric(Y)-1)-pred)
      pvalvec <- length(IN)
      for (ki in 1:length(IN)){
        pvalvec[ki] <- getpvalClassif(resid[IN[[ki]]],resid[-IN[[ki]]],test=test)
      }
      pval <- min(pvalvec)*(length(IN)-1)
      
    }else{
      useexact <- FALSE
      if( !is.function(test)){
        if( test=="exact") useexact <- TRUE
      }
      if(!useexact){
        
        K <- length(IN)
        
        resid <- residuals(linm)
        pvalvec <- numeric(length(IN))
        for (ki in 1:length(IN)){
          pvalvec[ki] <- pvalfunc( resid[IN[[ki]]], resid[-IN[[ki]]],test=test)
          if(!is.function(test)){
            if(test=="correlation"){
              nonzerosd <- which(apply(X,2,sd)>0)
              if(length(nonzerosd)>0){
                corpval <- numeric(length(nonzerosd))
                for (k in 1:length(corpval)) corpval[k] <- cor.test(X[IN[[ki]],nonzerosd[k]], resid[IN[[ki]]])$p.value
                pvalvec[ki] <- 2*min(pvalvec[ki], length(corpval)*min(corpval))
              }
            }
          }
        }
        pval <- min(pvalvec)*(length(IN)-1)
        
      }else{
        
        K <- length(IN)
        
        if(!is.matrix(X)) X <- as.matrix(X)
        n <- nrow(X)
        
        pvalvec <- numeric(length(IN))
        for (ki in 1:length(IN)){
          nk <- length(IN[[ki]])
          nko <- n-nk
          p <- ncol(X)
          linm <- lm.fit( X[-IN[[ki]], ,drop=FALSE] , Y[-IN[[ki]]]) ## fit a model on all other data
          pred <- as.numeric(X[IN[[ki]], ,drop=FALSE] %*% coefficients(linm))
          diff <- Y[IN[[ki]]] - pred
          
          selobs <-  if( nk>maxNoObs)  sample( IN[[ki]], maxNoObs) else IN[[ki]]
          if( nk>maxNoObs) diff <- diff[ IN[[ki]] %in% selobs]
          nk <- length(selobs)
          COV <- diag(length(diff)) + X[ selobs,] %*% solve(t(X[-IN[[ki]],])%*%X[-IN[[ki]],], t(X[ selobs,]))
          
          stat <- (t(diff)%*% solve(COV, diff)) / (nk * var(residuals(linm))*nko/(nko-p))
          pval <- 1-pf(stat, nk, n-nk - ncol(X))
          
          pvalvec[ki] <- pval
        }
        pval <- min(pvalvec) * (length(IN))
      }
      
    }
    pval <- min(1,pval)
    
    return(list(pval=pval,coefficients=coefficients,coefficientsvar=coefficientsvar))
  }





getblanketall <- function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
  p <- ncol(X)

  usevar <- 1:p
  
  testsets <- list()
  if(length(usevar)>0){
    for (ic in ((1:2^length(usevar))-1)){
      testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
    }
  }
  testsets <- unique(testsets)
  le <- sapply(testsets,length)
  testsets <- testsets[order(le)]
  return(testsets)
}


ICPAutoCor <- function (X, Y,ylagged, ExpInd, alpha = 0.01, test = "exact", selection = c("lasso", 
                                                                            "all", "stability", "boosting")[if (ncol(X) <= 8) 2 else 4], 
                 maxNoVariables = 8, maxNoVariablesSimult = 8, maxNoObs = 200, 
                 showAcceptedSets = FALSE, showCompletion = TRUE, stopIfEmpty = FALSE, gof= max(0.01,alpha)) 
{
###
# Adapted from the code of the package InvariantCausalPrediction which was published by Nicolai Meinshausen und GPL-3
#
###
  if (is.vector(X) & is.numeric(X)) 
    X <- matrix(X, ncol = 1)
  if (!is.matrix(X) & !is.data.frame(X)) 
    stop("'X' must be a matrix or data frame")
  if (!is.vector(Y) & !is.factor(Y)) 
    stop("'Y' must be a vector or factor")
  if (is.function(test)) {
    pval <- test((1:10) + 0.5, 1:10)
    if (!is.numeric(pval)) 
      stop("function 'test' has to return a numeric value")
    if (length(pval) > 1) 
      stop("function 'test' needs to return a scalar (the p-value of the null hypothesis test that 'x' and 'z' are sampled from the same distribution")
    if (pval < 0 | pval > 1) 
      stop("the p-value of function 'test' needs to be in [0,1]")
  }
  if (!is.list(ExpInd)) {
    if (length(ExpInd) != length(Y)) 
      stop("if `ExpInd' is a vector, it needs to have the same length as `Y'")
    uni <- unique(ExpInd)
    if (length(uni) == 1) 
      stop(paste("there is just one environment ('ExpInd'=", 
                 uni[1], " for all observations) and the method needs at least two distinct environments", 
                 sep = ""))
    if (min(table(ExpInd)) <= 2) {
      cat("\n out put of 'table(ExpInd)':\n ")
      print(table(ExpInd))
      stop("one environment has just one or two observations (as supplied by 'ExpInd'); there need to be at least 3 (and ideally dozens) of observations in each environment; the output of 'table(ExpInd)' is given below to show the number of observations in each unique environment as supplied by 'ExpInd'")
    }
    K <- length(uni)
    ExpIndNEW <- list()
    for (uc in 1:K) {
      ExpIndNEW[[uc]] <- which(ExpInd == uni[uc])
      attr(ExpIndNEW[[uc]], "value") <- uni[uc]
    }
    ExpInd <- ExpIndNEW
    rm(ExpIndNEW)
  } else {
    ran <- range(unlist(ExpInd))
    if (ran[1] < 1) 
      stop(paste("if `ExpInd' is a list with indicies of observations, \n minimal entry has to be at least 1 but is", 
                 ran[1]))
    if (ran[2] > length(Y)) 
      stop(paste("if `ExpInd' is a list with indicies of observations, \n maximal entry has to be at most equal \n to the length", 
                 length(Y), "of the observations but is", ran[2]))
    if (min(sapply(ExpInd, length) <= 2)) 
      stop("one environment has just one or two observations (as supplied by 'ExpInd'); there need to be at least 3 (and ideally dozens) of observations in each environment")
  }
  getblanket <- getblanketall
  if (selection == "lasso") {
    getblanket <- getblanketlasso
  }
  if (selection == "stability") {
    getblanket <- getblanketstability
  }
  if (selection == "boosting") {
    getblanket <- getblanketboosting
  }
  if (is.data.frame(X)) {
    if (any(sapply(X, class) == "factor")) {
      Z <- X
      X <- matrix(nrow = nrow(Z), ncol = sum(sapply(X, 
                                                    function(x) if (is.factor(x)) length(levels(x)) else 1)))
      cc <- 0
      colX <- character(0)
      for (k in 1:ncol(Z)) {
        if (is.numeric(Z[, k])) {
          cc <- cc + 1
          X[, cc] <- Z[, k]
          colX <- c(colX, colnames(Z)[k])
        }
        else {
          nf <- length(lev <- levels(Z[, k]))
          for (nfc in 1:(nf - 1)) {
            cc <- cc + 1
            X[, cc] <- as.numeric(Z[, k] == lev[nfc])
            colX <- c(colX, paste(colnames(Z)[k], "_", 
                                  lev[nfc], sep = ""))
          }
        }
      }
      X <- X[, 1:cc]
      colnames(X) <- colX
      X <- as.matrix(X)
    }
  }
  if (length(ucol <- unique(colnames(X))) < min(3, ncol(X))) 
    colnames(X) <- paste("Variable", 1:ncol(X), sep = "_")
  if (length(unique(Y)) == 2 & !is.factor(Y)) {
    warning("\n Y only has 2 unique values -- using classification")
    Y <- as.factor(Y)
  }
  K <- length(ExpInd)
  n <- nrow(X)
  p <- ncol(X)
  X <- cbind(ylagged,X)
  X <- cbind(rep(1, nrow(X)), X)
  ConfInt <- matrix(NA, nrow = 2, ncol = p)
  Coeff <- list()
  CoeffVar <- list()
  for (k in 1:p) {
    Coeff[[k]] <- numeric(0)
    CoeffVar[[k]] <- numeric(0)
  }
  pvall <- numeric(0)
  Pall <- numeric()
  cont <- TRUE
  pvalempty <- getpval(Y, X[, c(1), drop = FALSE], ExpInd, test = test, 
                       maxNoObs = maxNoObs)$pval
  if (pvalempty > alpha) {
    ConfInt <- matrix(0, nrow = 2, ncol = p)
    pvall <- c(pvall, pvalempty)
    if (showAcceptedSets) 
      cat(paste("\n accepted empty set"))
    if (stopIfEmpty) 
      cont <- FALSE
  }
  testsets <- if (any(!is.null(c(maxNoVariables, maxNoVariablesSimult)))){
    getblanket(X[, -c(1,2), drop = FALSE], Y, maxNoVariables = maxNoVariables, 
               maxNoVariablesSimult = maxNoVariablesSimult)
  }else{
    getblanket(X[, -c(1,2), drop = FALSE], Y)
  }
  len <- sapply(testsets, length)
  lcsingle <- sum(len == 1)
  usedvariables <- unique(unlist(testsets))
  lc <- 0
  printoutat <- ifelse(length(testsets) > 0, 2^(1:ceiling(log2(length(testsets)))), 
                       NA)
  if(cont & stopIfEmpty) intersection <- NULL
  acceptedSets <- list()
  while (cont && lc < length(testsets)) {
    if (showCompletion) {
      if (lc %in% printoutat) {
        cat(paste("\n *** ", round(100 * lc/length(testsets)), 
                  "% complete: tested ", lc, " of ", length(testsets), 
                  " sets of variables ", sep = ""))
      }
    }
    lc <- lc + 1
    usevariab <- testsets[[lc]]
    notusevariab <- (1:p)[-testsets[[lc]]]
    tmp <- getpval(Y, X[, c(1, 2,2 + usevariab), drop = FALSE], 
                   ExpInd, test = test, maxNoObs = maxNoObs)
    pval <- tmp$pval
    Pall <- c(Pall, pval)
    if (pval > alpha) {
      if(stopIfEmpty){
        if(is.null(intersection)){
          intersection <- usevariab
        }else{
          intersection <- intersect(intersection,usevariab)
        }
        if( length(intersect)==0) cont <- FALSE
      }
      acceptedSets[[length(acceptedSets)+1]] <- usevariab
      
      if (showAcceptedSets) 
        cat(paste("\n accepted set of variables ", paste(usevariab, 
                                                         collapse = ","), sep = ""))
      ConfInt[1, usevariab] <- pmax(ConfInt[1, usevariab, 
                                            drop = FALSE], tmp$coefficients + qnorm(1 - alpha/4) * 
                                      tmp$coefficientsvar, na.rm = TRUE)
      ConfInt[2, usevariab] <- pmin(ConfInt[2, usevariab, 
                                            drop = FALSE], tmp$coefficients - qnorm(1 - alpha/4) * 
                                      tmp$coefficientsvar, na.rm = TRUE)
      for (kc in usevariab) {
        Coeff[[kc]] <- c(Coeff[[kc]], tmp$coefficients[which(usevariab == 
                                                               kc)])
        CoeffVar[[kc]] <- c(CoeffVar[[kc]], tmp$coefficientsvar[which(usevariab == 
                                                                        kc)])
      }
      if (length(notusevariab) >= 1) {
        ConfInt[1, notusevariab] <- pmax(ConfInt[1, notusevariab, 
                                                 drop = FALSE], 0, na.rm = TRUE)
        ConfInt[2, notusevariab] <- pmin(ConfInt[2, notusevariab, 
                                                 drop = FALSE], 0, na.rm = TRUE)
      }
      pvall <- c(pvall, pval)
    }
  }
  colnames(ConfInt) <- colnames(X[, -c(1,2), drop = FALSE])
  if (is.null(colnames(ConfInt))) 
    colnames(ConfInt) <- paste("Variable", 1:ncol(ConfInt))
  sig <- apply(sign(ConfInt[2:1, , drop = FALSE]), 2, function(x) prod(x))
  sigo <- sign(ConfInt[1, ])
  maximin <- sigo * apply(abs(ConfInt[2:1, , drop = FALSE]), 
                          2, min) * (sig >= 0)
  pvalues <- rep(1, p)
  modelReject <- (!any(c(pvalempty,Pall) > gof))
  if (!modelReject) {
    for (k in usedvariables) {
      sel <- which(sapply(testsets, function(x, z) z %in% 
                            x, k))
      if(length(sel) > 0){
        pvalues[k] <- if(length(Pall[-sel])>0) max(Pall[-sel], pvalempty) else 1
      }
    }
  }else{
    ConfInt <-  NULL
    maximin <- NULL
  }
  
  retobj <- list(ConfInt = ConfInt, maximinCoefficients = maximin, 
                 alpha = alpha, colnames = colnames(ConfInt), factor = is.factor(Y), 
                 dimX = dim(X), Coeff = Coeff, CoeffVar = CoeffVar, modelReject = modelReject, acceptedSets = acceptedSets,
                 usedvariables = usedvariables, pvalues = pvalues, stopIfEmpty=stopIfEmpty, noEnv = length(ExpInd), gof=gof, bestModel=max(c(pvalempty,Pall)))
  class(retobj) <- "InvariantCausalPrediction"
  return(retobj)
}




gen_mat <- function(data, lag, timepoint,name){
  current_data <- data[data$Name == name, ]
  current_data <- current_data[current_data$Time ==timepoint, ]
  
  lagged_data <- data[data$Time ==timepoint- lag, ]
  ylagged <- lagged_data[lagged_data$Name == name,][['Value']]
  lagged_data <- lagged_data[lagged_data$Name != name, ]
  wide_df <- pivot_wider(lagged_data, id_cols = Country, names_from = Name, values_from = Value)
  wide_df <- as.matrix(wide_df[ , -1])
  y <- current_data[['Value']]
  
  ExpInd <- current_data[['ContCode']]
  return(list(wide_df = wide_df, y = y, ExpInd = ExpInd,ylagged=ylagged))
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
  stacked_ylagged <- do.call(c,lapply(result_list, function(x) x[[4]]))
  return(list(wide_df = stacked_wide_df, y = stacked_y, ExpInd = stacked_ExpInd,ylagged = stacked_ylagged))
}


test_all_icp <- function(data) {
  names <- unique(data$Name)
  all_results <- list()
  timepoints <-unique(data$Time)
  lag <-1
  for (effect in names) {
    icpdat <- stack_gen_mat(data,lag=lag,timepoints= (lag + min(as.numeric(timepoints))):max(as.numeric(timepoints)),effect)
    X <- as.matrix(icpdat$wide_df)
    Y <- as.numeric(icpdat$y)
    ylagged <- as.numeric(icpdat$ylagged)
    ExpInd <- as.numeric(icpdat$ExpInd)
    
    icpres <- ICPAutoCor(X, Y, ylagged, ExpInd, alpha = 0.05, selection = "lasso")
    
    all_results[[effect]] <- icpres
  }
  
  return(all_results)
}
#res <- test_all_icp(df)
res<- ICPAutoCor(data_simul$X,as.numeric(data_simul$Y),as.numeric(data_simul$Y_lagged),data_simul$ExpInd,alpha=0.05)