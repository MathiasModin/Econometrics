# This function provides an automated way of performing the augmented dickey fuller test.


# Main function -----------------------------------------------------------
require(urca)

auto.adf <- function(data, lags="AIC", maxlag=18, alpha="0.5"){
  
  # Pre-allocation of space
  type <- c("trend", "drift", "none") # defining all types
  biglist <- list()
  # Setting critical value corresponding to significance level
  # alpha is then the corresponding column in the distribution table
  if (alpha == "0.01"){
    alpha <- 2
  }else if (alpha == "0.1"){
    alpha <- 4
  }else{
    alpha <- 3
  }
  
  # Converting data to matrix and storing information
  data <- as.matrix(data)
  names <- colnames(data)
  urr <- matrix(NA, ncol(data), 2)
  urr[,1]<- names
  
  
  # looping over the different types of tests
  for (j in 1:ncol(data)){
    results <- list(trend=NA, drift=NA, none=NA) # 

    for (i in 1:3){
    model <- ur.df(data[,j], type=type[i], lags=maxlag, selectlags=lags)
    results[[i]] <- cbind(t(model@teststat), model@cval)      
    }
    
    # Decision algorithm for unit root
    if (results$trend[1,1]< results$trend[1,alpha]){          # if tau3 < crit
      ur <- "tau3 rejected => No unit root"
    }else if (results$trend[3,1]< results$trend[3, alpha] &&  # if Phi3 < crit &
              results$drift[1,1] > results$drift[1, alpha] && # if tau2 > crit &
              results$drift[2,1] < results$drift[2, alpha] && # if Phi1 < crit &
              results$none[1,1] > results$none[1, alpha]){    # if tau1 > crit
      ur <- "has a unit root"
    }else if (results$drift[1,1]< results$drift[1, alpha]){   # tau2 < crit
      ur <- "tau2 rejected => No unit root"
    }else if (results$none[1,1] < results$none[1, alpha]){    # tau1 < crit
      ur <- "tau1 rejected => No unit root"    
    }
    jthname <- names[j]
    biglist[[jthname]] <- results
    urr[j,2] <- ur
  }
  z <- "urr"
  biglist[[z]] <- urr
  return(biglist)
}

adft <- function(data, lags=12, selectlags="AIC"){
    types <- c('trend', 'drift', 'none')
    
    for (i in 1:3){
      print(paste0('################',' ADF test with ', types[i],'################'))
      model <- ur.df(data, type = types[i], lags = lags, selectlags = selectlags)
      print(cbind(t(model@teststat),model@cval))
      print('####################################################') 
    }
}
