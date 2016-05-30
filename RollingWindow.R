# Timer -------------------------------------------------------------------

timer <- function(seconds){
    m <- seconds %/% 60
    s <- round(seconds %% 60,3)
    h <- m %/% 60
    m <- m %% 60
    return(paste0(h,"h ",m,"m ",s,"s"))
}



# Arima -------------------------------------------------------------------

r.arima <- function(x, n.ahead=12, size=0.65, order=c(1,0,0), type=c("Arima", "auto.arima")){
    s <- proc.time()
    require(forecast)
    win.size <- floor(length(x)*size)
    n.win <- length(x) - win.size
    fc <- matrix(NA, n.ahead+1, n.win)
    act <- matrix(NA, n.ahead+1, n.win-n.ahead)
    hRMSE <- matrix(NA, 1, n.ahead)
    rRMSE <- ts(matrix(NA, n.win-n.ahead, 1), start=start(x)+c(0, win.size), freq=frequency(x))
    colnames(fc) <- paste0("win=", seq(1:n.win))
    rownames(fc) <- c("origin", seq(1:(n.ahead)))
    rownames(fc)[2:(n.ahead+1)] <- paste0("h=", seq(1:n.ahead))
    rownames(hRMSE) <- "RMSE"
    colnames(hRMSE) <- paste0("h=", seq(1:n.ahead))
    rownames(rRMSE) <- paste(1:(n.win-n.ahead))
    colnames(rRMSE) <- "RMSE"
    
    plot(x, lwd=2, col="darkblue", ylab="")
    legend("bottom",
           col = c("darkblue", "darkred"),
           c(paste0("Observed ", deparse(substitute(x))), paste0(n.ahead, "-step ahead forecasts")),
           bty = "n",lwd=c(2,1))
    
    for (i in 1:n.win){
        StartDate <- start(x) + c(0,i-1)
        EndDate <- StartDate + c(0, win.size-1)
        sample <- window(x, start=StartDate, end=EndDate)
        if (type=="auto.arima"){
            model <- auto.arima(sample, trace=F)
        }else{
            model <- Arima(sample, order=order) 
        }
        
        fc[2:(n.ahead+1),i]<- forecast(model, h=n.ahead)$mean
        fc[1,i] <-sample[length(sample)]
        lines(ts(fc[,i], start=EndDate, freq=frequency(x)), col="darkred")
        if (i <= n.win-n.ahead){
            act[,i] <- window(x, start=EndDate, end=EndDate + c(0,n.ahead))
            rRMSE[i,1] <- accuracy(fc[-1,i], act[,i])[2]
        }
    }
    
    for (i in 1:n.ahead){
        hRMSE[i]<- accuracy(fc[i+1,], act[i+1,])[2]
    }
    
    print(timer(proc.time()[3]-s[3]))
    return(list(forecasts=fc, act.vals=act, hRMSE=hRMSE, rRMSE=rRMSE))
}


# VAR ---------------------------------------------------------------------

r.var <- function(x, n.ahead=12, size=0.65, var=1, p="bic"){
    s <- proc.time()
    require(vars)
    require(forecast)
    
    if(class(p)=="character"){var.select=1}else{var.select=0}
    if(p=="aic"){p=1}
    if(p=="hq"){p=2}
    if(p=="bic"){p=3}
    if(p=="fpe"){p=4}
    
    win.size <- floor(nrow(x)*size)
    n.win <- nrow(x) - win.size
    fc <- matrix(NA, n.ahead+1, n.win)
    act <- matrix(NA, n.ahead+1, n.win-n.ahead)
    hRMSE <- matrix(NA, 1, n.ahead)
    rRMSE <- ts(matrix(NA, n.win-n.ahead, 1), start=start(x)+c(0, win.size), freq=frequency(x))
    colnames(fc) <- paste0("win=", seq(1:n.win))
    rownames(fc) <- c("origin", seq(1:(n.ahead)))
    rownames(fc)[2:(n.ahead+1)] <- paste0("h=", seq(1:n.ahead))
    rownames(hRMSE) <- "RMSE"
    colnames(hRMSE) <- paste0("h=", seq(1:n.ahead))
    rownames(rRMSE) <- paste(1:(n.win-n.ahead))
    colnames(rRMSE) <- "RMSE"
    
    plot(x[,var], lwd=2, col="darkblue", ylab="")
    legend("bottom",
           col = c("darkblue", "darkred"),
           c(paste0("Observed ", colnames(x)[var]), paste0(n.ahead, "-step ahead forecasts")),
           bty = "n",lwd=c(2,1))
    
    for (i in 1:n.win){
        StartDate <- start(x) + c(0,i-1)
        EndDate <- StartDate + c(0, win.size-1)
        sample <- window(x, start=StartDate, end=EndDate)
        if (var.select==1){
            model<- VAR(sample,p=VARselect(sample,lag.max=12)$selection[p])
        }else{
            model<- VAR(sample,p=p)
        }
        
        fc[2:(n.ahead+1),i]<- forecast(model, h=n.ahead)$mean[[var]]
        fc[1,i] <-sample[nrow(sample),var]
        lines(ts(fc[,i], start=EndDate, freq=frequency(x)), col="darkred")
        if (i <= n.win-n.ahead){
            act[,i] <- window(x[,var], start=EndDate, end=EndDate + c(0,n.ahead))
            rRMSE[i,1] <- accuracy(fc[-1,i], act[,i])[2]
        }
    }
    
    for (i in 1:n.ahead){
        hRMSE[i]<- accuracy(fc[i+1,], act[i+1,])[2]
    }
    print(timer(proc.time()[3]-s[3]))
    return(list(forecasts=fc, act.vals=act, hRMSE=hRMSE, rRMSE=rRMSE))
}
