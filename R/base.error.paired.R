'base.error.paired' <- function (x, design,est.A, estimator, q, data.type) 
{
    require(LPE)
    #AM transformation
    AM <- na.omit(am.trans.paired(x, design))
    A <- AM[, 1]; M <- AM[, 2]
     
    #Compute quantiles of A
    quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm=TRUE)
    if (sum(A == min(A)) > (q * length(A))) {
        tmpA <- A[!(A == min(A))]
        quantile.A <- c(min(A), quantile(tmpA, probs=seq(q,1,q), na.rm=TRUE))
    }
    quan.n <- length(quantile.A)-1
    
    #Estimate signals and noises     
    var.M <- rep(NA, length = quan.n)
    est.AA <- rep(NA, length = quan.n)
    for (i in 1:quan.n) {
        bin <- which(A > quantile.A[i] & A <= quantile.A[i+1])
        if(length(bin) < 2) bin <- which((A >= quantile.A[i]) & (A <= quantile.A[i+1]))
        n.i <- length(!is.na(M[bin]))
        mult.factor <- (n.i-0.5)/(n.i-1)
        if(estimator =="mean") {
           var.M[i] <- mult.factor * var(M[bin], na.rm=TRUE)
           est.AA[i] <- mean(A[bin], na.rm=TRUE)
        }
        if(estimator =="median") {
           var.M[i] <- mult.factor * var(M[bin], na.rm=TRUE)
           est.AA[i] <- median(A[bin], na.rm=TRUE)
        }

        if(estimator =="huber") {
           var.M[i] <- mult.factor * (huber(M[bin])$s)^2
           est.AA[i] <- huber(A[bin])$mu
        }
    }

    #Adjustments (do not anything if RMA or MS) #CHECK 
    #if(data.type=="mas5" | data.type=="mas4" |data.type=="dchip") 
    #   var.M[1:max(which(var.M == max(var.M)))] <- max(var.M)
    #if(data.type=="cdna") { #lower 5% of A Upper???
    #   i <- which(medianAs <= quantile(medianAs, probs=0.05, na.rm=TRUE)) 
    #   max.var <- max(var.M[i])
    #   k <- max(which(var.M[i] == max.var)) 
    #   var.M[1:k] <- max.var
    #}

    #Interpolation and estimation
    base.var <- cbind(A = est.AA, var.M = var.M)
    df <- min(10, round(nrow(base.var)/2))
    sm.spline <- smooth.spline(base.var[, 1], base.var[, 2], df=df)
    var.genes <- fixbounds.predict.smooth.spline(sm.spline, est.A)$y #CHECK

    #Adjustment for negative variances
    if(length(which(var.genes <=0)) > 0) {
       i <-  which(var.genes <=0) 
       var.genes[i] <- min(var.genes[-i]) 
    }       

    return(list(A=est.A, var.genes=var.genes))
}
