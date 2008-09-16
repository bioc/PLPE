'generate.null' <- function(x, design, q)
{

    n.pairs <- ncol(x)/2
    n.genes <- nrow(x)
    x.null <- matrix(NA, 1, 2*n.pairs)

    if(n.genes*q < 5) stop('q is too small.')
    if(q > 1) stop('q is too large.')
    
    i <- (1:n.pairs)*2
    x1 <- as.vector(x[,-i]) 
    x2 <- as.vector(x[,i])
    
    A <- (x1+x2)/2
    quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm=TRUE)
    n.quan <- length(quantile.A)-1
    
    for(j in 1:n.quan) {
    
        k <- which((A >= quantile.A[j]) & (A <= quantile.A[j+1]))
        x.vec <- cbind(x1[k], x2[k])
        rk1 <- rank(x.vec[,1])
        rk2 <- rank(x.vec[,2])
    
        M <- abs(rk1-rk2)
        median.M <- median(M)
        k0 <- which(M <= median.M)
    
        for(i in 1:round(length(k)/n.pairs)) {
            x0 <- x.vec[sample(k0, size=n.pairs),]
            x.null <- rbind(x.null, as.vector(t(x0)))
        }
    }
    x.null <- x.null[-1,]
    
    return(x.null)

}
