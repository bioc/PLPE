'lpe.paired.default' <- function (x, design, data.type, q=0.01, probe.ID = NULL, estimator="median",  
                                    w=0.5, w.estimator="fixed", iseed=1234, ...) 
{
    require(Biobase)

    #Set a seed
    set.seed(iseed)
 
    #Sort by conditions 
    i <- sort.list(design[,1])
    x <- x[,i]; design <- design[i,]

    #Sort by pairs 
    i <- sort.list(design[,2])
    x <- x[,i]; design <- design[i,]

    #Counts
    n <- ncol(x)
    n.paired <- round(n/2)
    n.genes <- nrow(x)
    u.cond <- unique(design[,1])
    
    #Check data and design
    if(length(unique(design[,1])) != 2) stop('The number of conditions should be 2.')
    if(length(unique(design[,2]))*2 != ncol(x)) stop('Check pair numbers.')
    if(n.genes*q < 5) stop('q is too small.')
    if(q > 0.5) stop('q is too large.')
    
    #####################################################################
    #Obtain M and A
    i <- which(design[,1]==u.cond[1])
    x1 <- x[,i]; x2 <- x[,-i]
    M <- x1-x2
    A <- (x1+x2)/2 
    
    if(estimator=="mean") {
       stats=mean 
       est.M <- apply(M, 1, mean, na.rm=TRUE)
       est.A <- apply(A, 1, mean, na.rm=TRUE)
    }
    
    if(estimator=="median") {
       stats=median
       est.M <- apply(M, 1, median, na.rm=TRUE)
       est.A <- apply(A, 1, median, na.rm=TRUE)
    }
    
    if(estimator=="huber") {
       stats=huber
       require(MASS)
       
       est.M <- apply(M, 1, median, na.rm=TRUE)
       est.A <- apply(A, 1, median, na.rm=TRUE)
       
       M2 <- M; A2 <- A
       i <- which(ncol(M2)-apply(is.na(M2),1, sum) < 2)
       if(length(i) > 0) {
          M2 <- M2[-i,]; A2 <- A2[-i,]
       }
       
       est.M2 <- t(matrix(unlist(apply(M2, 1, huber)), 2, nrow(M2)))[,1]
       est.A2 <- t(matrix(unlist(apply(A2, 1, huber)), 2, nrow(A2)))[,1]
       
       if(length(i) > 0) {
          est.M[-i] <- est.M2 
          est.A[-i] <- est.A2
       }

       if(length(i) == 0) {
          est.M <- est.M2
          est.A <- est.A2
       }      
    }

    #####################################################################
    #Paired t-test
    t.stat <- var.t <- p.value.t <- rep(NA, n.genes) 

    d <- M
    i <- which(ncol(d)-apply(is.na(d),1, sum) < 2)
    if(length(i) > 0) d <- d[-i,]
    i <- which(apply(d, 1, var)==0)
    if(length(i) > 0) d <- d[-i,] #CHECK

    paired.t <- apply(d, 1, t.test) 
    paired.t <- t(matrix(unlist(paired.t), 10, nrow(d)))

    t.stat2 <- as.numeric(paired.t[,1])
    var.t2  <- apply(d, 1, var, na.rm=TRUE) 
    p.value.t2 <- as.numeric(paired.t[,3]) #CHECK

    if(length(i) > 0) {
       t.stat[-i] <- t.stat2
       var.t[-i] <- var.t2
       p.value.t[-i] <- p.value.t2
    }
    if(length(i) ==0) {
       t.stat <- t.stat2
       var.t <- var.t2
       p.value.t <- p.value.t2
    }

    #####################################################################
    #LPEP test
    
    base.var <- base.error.paired(x, design,est.A, estimator, q, data.type)
    var.L <- base.var$var.genes  
    L.stat <- est.M/sqrt(var.L)
    p.value.L <- (1-pnorm(abs(L.stat)))*2

    #####################################################################
    #Weighted LPEP test    
    
    var.Lw <- (w * var.t + (1-w) * var.L)  
    Lw.stat <- est.M/sqrt(var.Lw)
    p.value.Lw <- (1-pnorm(abs(Lw.stat)))*2
    
    i <- which(is.na(var.Lw)==TRUE)
    if(length(i) >0) {
        var.Lw[i] <- var.L[i] 
        Lw.stat[i] <- L.stat[i]
        p.value.Lw[i] <- p.value.L[i]

    }

    #####################################################################
    #Estimate w
    
    if(w.estimator=="random") {
       div=100
       ww <- c(0, (1:div)/div)
       tmp <- rep(NA, (div+1))
       for(i in 1:(div+1)) {
           var.Lww <- (ww[i] * var.t + (1-ww[i]) * var.L)  
           Lww.stat <- est.M/sqrt(var.Lww)
           tmp[i] <- abs(lm(Lww.stat~A)$coef[2])
       }
       w <- (min(which(tmp==min(tmp)))-1)/div 

       p.value.Lw <- (1-pnorm(abs(Lw.stat)))*2
       i <- which(is.na(var.Lw)==TRUE)
       if(length(i) >0) {
           var.Lw[i] <- var.L[i] 
           Lw.stat[i] <- L.stat[i]
           p.value.Lw[i] <- p.value.L[i]
       }
    } 
        

    #####################################################################
    #Outputs
  
    data.out <- data.frame(A=est.A, M=est.M, 
                var.L=var.L,   L.stat =L.stat,  p.value.L =p.value.L, 
                var.t=var.t,   t.stat =t.stat,  p.value.t =p.value.t,
                var.Lw=var.Lw, Lw.stat=Lw.stat, p.value.Lw=p.value.Lw)
   
    #Add probe IDs                       
    if(length(probe.ID) == nrow(x)) row.names(data.out) <- probe.ID
    else row.names(data.out) <- 1:nrow(data.out)    

    # results
    res <- list(design=design, data.type=data.type, q=q, 
                estimator=estimator, w.estimator=w.estimator, w=w, test.out=data.out, iseed = iseed)
    class(res) <- "lpe.paired"
    res
}
