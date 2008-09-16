###################################################################################
#
# FDR computation
#
###################################################################################
lpe.paired.fdr.default <- function(x, obj, n.iter=5, lambda=0.9, ...) {

    #Save parameters
    design <- obj$design
    q <- obj$q
    w <- obj$w
    estimator <- obj$estimator 
    w.estimator <- obj$w.estimator 

    data.type <- obj$data.type

    #Sort by conditions 
    i <- sort.list(design[,1])
    x <- x[,i]; design <- design[i,]

    #Sort by pairs 
    i <- sort.list(design[,2])
    x <- x[,i]; design <- design[i,]

    n.gene <- nrow(x)
    ii <- rep(1, n.gene)

    L0 <- c(); Lw0 <- c() 
    for(i in 1:n.iter) {
        x0 <- generate.null(x=x, design=design, q=q)
        tmp <-  lpe.paired(x=x0, design=design, q=q, data.type=data.type, 
                           estimator=estimator, w=w, w.estimator=w.estimator)
        L0  <- c(L0,  abs(tmp$test.out$L.stat)) 
        Lw0 <- c(Lw0, abs(tmp$test.out$Lw.stat)) 
    }

    ##################################################################### 
    #Using L stat  
    L <- obj$test.out$L.stat 
    R   <- rep(NA, nrow(x))
    for(i in 1:nrow(x)) R[i] <- sum(ii[(abs(L) >= abs(L[i]))]) #more efficient?
     
    jj <- rep(1, length(L0))
    m.lambda <- quantile(L0, probs=lambda)
    pi0 <- sum(ii[(abs(L) <= m.lambda)])/(sum(jj[(L0 <= m.lambda)])/n.iter)

    R0  <- rep(NA, n.gene)
    for(i in 1:n.gene) R0[i]  <- sum(jj[(L0 >= abs(L[i]))])/n.iter 
  
    FDR.L <- R0/R
    FDR.L[FDR.L >1] <- 1
    FDR.L <- FDR.L*pi0

    ##################################################################### 
    #Using Lw stat    
    Lw <- obj$test.out$Lw.stat 
    Rw   <- rep(NA, nrow(x))
    for(i in 1:nrow(x)) Rw[i] <- sum(ii[(abs(Lw) >= abs(Lw[i]))]) #more efficient?

    jj <- rep(1, length(Lw0))
    m.lambda <- quantile(Lw0, probs=lambda)
    pi0 <- sum(ii[(abs(Lw) <= m.lambda)])/(sum(jj[(Lw0 <= m.lambda)])/n.iter)

    Rw0  <- rep(NA, n.gene)
    for(i in 1:n.gene) Rw0[i]  <- sum(jj[(Lw0 >= abs(Lw[i]))])/n.iter 

    FDR.Lw <- Rw0/Rw
    FDR.Lw[FDR.Lw >1] <- 1
    FDR.Lw <- FDR.Lw*pi0

    ##################################################################### 
    FDR <- data.frame(L, FDR.L, Lw, FDR.Lw)
    res <- list(FDR=FDR, pi0=pi0, design=design, data.type=data.type, estimator=estimator, w=w, 
                n.iter = n.iter, lambda = lambda, w.estimator=w.estimator)
    class(res) <- "lpe.paired.fdr"    
    res
}


###################################################################################END


