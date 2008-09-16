'am.trans.paired' <- function (y, design) 
{

    i <- sort.list(design[,2])
    y <- y[,i]; design <- design[i,]
    n <- length(unique(design[,2]))

    A <- c()
    M <- c()
    for(i in 1:n) {
        A <- c(A,(y[,(2*i-1)]+y[,(2*i)])/2)
        M <- c(M,(y[,(2*i-1)]-y[,(2*i)])  )

        A <- c(A,(y[,(2*i-1)]+y[,(2*i)])/2)
        M <- c(M,(y[,(2*i)]-y[,(2*i-1)])  )
    }

    return(cbind(A, M))
}
