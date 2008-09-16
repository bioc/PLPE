'summary.lpe.paired' <- function(object, ...)
{
    cat("\n Local Pooled Error Test for Paired Data \n")

    cat("\n data type :", format(object$data.type))
    cat("\n specification for the estimator :", format(object$estimator))
    cat("\n quantile for intercals of intensities :", object$q)
    cat("\n weight parameter :", object$w)
    cat("\n estimate the weight :", object$w.estimator)
    cat("\n seed numner :", object$iseed)
    cat("\n \n design matrix \n")
    print(object$design)
    cat("\n \n matrix for test results \n")
    print(head(object$test.out))
    cat(ifelse(dim(object$test.out)[1] > 6, "...", ""), "\n")
    invisible(object)
}
