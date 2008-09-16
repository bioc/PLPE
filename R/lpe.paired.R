##########################################################################
#
#        
#    Local Pooled Error (LPE) Test for Paired High-throughput Data
#
#                    by
#
#          HyungJun Cho and Jae K. Lee
#
#
##########################################################################

.First.lib <- function(lib, pkg) { 
   cat("LPEP version 1.0.0 \n") 
   invisible()
   if(.Platform$OS.type=="windows" && require(Biobase) && interactive() 
   && .Platform$GUI=="Rgui") { addVigs2WinMenu("LPEP") }
   if( !require(methods) ) stop("we require methods for package LPEP")
   
}

'lpe.paired' <- function(x,...)
{
    UseMethod("lpe.paired")
}
