## return periods

kendallRP <- function(kendallFun=NULL, cl=c(.99,.999), mu=1, copula=NULL) {
  if(is.null(kendallFun) & is.null(copula)) stop("Either the kendall distribution function or the copula must be provided. Note that the calculation of the kendall distribution function from the copula is pretty time consuming. Saving them separately might be advantougous.")
  if(is.null(kendallFun)) kendallFun <- genEmpKenFun(copula)
  if(length(mu)>1 & length(cl) > 1) stop("Either the critial level (cl) or mu may be of length larger than 1!")
  return(mu/(1-kendallFun(cl)))
}   

criticalLevel <- function(kendallFun=NULL, KRP=c(100,1000), mu=1, copula=NULL) {
  if(is.null(kendallFun) & is.null(copula)) stop("Either the kendall distribution function or the copula must be provided. Note that the calculation of the kendall distribution function from the copula is pretty time consuming. Saving them separately might be advantougous.")
  if(is.null(kendallFun)) kendallFun <- genEmpKenFun(copula)
  if(length(mu)>1 & length(KRP) > 1) stop("Either the kendall return period or mu may be of length larger than 1!")
  invKenFun <- genInvKenFun(kendallFun)
  return(invKenFun(1-mu/KRP))
}
## next: calculating critical layer, sampling from the layer, selecting "typical" points