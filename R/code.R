#### regularized projection depth ----
#' Regularized Projection Depth
#'
#' Computation of the regularized projection depth for multivariate and 
#' functional data.
#' 
#' @param mu A single vector of dimension \code{d}, or a numerical matrix
#' of dimension \code{m}-times-\code{d} of \code{m} vectors in dimension 
#' \code{d}. Each row of the matrix corresponds to one point in which
#' the regularized projection outlyingness/depth should be calclulated.
#' 
#' @param X A numerical matrix with the dataset of dimension 
#' \code{n}-times-\code{d}, where \code{n} is the size of the dataset and 
#' \code{d} is its dimension.
#' 
#' @param beta A non-negative number specifying the amount of regularization 
#' to be applied. If \code{quant=TRUE}, must be a number in the interval 
#' \code{[0,1]}. If \code{quant=FALSE}, can be any non-negative number. Choice
#' \code{beta=0} in any case means that no regularization is applied. By 
#' default \code{beta=0}.
#' 
#' @param quant An indicator of whether the parameter \code{beta} stands for
#' regularization in terms of quantiles (\code{TRUE}, default), or in terms
#' of the nominal value of the MAD (\code{FALSE}).
#' 
#' @param ndir A number of randomly chosen directions in the unit sphere to 
#' approximate the projection depth. By default \code{ndir=1e4}.
#' 
#' @param ndir2 A number of randomly chosen directions in the unit sphere to 
#' find the threshold for MAD based on \code{beta} if \code{quant=TRUE}. This
#' parameter is ignored if \code{quant=FALSE}. 
#' 
#' @return A numerical vector \code{res} of approximate (regularized) projection 
#' outlyingness values. Points with lower values correspond to deeper points in
#' the sample. Convert to approximate (regularized) depth by transforming
#' \code{1/(1+res)}.
#' 
#' @seealso \link[ddalpha:depth.projection]{depth.projection}
#'
#' @examples
#' n = 500
#' X = cbind(rnorm(n),5*rnorm(n))
#' 
#' B = 1.5
#' evl = 1e2
#' xg = B*seq(min(X[,1]),max(X[,1]),length=evl)
#' yg = B*seq(min(X[,2]),max(X[,2]),length=evl)
#' mu = as.matrix(expand.grid(xg,yg))
#' m = nrow(mu)
#' 
#' # (non-regularized) projection outlyingness and depth
#' O = RPDepth(mu, X)
#' D = matrix(1/(1+O),ncol=evl)
#' 
#' contour(xg, yg, D, col="orange", lwd=2)
#' points(X, cex=.15, pch=16)
#' 
#' # regularized projection outlyingness and depth
#' beta = .5
#' O2 = RPDepth(mu, X, beta)
#' D2 = matrix(1/(1+O2),ncol=evl)
#' 
#' contour(xg, yg, D2, col="navy", lwd=2, add=TRUE)

RPDepth = function(mu, X, beta = 0, quant = TRUE, ndir = 1e4, ndir2 = 1e3){
  
  n = nrow(X)
  d = ncol(X)
  if((!is.matrix(mu))&(d==1)) mu = matrix(mu,ncol=1)
  if((!is.matrix(mu))&(d>1))  mu = matrix(mu,nrow=1)
  m = nrow(mu)
  if(ncol(mu)!=d) stop("The dimensions of X and mu do not match.")
  
  res = RPD_outl_C(mu, X, beta, quant, n, d, m, ndir, ndir2)
  
  return(res)
}