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
#' @param alpha A non-negative (possibly \code{Inf}) number specifying the 
#' amount of regularization to be applied to the operator part. Always an
#' absolute value.
#' 
#' @param beta A non-negative number specifying the amount of regularization 
#' to be applied in the MAD part. If \code{quant=TRUE}, must be a number 
#' in the interval \code{[0,1]}. If \code{quant=FALSE}, can be any non-negative 
#' number. Choice \code{beta=0} in any case means that no regularization is 
#' applied to the MAD part. By default \code{beta=0}.
#' 
#' @param quant An indicator of whether the parameter \code{beta} stands for
#' regularization in terms of quantiles (\code{TRUE}, default), or in terms
#' of the nominal value of the MAD (\code{FALSE}).
#' 
#' @param operator The operator specifying the penalization for the operator
#' part. Must be a list of two components: (i) \code{values}, a numerical vector
#' of length \code{K} with positive eigenvalues of the operator Gamma, and (ii)
#' \code{vectors}, a numerical matrix of dimension \code{d}-times-\code{K} whose
#' columns specify the eigenvectors of Gamma. By default \code{NULL}, which
#' means that no regularization is applied to the operator part. 
#' 
#' @param ndir A number of randomly chosen directions in the unit sphere to 
#' approximate the projection depth. By default \code{ndir=1e4}.
#' 
#' @param ndir2 A number of randomly chosen directions in the unit sphere to 
#' find the threshold for MAD based on \code{beta} if \code{quant=TRUE}. This
#' parameter is ignored if \code{quant=FALSE}. 
#' 
#' @param Jmax The maximum number of unit directions to be sampled and checked
#' for whether they satisfy the conditions imposed on the MAD part and the 
#' operator part. If \code{Jmax} directions are considered but \code{ndir} of 
#' them do not satisfy both conditions, gives an error. By default, 
#' \code{Jmax=1e6}.
#' 
#' @param echo A boolean value that should be used only for diagnostic purposes. If 
#' \code{echo=TRUE}, as a part of the output also diagnostic messages about the 
#' process of searching for regularized directions is given. By default
#' \code{echo=FALSE}.
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
#' 
#' # regularized depth: regularization with the operator
#' # (A) custom operator
#' operator = list()
#' operator$values = c(0.01,100)
#' operator$vectors = cbind(c(1,0),c(1,1))
#' alpha = 1.5
#' D = RPDepth(mu, X, alpha = alpha, beta = 0, operator = operator)
#' D1 = matrix(1/(1+D),ncol=evl)
#' 
#' # (B) variance operator of X
#' Sigma = var(X)
#' operator = eigen(Sigma)
#' alpha = 10
#' D = RPDepth(mu, X, alpha = alpha, beta = 0, operator = operator)
#' D2 = matrix(1/(1+D),ncol=evl)
#' 
#' contour(xg, yg, D1, col="orange", lwd=2)
#' points(X, cex=.15, pch=16)
#' contour(xg, yg, D2, add=TRUE, col=2, lwd=.5)

RPDepth = function(mu, X, alpha = 0, beta = 0, quant = TRUE, operator = NULL, 
                   ndir = 1e4, ndir2 = 1e3, Jmax = 1e6, 
                   echo = FALSE){
  
  n = nrow(X)
  d = ncol(X)
  if((!is.matrix(mu))&(d==1)) mu = matrix(mu,ncol=1)
  if((!is.matrix(mu))&(d>1))  mu = matrix(mu,nrow=1)
  m = nrow(mu)
  if(ncol(mu)!=d) stop("The dimensions of X and mu do not match.")
  
  if(!is.null(operator)){
    gammas = operator$values
    basis = operator$vectors
    ngammas = length(gammas)
  } else {
    alpha = 0
    ngammas = 1
    gammas = 1
    basis = matrix(rep(0,d),ncol=1)
  }
  if(any(gammas<=0)) 
    stop("The eigenvalues in operator$vectors must be all positive.")
  gammas = 1/sqrt(gammas)
  
  res = RPD_outl_C(mu, X, alpha, beta, quant, 
                   gammas, basis, n, d, m, ngammas, 
                   ndir, ndir2, Jmax,
                   echo)
  
  return(res)
}

#### MAD ----
#' Median Absolute Deviation
#'
#' Standard median absolute deviation (MAD) of a numerical vector. Compared to
#' function \link[stats:mad]{mad}, no additional adjustment factor for 
#' consistency with normal distribution is used.
#' 
#' @param x A numerical vector whose MAD we compute.
#' 
#' @return A single numerical value, median of the absolute distances from the 
#' median of \code{x}.
#' 
#' @seealso \link[stats:mad]{mad}
#'
#' @examples
#' n = 10
#' X = rnorm(n)
#' 
#' MAD_R(X)
#' mad(X)/1.4826 # mad is using an adjustment factor 1.4826

MAD_R = function(x) median(abs(x - median(x)))

#### projection pursuit PCA ----
#' Projection Pursuit Principal Component Analysis
#'
#' The PCA using a projection pursuit algorithm from Croux and Gazen (2005).
#' 
#' @param X A numerical matrix with the dataset of dimension 
#' \code{n}-times-\code{d}, where \code{n} is the size of the dataset and 
#' \code{d} is its dimension.
#' 
#' @param q Number of components to obtain. By default set to \code{d}, the 
#' dimension of the dataset.
#' 
#' @param robust An indicator of whether robust estimators of location and scale
#' should be used (\code{TRUE}, default), or not (\code{FALSE}). If 
#' \code{robust=TRUE}, location is estimated using the function 
#' \link[ICSNP:spatial.median]{spatial median} and scale is estimated using
#' the (unadjusted) median absolute deviation \link{MAD_R}. If 
#' \code{robust=FALSE}, mean is used for location, and the standard deviation 
#' for scale.
#' 
#' @param eps A small constant to check whether a vector is numerically zero or
#' not. By default \code{1e-16}.
#' 
#' @param vrs The program is given is two versions: (i) \code{vrs="C"} for the 
#' fast \code{C++} implementation (default), and (ii) \code{vrs="R"} for the
#' slower \code{R} version of the code. Up to small numerical errors, both
#' versions give identical results. 
#' 
#' @return A list composed of:
#' \itemize{
#'  \item{"values"}{ A numerical vector of length \code{q} with the 
#'  generalized eigenvalues. This output is analogous to the 
#'  output of function \link[base:eigen]{eigen}.}
#'  \item{"vectors"}{ A numerical matrix of dimension \code{d}-times-\code{q},
#'  each column for one generalized eigenvector corresponding to the 
#'  corresponding entry in \code{values}. This output is analogous to the 
#'  output of function \link[base:eigen]{eigen}.}
#'  \item{"PCA"}{ A numerical matrix of dimension \code{n}-times-\code{q} with
#'  coordinates of the data points from \code{X} in the basis of vectors from
#'  \code{vectors}.}
#'  \item{"location"}{ A numerical vector of length \code{d} with the location
#'  estimator for \code{X}.}
#' }
#' 
#' @seealso \link[rrcov:PcaProj]{PcaProj}, \link[rrcov:PcaClassic]{PcaClassic} 
#'
#' @examples
#' n = 500
#' d = 2
#' X = cbind(rnorm(n), rnorm(n))
#' 
#' ### Situation without contamination
#' 
#' res = PP_PCA(X, robust = FALSE) # non-robust version
#' resRob = PP_PCA(X) # robust version
#' 
#' plot(X, cex=.4, asp=1)
#' for(i in 1:d) arrows(res$location[1],res$location[2],
#'                      res$location[1] + res$vectors[1,i],
#'                      res$location[2] + res$vectors[2,i],
#'                      length = 0.1, col="orange", lwd=4)
#' for(i in 1:d) arrows(resRob$location[1],resRob$location[2],
#'                      resRob$location[1] + resRob$vectors[1,i],
#'                      resRob$location[2] + resRob$vectors[2,i],
#'                      length = 0.1, col="magenta", lwd=4)
#' legend("bottomleft",c("PCA","robPCA"),col=c("orange","magenta"),lwd=4)
#' 
#'  ### Situation with contamination
#'  alpha = .1
#'  m = floor(n*alpha)
#'  X[1:m,] = cbind(rnorm(m, 5, 1/10), rnorm(m, 5, 1/10))
#'  
#'  res = PP_PCA(X, robust = FALSE)
#'  resRob = PP_PCA(X)
#'  
#'  plot(X, cex=.4, asp=1)
#' for(i in 1:d) arrows(res$location[1],res$location[2],
#'                      res$location[1] + res$vectors[1,i],
#'                      res$location[2] + res$vectors[2,i],
#'                      length = 0.1, col="orange", lwd=4)
#' for(i in 1:d) arrows(resRob$location[1],resRob$location[2],
#'                      resRob$location[1] + resRob$vectors[1,i],
#'                      resRob$location[2] + resRob$vectors[2,i],
#'                      length = 0.1, col="magenta", lwd=4)
#' legend("bottomleft",c("PCA","robPCA"),col=c("orange","magenta"),lwd=4)

PP_PCA = function(X, q = NULL, robust = TRUE, 
                  eps = 1e-16, 
                  vrs = c("C", "R")){
  if(is.null(q)) q = ncol(X)
  n = nrow(X)
  d = ncol(X)
  vrs = match.arg(vrs, c("C", "R"))
  
  if(robust){
    center = ICSNP::spatial.median
    sfun = MAD_R
  } else {
    center = colMeans
    sfun = sd
  }
  
  mu = center(X)
  X = t(X)
  X = X - mu # centering
  
  if(vrs=="C"){
    res = PP_PCA_C(X, n, d, q, robust, eps = eps, echo = FALSE)
    lambda = c(res$values)
    v = res$vectors
    y = res$PCA
  }
  
  if(vrs=="R"){
    y = matrix(nrow=n,ncol=q)   # scores
    lambda = rep(NA,q)          # values
    v = matrix(nrow=d, ncol=q)  # vectors
  
    for(k in 1:q){
      if(k>1) X = X - outer(v[,k-1],y[,k-1])
      Xnorms = apply(X,2,function(x) sqrt(sum(x^2)))
      Xinds = (Xnorms<eps)
      Xmult = 1/Xnorms
      Xmult[Xinds] = 0 # if the norm is 0, project to origin
      A = t(X)*Xmult
      # A = t(apply(X,2,function(x) x/c(sqrt(sum(x^2)))))
      AX = A%*%X
      madk = apply(AX,1,sfun)
      kmax = which.max(madk)
      y[,k] = AX[kmax,]            # scores at max
      v[,k] = A[kmax,]             # maximum direction
      lambda[k] = madk[kmax]^2     # max value of mad
    }
  }
  return(list(values = lambda, vectors = v, PCA = y, location = mu))
}

#### all.equal with sign alignment ----
#' all.equal with Sign Alignment for Matrices
#'
#' Function \link[base:all.equal]{all.equal} that is adjusted so that possible
#' differences in sign between the columns of the two matrices are not taken 
#' into account.
#' 
#' @param A A numerical matrix of dimension \code{n}-times-\code{d}.
#' 
#' @param B A numerical matrix of dimension \code{n}-times-\code{d}.
#' 
#' @return A result of function \link[base:all.equal]{all.equal} applied to
#' aligned matrices \code{A} and \code{B}. Gives \code{TRUE} if the two matrices
#' are almost equal, up to the signs of their columns.
#' 
#' @seealso \link[base:all.equal]{all.equal}
#'
#' @examples
#' A = matrix(c(1,0,0, 0,1,0, 0,0,1),ncol=3, byrow=TRUE)
#' B = matrix(c(-1,0,0, 0,-1,0, 0,0,1),ncol=3, byrow=TRUE)
#' 
#' all.equal(A, B)
#' all_equal_S(A, B)

all_equal_S = function(A, B){
  q = ncol(A)
  s1 = apply(A,2,function(x) which.max(abs(x)))
  s2 = apply(B,2,function(x) which.max(abs(x)))
  s = sign(A[cbind(s1,1:q)]*B[cbind(s2,1:q)])
  B = B%*%diag(s)
  all.equal(A,B)
}
