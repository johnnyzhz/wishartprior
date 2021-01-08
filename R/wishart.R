#' Obtain the mean and variance of a Wishart distribution
#'
#'#'This function calculates the mean and variance of a
#' Wishart distribution based on the input of a scale matrix
#' and the degrees of freedom.
#'
#' @param sigma The scale matrix.
#' @param m Degrees of freedom.
#' @return The mean and variance.
#' @examples
#' sigma <- array(c(2,0.5,0.5,1), dim=c(2,2))
#' stat.Wishart(sigma, 10)
stat.Wishart <- function(sigma, m){
  m.wishart <- m*sigma
  diag.sigma <- diag(sigma)
  v.wishart <- m*sigma^2 + m*diag.sigma%*%t(diag.sigma)
  list(mean=m.wishart, variance=v.wishart)
}

#' Obtain the mean and variance of an inverse Wishart distribution
#'
#'This function calculates the mean and variance of an inverse
#' Wishart distribution based on the input of a scale matrix
#' and the degrees of freedom.
#'
#' @param sigma The scale matrix.
#' @param m Degrees of freedom.
#' @return A 3d array of generated matrices.
#' @examples
#' sigma <- array(c(2,0.5,0.5,1), dim=c(2,2))
#' stat.iWishart(sigma, 10)
stat.iWishart <- function(sigma, m){
  p <- nrow(sigma)
  m.iwishart <- sigma/(m-p-1)
  diag.sigma <- diag(sigma)
  v.iwishart <- (m-p+1)*sigma^2 + (m-p-1)*diag.sigma%*%t(diag.sigma)
  v.iwishart <- v.iwishart/((m-p)*(m-p-1)^2*(m-p-3))
  list(mean=m.iwishart, variance=v.iwishart)
}

#' Generate random numbers from an inverse Wishart distribution
#'
#'This function generates random matrix from an inverse
#' Wishart distribution.
#'
#' @param n Sample size.
#' @param Sigma The scale matrix.
#' @param df Degrees of freedom.
#' @return The mean and variance.
#' @seealso \code{\link{rWishart}} for random matrix from a Wishart distribution.
#' @examples
#' n <- 1000
#' df <- 10
#' Sigma <- array(c(2,0.5,0.5,1), dim=c(2,2))
#' dset1 <- rWishart(n, df, Sigma)
#'
#' apply(dset1, 1:2, mean)  ## mean
#' apply(dset1, 1:2, var)   ## variance
riWshart <- function(n, df, Sigma){
  Sigma.inv <- solve(Sigma)
  temp <- rWishart(n, df, Sigma.inv)
  temp <- lapply(1:n, function(x) solve(temp[ , , x]))
  array(unlist(temp), dim = c(dim(temp[[1]]), length(temp)))
}

#' Posterior based on an inverse Wishart distribution for covariance matrix
#'
#'Combine data and inverse Wishart prior of covariance matrix.
#'
#' @param n Sample size.
#' @param S Sample covariance matrix (MLE).
#' @param m0 Degrees of freedom of the inverse Wishart prior.
#' @param V0 Scale matrix of the inverse Wishart prior.
#' @return The mean and variance of the inverse Wishart distribution.
#' @seealso \code{\link{posterior.Wishart}} precision matrix.
posterior.iWishart <- function(n, S, m0, V0){
  m1 <- n + m0
  V1 <- n*S + V0
  stat.iWishart(V1, m1)
}

#' Posterior based on a Wishart distribution for precision matrix
#'
#'Combine data and Wishart prior of precision matrix.
#'
#' @param n Sample size.
#' @param S Sample covariance matrix (MLE).
#' @param w0 Degrees of freedom of the Wishart prior.
#' @param U0 Scale matrix of the Wishart prior.
#' @return The mean and variance of both Wishart and inverse Wishart distribution.
#' @seealso \code{\link{posterior.iWishart}} covariance matrix.
posterior.Wishart <- function(n, S, w0, U0){
  w1 <- n + w0
  U1 <- solve(n*S + solve(U0))
  V1 <- n*S + solve(U0)
  list(Wishart=stat.Wishart(U1, w1),
       iWishart=stat.iWishart(V1, w1))
}

#' Calculate the degrees of freedom based on the mean and variance
#'
#'One can input the expectation and variance of a diagonal element of the matrix.
#'
#' @param m Mean.
#' @param var Variance.
#' @param p The dimension of the matrix.
#' @return The degrees of freedom of an inverse Wishart distribution.
df <- function(m, var, p){
  2*m^2/var + p + 3
}
