#' Find c2 for a Bivariate Normal Tail Probability
#'
#' Given a correlation \eqn{\rho}, this function finds \eqn{c_2} such that
#' \eqn{P(Z_1 > c_1 \text{ or } Z_2 > c_2) = \alpha}, where
#' \eqn{(Z_1, Z_2)} follows a bivariate normal distribution with mean 0,
#' unit variances, and correlation \eqn{\rho}.
#'
#' Equivalently, it solves
#' \eqn{1 - \Phi_2(c_1, c_2; \rho) - \alpha = 0},
#' where \eqn{\Phi_2} is the bivariate normal CDF.
#'
#' @param rho Correlation parameter. Must satisfy \eqn{|\rho| < 1}.
#' @param c1 Fixed cutoff for the first normal variable. Default is 2.390.
#' @param alpha Target upper-tail probability. Default is 0.025.
#'
#' @return A numeric scalar giving the value of \eqn{c_2}.
#'
#' @examples
#' find_c2(rho = 0.3)
#' find_c2(rho = 0.0, c1 = 2.390, alpha = 0.025)
#'
#' @importFrom mvtnorm pmvnorm
#' @export


find_c2 <- function(rho, c1 = 2.390, alpha = 0.025) {
  stopifnot(abs(rho) < 1)
  
  f <- function(c2) {
    # Phi_2(c1, c2; rho) = P(Z1 <= c1, Z2 <= c2)
    Phi2 <- mvtnorm::pmvnorm(
      upper = c(c1, c2),
      corr = matrix(c(1, rho, rho, 1), nrow = 2),
      algorithm = mvtnorm::TVPACK()
    )[1]
    
    1 - Phi2 - alpha
  }
  
  # widen interval if needed
  uniroot(f, interval = c(-10, 10), tol = 1e-10)$root
}
