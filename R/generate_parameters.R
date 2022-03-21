
#' Pick parameter values for simulation, from a normal distribution
#'
#' @param Tt Final period (t=0,1,...,Tt)
#' @param mu_Beta_L numeric. mean for parameters in covariate models
#' @param mu_Beta_A numeric. mean for parameters in treatment models
#' @param mu_Beta_Y numeric. mean for parameters in outcome models
#' @param sd_Beta_L numeric. standard deviation for parameters in covariate models
#' @param sd_Beta_A numeric. standard deviation for parameters in treatment models
#' @param sd_Beta_Y numeric. standard deviation for parameters in outcome models
#' @param range_ymeans
#' @param constant_dydu
#' @param ylink
#'
#' @return
#' @export
#'
#'
#' @examples
generate_parameters <- function(Tt,
                                mu_Beta_L=0.2, mu_Beta_W=0.2, mu_Beta_A=0.2, mu_Beta_Y=0.2,
                                sd_Beta_L=0.2, sd_Beta_W=0.2, sd_Beta_A=0.2, sd_Beta_Y=0.2,
                                range_ymeans = c(-5, 5),
                                constant_dydu = 0.1,
                                check_CDF=FALSE) {

  Beta_U = -log(1 / 0.5 - 1)
  mean_U = plogis(Beta_U)

  # matrices of parameters
  Beta_L = matrix( rnorm(n = 2*(Tt+1), mean =mu_Beta_L, sd = sd_Beta_L), nrow = Tt + 1, byrow = TRUE)# 2 because terms for A + intercept
  Beta_W = matrix( rnorm(n = 2*(Tt+1), mean =mu_Beta_W, sd = sd_Beta_W), nrow = Tt + 1, byrow = TRUE)# 2 because terms for A + intercept
  Beta_A = matrix( rnorm(n = 5*(Tt+1), mean =mu_Beta_A, sd = sd_Beta_A), nrow = Tt + 1, byrow = TRUE)# 5 because terms for U, L, W, W^2, + intercept
  Beta_Y = matrix( rnorm(n = 6*(Tt+1), mean =mu_Beta_Y, sd = sd_Beta_Y), nrow = Tt + 1, byrow = TRUE)# 6 because terms for U, L, W, W^2, A + intercept

  Beta_Y[,2] = constant_dydu #constant effect of U0 = part of parallel trends

  # balancing intercepts
  mean_Y = runif(Tt + 1, min = range_ymeans[1], max = range_ymeans[2])
  mean_L = runif(Tt + 1, min = 0.1, max = 0.5)
  mean_W = 0
  mean_A_conditional = runif(Tt + 1, min = 1/(Tt^2), max = 1/(Tt^2) + 1/Tt) # this is a function of T to avoid very small pr(A_t=0) for large t
  mean_A_marginal = 1 - cumprod(1 - mean_A_conditional)  # basically a kaplan meier risk function

  Beta_Y[, 1] = mean_Y - Beta_Y[, 2]*mean_U - Beta_Y[, 3]*mean_L - Beta_Y[,4]*mean_W - Beta_Y[,5]*mean_W - Beta_Y[,6]*mean_A_marginal
  Beta_A[, 1] = -log(1 / mean_A_conditional - 1) - Beta_A[, 2]*mean_U - Beta_A[, 3]*mean_L
  Beta_L[, 1] = -log(1 / mean_L - 1) - Beta_L[, 2]*mean_A_marginal #fix - should depend on At-1 not At
  Beta_W[, 1] = mean_W - Beta_W[, 2]*mean_A_marginal #fix - should depend on At-1 not At

  colnames(Beta_Y) = c('intercept', 'U0', 'Lt', 'Wt','Wt^2','At')
  colnames(Beta_A) = c('intercept', 'U0', 'Lt', 'Wt','Wt^2')
  colnames(Beta_L) = colnames(Beta_W) = c('intercept', 'At-1')


  if(check_CDF) check_cdf(Beta_Y)

  return ( list(U=Beta_U, L=Beta_L, W=Beta_W, A=Beta_A, Y=Beta_Y))

}

check_cdf <- function(Beta_Y) {
  # for simulating based on hazards, check that the generated probabilities sum to less than one
  max_f <- Beta_Y
  max_f[,-1][max_f[,-1] < 0] = 0
  max_F = sum(plogis(rowSums(max_f)))
  if (max_F > 1) warning(glue::glue('If simulating time-to-event outcomes, you should reduce range_ymeans because parameters imply max(CDF)={max_F}.'))
}

check_rbinom_identity <- function(Beta_Y, silently=FALSE) {
  #throw error if any pr(y|x) outside (0,1)
  Beta_Y_min = Beta_Y_max = Beta_Y
  Beta_Y_min[,-1][Beta_Y_min[,-1] > 0] = 0
  Beta_Y_max[,-1][Beta_Y_max[,-1] < 0] = 0

  below_0 = any(rowSums(Beta_Y_min) < 0)
  exceeds_1 = any(rowSums(Beta_Y_max) > 1)

  beta_status = below_0 + exceeds_1

  if(silently) {
    return (beta_status == 0)
  } else if(below_0 & exceeds_1) {
    stop('Beta_Y implies probabilities both <0 and >1')
  } else if (below_0 & !exceeds_1) {
    stop('Beta_Y implies probabilities <0')
  } else if (!below_0 & exceeds_1) {
    stop('Beta_Y implies probabilities >1')
  } else {
    return (TRUE)
  }
}

generate_parameters_until_noerror <- function(..., maxiter=20, original_maxiter=maxiter) {

  Beta = generate_parameters(..., check_rbinom_identity=FALSE)

  if( maxiter == 0 )  {
    stop('Maximum number of iterations reached without a solution')
  } else if(!check_rbinom_identity(Beta$Y, silently=TRUE)) { #solution not found yet, so iterate
      recursive_args <- as.list(match.call()[-1])
      recursive_args[['original_maxiter']] = original_maxiter
      recursive_args[['maxiter']] <- maxiter - 1
      return ( do.call(generate_parameters_until_noerror, recursive_args) )
  } else { #solution found
    return ( c(Beta, list(iterations = original_maxiter - maxiter + 1)) )
  }
}

