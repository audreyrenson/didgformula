
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
generate_parameters <- function(Tt, mu_Beta_L=0.2, mu_Beta_A=0.2,
                                mu_Beta_Y=0.2, sd_Beta_L=0.2, sd_Beta_A=0.2,
                                sd_Beta_Y=0.2, range_ymeans = c(-5, 5),
                                constant_dydu = 0.1, check_CDF=FALSE) {

  Beta_U = -log(1 / 0.5 - 1)
  mean_U = plogis(Beta_U)

  # matrices of parameters
  Beta_L = matrix( rnorm(n = 2*(Tt+1), mean =mu_Beta_L, sd = sd_Beta_L), nrow = Tt + 1)# 2 because terms for A + intercept
  Beta_A = matrix( rnorm(n = 3*(Tt+1), mean =mu_Beta_A, sd = sd_Beta_A), nrow = Tt + 1)# 3 because terms for U, L + intercept
  Beta_Y = matrix( rnorm(n = 5*(Tt+1), mean =mu_Beta_Y, sd = sd_Beta_Y), nrow = Tt + 1)# 5 because terms for U, L, A, U*L + intercept

  Beta_Y = fix_dydu(Beta_Y, constant_dydu) #constant effect of U0 = part of parallel trends

  # balancing intercepts
  mean_Y = runif(Tt + 1, min = range_ymeans[1], max = range_ymeans[2])
  mean_L = runif(Tt + 1, min = 0.1, max = 0.5)
  mean_A_conditional = runif(Tt + 1, min = 1/(Tt^2), max = 1/(Tt^2) + 2/Tt) # this is a function of T to avoid very small pr(A_t=0) for large t
  mean_A_marginal = 1 - cumprod(1 - mean_A_conditional)  # basically a kaplan meier risk function

  Beta_Y[, 1] = mean_Y - Beta_Y[, 2]*mean_U - Beta_Y[, 3]*mean_L - Beta_Y[,4]*mean_A_marginal - Beta_Y[,5]*mean_U*mean_L #last piece only works for binary independent  L & U
  Beta_A[, 1] = -log(1 / mean_A_conditional - 1) - Beta_A[, 2]*mean_U - Beta_A[, 3]*mean_L
  Beta_L[, 1] = -log(1 / mean_L - 1) - Beta_L[, 2]*mean_A_marginal

  if(check_CDF) check_cdf(Beta_Y)

  return ( list(U=Beta_U, L=Beta_L, A=Beta_A, Y=Beta_Y))

}

fix_dydu <- function(Beta_Y, constant_dydu) {

    Beta_Y[, 2] = constant_dydu
    Beta_Y[, 5] = 0

    return ( Beta_Y )
}

check_cdf <- function(Beta_Y) {
  # for simulating based on hazards, check that the generated probabilities sum to less than one
  max_f <- Beta_Y
  max_f[,-1][max_f[,-1] < 0] = 0
  max_F = sum(plogis(rowSums(max_f)))
  if (max_F > 1) warning(glue::glue('If simulating time-to-event outcomes, you should reduce range_ymeans because parameters imply max(CDF)={max_F}.'))
}

#   } else if (ylink == 'rbinlogit') {
#
#
#     # find the conditional ORs for Y by U that yeild a constant risk difference
#     # this doesn't quite work and I need to put this problem away for the moment, my brain isn't working
#     # currently there are a lot of hard-coded indices here, make this more dynamic/explicit
#     combos <- expand.grid(t = 1:nrow(Beta_Y) - 1, int = 1,  L=c(0,1)) #don't need to include A because we only care about when A==0
#     combos$prY_Uequals0 = plogis( Beta_Y[combos$t+1, 1] + Beta_Y[combos$t+1, 2]*combos$L )
#     combos$prY_Uequals1 = combos$prY_Uequals0 + constant_dydu
#
#     #then figure out what odds ratios are needed
#     log_odds_ratio <- function(p1, p0) log( (p1/(1-p1)) / (p0/(1-p0)) )
#     combos$log_or_YU = log_odds_ratio(combos$prY_Uequals1, combos$prY_Uequals0)
#
#     #then fill in the right terms
#     Beta_Y[, 2] = combos$log_or_YU[combos$L == 0]
#     Beta_Y[, 5] = combos$log_or_YU[combos$L == 1] - combos$log_or_YU[combos$L == 0]
#
#     return(Beta_Y)
#   } else {
#     stop("ylink must be either 'identity' or 'logit'")
#   }
#
# }
