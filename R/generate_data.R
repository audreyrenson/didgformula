#' Generate simulated data that meets the parallel trends assumption
#'
#' @param N int. Number of independent observations
#' @param Tt int. Number of periods, minus 1. I.e. there are Tt + 1 periods.
#' @param Beta list of length 4. Output of generate_parameters().
#' @param potential_outcomes logical. Should outcomes and covariates be generated with exposure set to 0 at all times?
#'
#' @return Data frame with N rows and (Tt+1)*3 + 2 columns - 'uid' is a unique identifier, 'U0' is an 'unmeasured' baseline covariate, L{t},A{t},Y{t} are covariates, exposures, and outcomes, respectively.
#' @export
#'
#' @examples
generate_data <- function(N, Tt, Beta, potential_outcomes=FALSE){

  df = data.frame(uid       = seq_len(N),
                  intercept = 1,
                  zeros     = 0,
                  U0        = rbinom_logit(X=1, Beta=Beta$U, N=N))

  #what variables does each variable depend on?
  vars_L = function(t) c('intercept', if(t<1 | potential_outcomes) 'zeros' else glue::glue('A{t-1}'))
  vars_A = function(t) c('intercept', 'U0', glue::glue('L{t}'))
  vars_Y = function(t) c('intercept', 'U0', glue::glue('L{t}'), if(potential_outcomes) 'zeros' else glue::glue('A{t}'), 'U0timesLt')


  for(t in 0:Tt) {
    df[[glue::glue('L{t}')]] = rbinom_logit(X=df[vars_L(t)], Beta=Beta$L[t+1, ], N=N)
    df[['U0timesLt']] = df[[glue::glue('L{t}')]] * df$U0

    if (t<1) {
      df[[glue::glue('A{t}')]] = 0
    } else {
      df[[glue::glue('A{t}')]] = rbinom_logit(X=df[vars_A(t)], Beta=Beta$A[t+1, ], N=N)^(1 - df[[glue::glue('A{t-1}')]] == 1) #monotonic treatment assignment
    }

    df[[glue::glue('Y{t}')]] = rnorm_identity(X=df[vars_Y(t)], Beta=Beta$Y[t+1, ], N=N)

  }
  return( df[, -which( names(df) %in% c('intercept','zeros', 'U0timesLt') ) ] )
}

# generate data from a binomial distribution based on a linear-logistic model
rbinom_logit <- function(X, Beta, N, n=1) rbinom(n=N, p=plogis(as.matrix(X) %*% Beta), size=n)
# generate data from a normal distribution based on a linear-identity model
rnorm_identity <- function(X, Beta, N, sd=1) rnorm(n=N, mean = as.matrix(X) %*% Beta, sd)
