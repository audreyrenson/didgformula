#' Generate simulated data that meets the parallel trends assumption
#'
#' @param N int. Number of independent observations
#' @param Tt int. Number of periods, minus 1. I.e. there are Tt + 1 periods.
#' @param Beta list of length 4. Output of generate_parameters().
#' @param potential_outcomes logical. Should outcomes and covariates be generated with exposure set to 0 at all times?
#' @param ylink chr. One of "rnorm_identity", "rbinom_logit", or "rbinom_logit_hazard".
#' @param binomial_n int length N. Defaults to all 1's. If ylink is rbinom_logit, you can optionally pass a vector of group sizes to generate aggregate binomial data. In this case, treatments and covariates will be constant at the group level for a given time period.
#'
#' @return Data frame with N rows and (Tt+1)*3 + 2 columns - 'uid' is a unique identifier, 'U0' is an 'unmeasured' baseline covariate, L{t},A{t},Y{t} are covariates, exposures, and outcomes, respectively. If binomial_n != 1, an additional column binomial_n is also included.
#' @export
#'
#' @examples
generate_data <- function(N,
                          Tt,
                          Beta,
                          potential_outcomes=FALSE,
                          ylink = "rnorm_identity",
                          binomial_n=1){

  df = data.frame(uid       = seq_len(N),
                  intercept = 1,
                  zeros     = 0,
                  U0        = rbinom_logit(X=1, Beta=Beta$U, N=N))

  #what variables does each variable depend on?
  vars_L = function(t) c('intercept', if(t<1 | potential_outcomes) 'zeros' else glue::glue('A{t-1}'))
  vars_A = function(t) c('intercept', 'U0', glue::glue('L{t}'))
  vars_Y = function(t) c('intercept', 'U0', glue::glue('L{t}'), if(potential_outcomes) 'zeros' else glue::glue('A{t}'), glue::glue('U0timesL{t}'))


  for(t in 0:Tt) {
    df[[glue::glue('L{t}')]] = rbinom_logit(X=df[vars_L(t)], Beta=Beta$L[t+1, ], N=N)
    df[[glue::glue('U0timesL{t}')]] = df[[glue::glue('L{t}')]] * df$U0

    if (t<1) {
      df[[glue::glue('A{t}')]] = 0
    } else {
      df[[glue::glue('A{t}')]] = rbinom_logit(X=df[vars_A(t)], Beta=Beta$A[t+1, ], N=N)^(1 - df[[glue::glue('A{t-1}')]] == 1) #monotonic treatment assignment
    }
  }

  df = cbind(df, simulate_y(df=df, Tt=Tt, vars_Y=vars_Y, Beta_Y=Beta$Y, ylink=ylink, binomial_n=binomial_n))

  if(length(binomial_n) > 1) df$binomial_n = binomial_n

  return( df[, -grep( c('intercept|zeros|U0times'), names(df))] )
}

# generate data from a binomial distribution based on a linear-logistic model
rbinom_logit <- function(X, Beta, N, n=1) rbinom(n=N, p=plogis(as.matrix(X) %*% Beta), size=n)

simulate_y <- function(df, Tt, vars_Y, Beta_Y,
                       ylink='rnorm_identity', binomial_n, ...) {

  eta_ti = sapply(0:Tt, function(t) as.matrix(df[vars_Y(t)]) %*% Beta_Y[t+1, ])

  if(ylink == 'rnorm_identity') {
    Y_ti = apply(eta_ti, 2, function(mu) rnorm(nrow(df), mean=mu, ...))
  } else if (ylink=='rbinom_logit') {
    Y_ti = apply(eta_ti, 2, function(expitp) rbinom(n=nrow(df), prob = plogis(expitp), size=binomial_n))
  } else if (ylink=="rbinom_identity") {
    Y_ti = apply(eta_ti, 2, function(p) rbinom(n=nrow(df), prob = p, size=binomial_n))
  } else if (ylink=='rbinom_logit_hazard') {
    p_ti = plogis(eta_ti)
    S_ti = 1 - matrixStats::rowCumsums(p_ti)
    h_ti = p_ti / S_ti
    Y_ti = apply(h_ti, MARGIN=2, FUN=function(h) rbinom(nrow(df), prob=h, size=binomial_n))
    already_had_event = cbind(FALSE, matrixStats::rowCummaxs(Y_ti)[, -ncol(Y_ti)] == 1 )
    Y_ti [already_had_event] = 0
  } else {
    stop('ylink must be one of "rnorm_identity", "rbinom_identity", "rbinom_logit", or "rbinom_logit_hazard"')
  }
  colnames(Y_ti) = sapply(0:Tt, function(t) glue::glue('Y{t}'))
  return (as.data.frame(Y_ti))
}


