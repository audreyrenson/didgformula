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

rbinom_logit <- function(X, Beta, N, n=1) rbinom(n=N, p=plogis(as.matrix(X) %*% Beta), size=n)
rnorm_identity <- function(X, Beta, N, sd=1) rnorm(n=N, mean = as.matrix(X) %*% Beta, sd)
