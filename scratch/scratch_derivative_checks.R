Beta <- generate_parameters(Tt, ylink='logit')
all_combos <- expand_grid(tt = 1:Tt, U=c(0,1), Ltmin1=c(0,1), Lt=c(0,1))
all_combos %>%
  mutate(eta_tmin1 =  Beta[tt, 1] + U*Beta[tt, 2] + Ltmin1*Beta[tt, 3]  + U*Ltmin1*Beta[tt, 5],
         eta_t = Beta[tt + 1, 1] + U*Beta[tt+1, 2] + Lt*Beta[tt+1, 3]  + U*Lt*Beta[tt+1, 5],
         EYtmin1 = plogis(eta_tmin1),
         EYt = plogis(eta_t))
  arrange(Lt, Ltmin1, tt) %>%
  group_by(Lt, Ltmin1, tt) %>%
  summarise(dydu_t = diff(EYt), dydu_tmin1 = diff(EYtmin1)) %>% View()




expand_grid(tt = 1:Tt, Ltmin1=c(0,1), Lt=c(0,1)) %>%
  mutate(pYt_U0 =plogis( Beta_Y[tt+1, 1]  + Beta_Y[tt+1, 2]*0  + Beta_Y[tt+1, 3]*Lt + Beta_Y[tt+1, 5]*0 ),
         pYt_U1 = plogis( Beta_Y[tt+1, 1]  + Beta_Y[tt+1, 2]*1  + Beta_Y[tt+1, 3]*Lt + Beta_Y[tt+1, 5]*0 ),
         pYtmin1_U0 = plogis( Beta_Y[tt, 1]  + Beta_Y[tt, 2]*0  + Beta_Y[tt, 3]*Ltmin1 + Beta_Y[tt, 5]*0 ),
         pYtmin1_U1 = plogis( Beta_Y[tt, 1]  + Beta_Y[tt, 2]*1  + Beta_Y[tt, 3]*Ltmin1 + Beta_Y[tt, 5]*Ltmin1),
         difft = pYt_U1 - pYt_U0,
         difftmin1 = pYtmin1_U1 - pYtmin1_U0)



  as_tibble(Beta_Y) %>%
    magrittr::set_colnames(c('int','bU','bL','bA','bUL')) %>%
    rownames_to_column(var='t') %>%
    mutate(t=as.numeric(t)-1) %>%
    mutate(prY0_U0_L0 = plogis(int),
           prY0_U0_L1 = plogis(int + bL),
           prY0_U1_L0 = plogis(int + bU),
           prY0_U1_L1 = plogis(int + bU + bL + bUL),
           dydu_L1    = prY0_U1_L1 - prY0_U0_L1,
           dydu_L0    = prY0_U1_L0 - prY0_U0_L0)

  #so we have that  plogis(int + bU) - plogis(int) is the same for all t
  #we also need plogis(int + bU + bL + bUL) - plogis(int + bL) is the same.







  fix_betat <- function(beta_t, beta_tmin1) {

    beta_t[5] = - beta_t[3] - beta_t[1] +
      log(  exp(beta_t[3]) + plogis(beta_tmin1[1])/plogis(beta_t[1]) *
              ( exp(beta_tmin1[2] + beta_tmin1[3] + beta_tmin1[5]) - exp(beta_tmin1[2])  )  )

    beta_t[2] = beta_tmin1[2]

  }

  Beta_Y[2, ] = fix_betat(Beta_Y[2, ], Beta_Y[1, ])
  Beta_Y[3, ] = fix_betat(Beta_Y[3, ], Beta_Y[2, ])
  Beta_Y[4, ] = fix_betat(Beta_Y[4, ], Beta_Y[3, ])
  Beta_Y[5, ] = fix_betat(Beta_Y[5, ], Beta_Y[4, ])
