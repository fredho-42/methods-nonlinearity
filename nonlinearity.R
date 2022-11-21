library(tidyverse)
library(ggplot2)
library(mgcv)
library(splines)


# Data prep ---------------------------------------------------------------

dat <- read_tsv('../UKDA-8715-tab/tab/longitudinal_td.tab') %>% 
  select(pidp, indscus_lw_9, strata, psu,             # select variables 
         wave, age_dv, sex_dv, ethn_dv, hiqual_dv, 
         bmi_dv, sf12pcs_dv, sf12mcs_dv) %>% 
  filter(wave %in% c(1, 9)) %>%                       # include only Waves 1 and 9
  pivot_wider(id_cols=c(pidp, indscus_lw_9,           # wide data for regression
                        strata, psu),                 
              names_from=wave, values_from=age_dv:sf12mcs_dv) %>% 
  select(-c(age_dv_9, sex_dv_9, ethn_dv_9,            # remove unused variables 
            hiqual_dv_9, bmi_dv_9)) %>% 
  na.omit()                                           # remove missing data

dat <- dat %>% 
  mutate(sex_dv_1=factor(sex_dv_1),                     # sex as factor
         ethn_dv_1=factor(if_else(ethn_dv_1==1, 1, 2)), # recode ethnicity because of small subsample
         hiqual_dv_1=factor(hiqual_dv_1),               # education as factor 
         obese=factor(if_else(                          # dichomtomised obesity
           bmi_dv_1>=30, 'Obese', 'Non-obese')), 
         bmi_cat=factor(case_when(                      # categorised obesity
           bmi_dv_1<18.5~ 'Underweight',     
           bmi_dv_1<25  ~ 'Normal weight', 
           bmi_dv_1<30  ~ 'Overweight', 
           bmi_dv_1<99  ~ 'Obese', 
           TRUE         ~ NA_character_, 
         )))


# Models ------------------------------------------------------------------
# adjusted for age, baseline MCS and PCS (P-splines)
# sex, ethnicity, education, 
# strata and PSU (random effects)
# sampling weights provided by US

m.1a  <- gam(sf12mcs_dv_9 ~ bmi_dv_1 +                                                      # Linear
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat)

m.1b  <- gam(sf12mcs_dv_9 ~ obese +                                                         # Binary 
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat)

m.1c  <- gam(sf12mcs_dv_9 ~ bmi_cat +                                                       # Categorical 
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat)

m.2a  <- gam(sf12mcs_dv_9 ~ bmi_dv_1 + I(bmi_dv_1^2) + I(bmi_dv_1^3) +                      # Cubic
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat) 

m.2b1 <- gam(sf12mcs_dv_9 ~ ns(bmi_dv_1, knots=seq(20, 40, 10), Boundary.knots=c(15, 55)) + # NCS (3 knots)
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat) 

m.2b2 <- gam(sf12mcs_dv_9 ~ ns(bmi_dv_1, knots=seq(15, 54, 3), Boundary.knots=c(15, 55)) + # NCS (14 knots)
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat) 

m.2c  <- gam(sf12mcs_dv_9 ~ s(bmi_dv_1, bs='ps') +                                         # P-spline
               s(sf12pcs_dv_1, bs='ps') + s(sf12mcs_dv_1, bs='ps') + s(age_dv_1, bs='ps') + 
               sex_dv_1 + ethn_dv_1 + hiqual_dv_1 + 
               s(strata, bs='re') + s(psu, bs='re'), weights=indscus_lw_9, data=dat) 


# Data prep for plots -------------------------------------------------------

dat.plot <- tibble(                              # create data for plotting
  
  m=rep(c('poly', 'ncs1', 'ncs2'), each=337),    # model labels
  
  x=rep(termplot(m.2a, se=T, plot=F)[[1]]$x, 3), # x values
  
  y=c(termplot(m.2a , se=T, plot=F)[[1]]$y +     # y values from linear term
        termplot(m.2a , se=T, plot=F)[[2]]$y +   # y values from quadratic term
        termplot(m.2a , se=T, plot=F)[[3]]$y,    # y values from cubic term
      termplot(m.2b1, se=T, plot=F)[[1]]$y,      # y values from NCS 1
      termplot(m.2b2, se=T, plot=F)[[1]]$y))     # y values from NCS 2

dat.plot <- dat.plot %>% 
  bind_rows(
    with(plot(m.2c, select=1)[[1]], 
         tibble(m='pspl', x, y=fit[, 1]))        # x and y values from p splines
  ) %>% 
  group_by(m) %>% 
  mutate(y=y-y[which(x>=20 & x<=20.05859)]) %>% # centre y at BMI=20
  ungroup()


dat.plot <- dat.plot %>% bind_rows(
  dat.plot, 
  termplot(m.1a, plot=F)$bmi_dv_1 %>%               # y values from assumed linear
    mutate(m='lin', y=y-y[which(x>=20 & x<=20.05859)]),
  termplot(m.1b, plot=F)$obese[c(1:2, 2), ] %>% 
    mutate(m='bin', x=c(12, 30, 50)),               # y values from dichotomised
  termplot(m.1c, plot=F)$bmi_cat[c(4, 1, 3, 2, 2), ] %>% 
    mutate(m='cat', x=c(12, 18.5, 25, 30, 50))) %>% # y values from categorised
  
  mutate(m=factor(m, levels=c('lin', 'bin', 'cat', 
                              'poly', 'ncs1', 'ncs2', 'pspl')))


# Plots -------------------------------------------------------------------

fig1 <- ggplot(data=dat.plot %>% 
                 filter(m %in% c('lin', 'bin', 'cat')) , 
               aes(x=x, y=y, linetype=m)) + 
  # geom_hline(yintercept =0, linetype='dashed') + 
  geom_step(data=dat.plot %>% filter(m %in% c('bin', 'cat')), aes(colour=m), size=0.8, alpha=0.7) + 
  geom_line(data=dat.plot %>% filter(m %in% c('lin')), aes(colour=m), size=0.8, alpha=0.7) + 
  scale_y_continuous(name='Coefficient for mental wellbeing', limits=c(-8, 2)) + 
  scale_x_continuous(name=bquote('BMI (kg/'*m^2*')'), limits=c(10, 50)) + 
  scale_color_brewer(palette = "Dark2", 
                     name='Method', 
                     limits=c('lin', 'bin', 'cat'), 
                     labels=c('Assumed linear', 
                              'Dichotomised', 
                              'Categorised')) + 
  scale_linetype_discrete(name='Method', 
                          limits=c('lin', 'bin', 'cat'), 
                          labels=c('Assumed linear', 
                                   'Dichotomised', 
                                   'Categorised')) + 
  theme_bw() + 
  theme(legend.position=c(.99, .01), legend.justification=c(1, 0))


fig2 <- ggplot(data=dat.plot %>% filter(!m %in% c('lin', 'bin', 'cat')), 
               aes(x=x, y=y, linetype=m)) + 
  # geom_hline(yintercept =0, linetype='dashed') + 
  geom_line(aes(colour=m)) + 
  scale_y_continuous(name='Coefficient for mental wellbeing', limits=c(-8, 2)) + 
  scale_x_continuous(name=bquote('BMI (kg/'*m^2*')'), limits=c(10, 50)) + 
  scale_color_brewer(palette = "Dark2", 
                     name='Method', 
                     labels=c('Cubic polynomial', 
                              'Natural cubic spline (3 knots)', 
                              'Natural cubic spline (14 knots)', 
                              'P-spline')) + 
  scale_linetype_discrete(name='Method', 
                          labels=c('Cubic polynomial', 
                                   'Natural cubic spline (3 knots)', 
                                   'Natural cubic spline (14 knots)', 
                                   'P-spline')) + 
  theme_bw() + 
  theme(legend.position=c(.99, .01), legend.justification=c(1, 0))



ggsave('./fig1.tiff', fig1, height=1250, width=1250, units='px', compression='lzw')
ggsave('./fig2.tiff', fig2, height=1250, width=1250, units='px', compression='lzw')
