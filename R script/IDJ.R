# This is the R script belonging to the article + Supplement of " Kinetics of heat-induced changes in dairy products: developments in data analysis and modeling techniques" by M.A.J.S. van Boekel, published in the International Dairy Journal, 2021, doi 10.1016/j.idairyj.2021.105187.

# Packages used:

library(brms)
library(tidyverse)
library(ggExtra)
library(papaja)
library(tidybayes)
library(GGally)
library(ggridges)
library(rstan)
library(patchwork)
library(here)
library(broom)
library(forcats)
library(modelr)
library(qqplotr)
library(bookdown)

theme_set(theme_bw())

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

# settings:
knitr::opts_chunk$set(echo = FALSE, 
                      fig.align = "center",
                      fig.width = 6, 
                      fig.height = 4, 
                      dev = "png",
                      results="asis",
                      warning = FALSE,
                      message = FALSE,
                      comment = NA)

# Code for Figure 1

t_sim <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
df_fo1 <- data.frame(t_sim) %>% mutate(ct1=exp(-0.2*t_sim)) %>% mutate(ct2=0.8*exp(-0.2*t_sim)) %>% mutate(ct3=0.6*exp(-0.2*t_sim)) %>% mutate(ct4=0.4*exp(-0.2*t_sim)) %>% mutate(ct5=0.2*exp(-0.2*t_sim))

df_fo1_plot <- df_fo1 %>% ggplot(aes(x=t_sim))+
  geom_line(aes(y=ct1))+
  geom_line(aes(y=ct2))+
  geom_line(aes(y=ct3))+
  geom_line(aes(y=ct4))+
  geom_line(aes(y=ct5))+
  labs(x="time", y=expression(paste(c[t])), subtitle = "A")+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))

df_fo2 <- data.frame(t_sim) %>% mutate(ct1=exp(-0.02*t_sim)) %>% mutate(ct2=exp(-0.07*t_sim)) %>% mutate(ct3=exp(-0.12*t_sim)) %>% mutate(ct4=exp(-0.2*t_sim)) %>% mutate(ct5=exp(-0.4*t_sim))

df_fo2_plot <- df_fo2 %>% ggplot(aes(x=t_sim))+
  geom_line(aes(y=ct1))+
  geom_line(aes(y=ct2))+
  geom_line(aes(y=ct3))+
  geom_line(aes(y=ct4))+
  geom_line(aes(y=ct5))+
  labs(x="time", y=expression(paste(c[t])), subtitle="B")+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))

df_fo1_plot+df_fo2_plot

# Code for Figure 2:

t_sim <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
df_order <- data.frame(t_sim) %>% mutate(c1=1-0.2*t_sim) %>% mutate(c2=exp(-0.2*t_sim)) %>% mutate(c3=1/(1+0.2*t_sim))
df_order %>% ggplot(aes(x=t_sim))+
  geom_line(aes(y=c1), lty=1)+
  geom_line(aes(y=c2),lty=2)+
  geom_line(aes(y=c3), lty=3)+
  labs(x="time, s", y="c, mol/L")+
  ylim(c(0,1))+
  annotate("text",x=5, y=0.56, label=expression(paste(n[t],"=2")))+
  annotate("text",x=5, y=0.43, label=expression(paste(n[t],"=1")))+
  annotate("text",x=5, y=0.2, label=expression(paste(n[t],"=0")))+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  ylim(0,1)

# HMF CASE STUDY 1

# Code for Figure 3

HMF <- read.csv(here("data", "HMF_Fink.csv"),header=TRUE, sep=";")
#HMF$serial <- as.factor(HMF$serial)
HMF_plot <- HMF %>% mutate(serial=str_c("Experiment ", serial))
(HMF_plot_facet <- HMF_plot %>% ggplot(aes(x=time, y=conc))+
    geom_point(shape=21, fill="red", size=1.5, stroke=1)+
    facet_wrap(~serial)+
    labs(x="time in days", y=expression(paste("HMF in  ", mu,"mol L"^-1)))+
    theme(strip.background = element_rect(color="black", fill="lightblue", size=1.5, linetype="solid"))
)

# Code for Figure S1:

set.seed(1)
sample_c0 <- rnorm(100,0,5)
sample_kr <-  rnorm(100,0.05,0.5)
sim_time <- seq(from = 0, to = 300, length.out = 100)
prior_sim <- data.frame(sim_time, sample_c0, sample_kr)

ggplot(data=prior_sim, aes(x=sim_time))+
  geom_abline(data=prior_sim, aes(intercept = sample_c0, slope=sample_kr, group=sim_time), colour="blue")+
  xlim(c(0,200))+
  ylim(c(-100,300))+
  labs(x="time", y="predicted HMF concentration")

# Code for Figure S2:

set.seed(5)
sample_c0 <- rnorm(100,0,5)
sample_kr <-  truncnorm::rtruncnorm(100,a=0, mean=0.05,sd=0.5)
sim_time <- seq(from = 0, to = 206, length.out = 100)
prior_sim <- data.frame(sim_time, sample_c0, sample_kr)

ggplot(data=prior_sim, aes(x=sim_time))+
  geom_abline(data=prior_sim, aes(intercept = sample_c0, slope=sample_kr, group=sim_time), colour="blue")+
  xlim(c(0,200))+
  ylim(c(-100,300))+
  labs(x="time", y="predicted HMF concentration")+
  theme_bw()

# Code for Table S1 (Latex code):
\begin{table}[ht]
\centering
\begin{tabular}{lrr}
\toprule
Parameter & Rhat & n\_eff \\ 
\midrule
b\_Intercept & 1.00 & 4969 \\ 
b\_time & 1.00 & 4985 \\ 
sigma & 1.00 & 2931 \\ 
\bottomrule
\end{tabular}
\caption{Rhat and n\_eff for single-level Bayesian regression for HMF experiment 1; 8000 post-warmup MCMC samples.} 
\end{table}

# Code for Figure S3:

HMF_1 <- subset(HMF, serial==1, select=time:conc)

HMFfit1 <- 
  brm(data = HMF_1, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0,10)", class = "Intercept"),
                set_prior("normal(0.14, 0.1)", class = "b"),
                set_prior("cauchy(0,10)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      file=here("fits", "HMFfit1"))
postHMF1 <- posterior_samples(HMFfit1) %>% select(-lp__)

# Code for Figure S4, analyzes the repetitions one by one

HMF_2 <- subset(HMF, serial==2, select=time:conc)

HMFfit2 <- 
  brm(data = HMF_2, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      seed=15, file=here("fits", "HMFfit2"))
post2 <- posterior_samples(HMFfit2)

HMF_3 <- subset(HMF, serial==3, select=time:conc)

HMFfit3 <- 
  brm(data = HMF_3, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      seed=15, file=here("fits", "HMFfit3"))
post3 <- posterior_samples(HMFfit3)

HMF_4 <- subset(HMF, serial==4, select=time:conc)

HMFfit4 <- 
  brm(data = HMF_4, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      seed=15, file=here("fits", "HMFfit4"))
post4 <- posterior_samples(HMFfit4)

HMF_5 <- subset(HMF, serial==5, select=time:conc)

HMFfit5 <- 
  brm(data = HMF_5, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      seed=15, file=here("fits", "HMFfit5"))
post5 <- posterior_samples(HMFfit5)

HMF_6 <- subset(HMF, serial==6, select=time:conc)

HMFfit6 <- 
  brm(data = HMF_6, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      seed=15, file=here("fits", "HMFfit6"))
post6 <- posterior_samples(HMFfit6)

# building df with no-pooling results, ic_np=intercepts no pooling, sl_np=slope no pooling

serial <- c(1,2,3,4,5,6)
ic_np <- c(fixef(HMFfit1)[1,1], fixef(HMFfit2)[1,1],fixef(HMFfit3)[1,1],fixef(HMFfit4)[1,1],fixef(HMFfit5)[1,1],fixef(HMFfit6)[1,1])
sl_np <- c(fixef(HMFfit1)[2,1],fixef(HMFfit2)[2,1],fixef(HMFfit3)[2,1],fixef(HMFfit4)[2,1],fixef(HMFfit5)[2,1],fixef(HMFfit6)[2,1])
df_no_pooling <- data.frame(serial, ic_np,sl_np) %>% mutate(serial=str_c("Experiment ", serial))


HMF_plot <- HMF %>% mutate(serial=str_c("Experiment ", serial))
df_no_pooling %>% 
  ggplot()+
  geom_point(data=HMF_plot,aes(x=time, y=conc),color="red", size=2)+
  geom_abline(data=df_no_pooling, aes(intercept=ic_np, slope=sl_np), size=1.2)+
  facet_wrap(~serial)+
  labs(x="time, days", y=expression(paste("HMF concentration, ", mu,"mol/L")))+
  theme(strip.background = element_rect(color="black", fill="lightblue", size=1.5, linetype="solid")   )

# Code for Figure 4:

#first extract a file 'fit' with posterior_samples, then select parameters of interest in the dataframe to plot
cor_HMF1 <- dplyr::select(post1,b_Intercept:sigma)
# change the names of the columns to be displayed in the panels
cor_HMF1 <- setNames(cor_HMF1, c(expression(paste(italic("c")[0])), 
                                 expression(paste(italic("k")[r])), expression(sigma[e])))

# use ggally for a pairs plot
(cor_HMF1_plot <-cor_HMF1  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")),
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))+
    theme(axis.text.x = element_text(angle=70, hjust=1))
  
)

# Code for Figure 5:

p_HMF <-
  bind_rows(
    post1,
    post2,
    post3,
    post4,
    post5,
    post6
  )
iter <- 8000

p_HMF <- 
  p_HMF %>% 
  mutate(serial = rep(c("experiment 1","experiment 2","experiment 3","experiment 4","experiment 5","experiment 6"),each = iter)) 

plot_c0 <- p_HMF %>% 
  ggplot(aes(x = b_Intercept, y = serial)) +
  geom_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(italic("c")[0]), y="")


plot_kr <- p_HMF %>% 
  ggplot(aes(x = b_time, y = serial)) +
  stat_halfeyeh(fill = "lightblue", point_interval = mean_qi, .width = .95) +
  labs(x=expression(italic("k")[r]), y="")

plot_c0_kr=plot_c0 + plot_kr
plot_c0_kr[[2]]=plot_c0_kr[[2]]+theme(axis.text.y=element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.title.y = element_blank() )
plot_c0_kr

# Code for Figure S5:

HMFfit_all <- 
  brm(data = HMF, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000,  control = list(adapt_delta = 0.9),
      seed=15, file=here("fits", "HMFfit_all"))
HMFfit_all_post <- posterior_samples(HMFfit_all) %>% select(-lp__)
bayesplot::mcmc_trace(HMFfit_all_post)

# Code for Table S2 (Latex code):

\begin{table}[ht]
\centering
\begin{tabular}{lrr}
\toprule
Parameter & Rhat & n\_eff \\ 
\midrule
b\_Intercept & 1.00 & 7130 \\ 
b\_time & 1.00 & 7386 \\ 
sigma & 1.00 & 6564 \\ 
\bottomrule
\end{tabular}
\caption{Rhat and n\_eff for single-level Bayesian regression for the HMF data with complete pooling; 8000 post-warmup MCMC samples} 
\end{table}

# Code for Figure S6:

HMF_cor2 <- dplyr::select(HMFfit_all_post,b_Intercept:sigma)

HMF_cor2 <- setNames(HMF_cor2, c(expression(paste("c"[0])), expression(paste("k"[r])), expression(sigma[e])))

(HMF_corplot2 <-HMF_cor2 %>%  ggpairs(diag=list(continuous="densityDiag"),
                                      mapping=aes(fill="red"),
                                      upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")),
                                      labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))+
    theme(axis.text.x = element_text(angle=70, hjust=1))
)

# Code for Figure S7:

HMF_resid <- HMF %>%
  add_residual_draws(HMFfit_all) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval()+
  labs(subtitle = "A")

HMF_QQ <- HMF %>%
  add_residual_draws(HMFfit_all) %>%
  mean_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()+
  labs(subtitle = "B")

resid1 = resid(HMFfit_all)[, "Estimate"]
lagx1 <- resid1[-1]
lagy1 <- resid1[-length((resid1))]
lag1 <- data.frame(lagx1,lagy1)
lag_HMF <- ggplot(lag1, aes(x=lagx1, y=lagy1))+
  geom_point()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(subtitle = "C")

HMF_ppc <- pp_check(HMFfit_all, nsamples = 100)+labs(x="HMF values", subtitle = "D")
HMF_resid + HMF_QQ + lag_HMF + HMF_ppc 

# Code for complete pooling regression:

HMFfit_all <- 
  brm(data = HMF, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.05, 0.5)", class = "b"),
                set_prior("cauchy(0,10)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000,  control = list(adapt_delta = 0.9),
      seed=15, file=here("fits", "HMFfit_all"))
HMFfit_all_post <- posterior_samples(HMFfit_all)

# Code for Table 1:

HMFfit_pooled_summary <- summary(HMFfit_all)

HMFfit_pooled_summary <- rbind(data.frame(HMFfit_pooled_summary$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS),  data.frame(HMFfit_pooled_summary$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

rownames(HMFfit_pooled_summary) <- c("$c_0 in  \\mu \\text {mol L}^{-1}$", "$k_r in  \\mu \\text {mol L}{^-1} \\text {day}^{-1})$", "$\\sigma_{e}$")

colnames(HMFfit_pooled_summary) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  HMFfit_pooled_summary,
  placement = "H",
  align = c("c", "c", "c","c"),
  caption = "(ref:HMF-summary-pooled)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,3,2,2)
  ),
  escape = FALSE
)

# Code for Figure 6, fit resulting from completely pooled data:

newavg <- data.frame(time = seq(from = 0, to = 310, by = 10))

# credible interval
fitavg1 <- cbind(newavg, fitted(HMFfit_all, newdata = newavg, re_formula = NA)[,-2])
names(fitavg1) <- c("time", "conc", "lower", "upper")

# prediction interval
predavg1 <- cbind(newavg, predict(HMFfit_all, newdata = newavg, re_formula = NA)[,-2])
names(predavg1) <- c("time", "conc", "lower", "upper")

(HMF_plot_all <- ggplot(HMF,aes(x=time, y=conc)) +
    geom_point(color="red")+  # data
    geom_line(data=fitavg1, (aes(y=conc)), size=1.5)+ # regression line added
    geom_ribbon(data = fitavg1, aes(ymin = lower, ymax = upper), fill = "blue", alpha =0.5)+ # CI added
    geom_ribbon(data = predavg1, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha =0.5)+ # PI added
    labs(x="time in days", y=expression(paste("HMF in  ", mu,"mol L"^-1)))
)


# Code for partial pooling regression:

HMFfit_multi <- 
  brm(data = HMF, family = gaussian,
      formula = conc ~ 1 + time +(1+time|serial),
      prior = c(set_prior("normal(0, 10)", class = "Intercept"),
                set_prior("normal(0.1, 0.05)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000,  control = list(adapt_delta = 0.9),
      seed=15, file=here("fits", "HMFfit_multi"))
HMFfit_multi_post <- posterior_samples(HMFfit_multi)

# Code for Figure S8:

bayesplot::mcmc_trace(HMFfit_multi_post)+
  theme(axis.text.x = element_text(angle=70, hjust=1))

# Code for Table S3 (latex code):

\begin{table}[ht]
\centering
\begin{tabular}{lrr}
\toprule
Parameter & Rhat & n\_eff \\ 
\midrule
b\_Intercept & 1.00 & 2487 \\ 
b\_time & 1.00 & 2617 \\ 
sd\_serial\_\_Intercept & 1.00 & 2449 \\ 
sd\_serial\_\_time & 1.00 & 2178 \\ 
cor\_serial\_\_Intercept\_\_time & 1.00 & 2616 \\ 
sigma & 1.00 & 4487 \\ 
r\_serial[1,Intercept] & 1.00 & 3296 \\ 
r\_serial[2,Intercept] & 1.00 & 3143 \\ 
r\_serial[3,Intercept] & 1.00 & 3641 \\ 
r\_serial[4,Intercept] & 1.00 & 2953 \\ 
r\_serial[5,Intercept] & 1.00 & 3100 \\ 
r\_serial[6,Intercept] & 1.00 & 2744 \\ 
r\_serial[1,time] & 1.00 & 3473 \\ 
r\_serial[2,time] & 1.00 & 3319 \\ 
r\_serial[3,time] & 1.00 & 4131 \\ 
r\_serial[4,time] & 1.00 & 3226 \\ 
r\_serial[5,time] & 1.00 & 3100 \\ 
r\_serial[6,time] & 1.00 & 2837 \\ 
\bottomrule
\end{tabular}
\caption{Rhat and n\_eff for multilevel Bayesian regression for the HMF data with partial pooling (8000 post-warmup MCMC samples); serial refers to the experiment number} 
\end{table}

# Code for Figure S9:

HMFfit_cor2 <- dplyr::select(HMFfit_multi_post,b_Intercept:sigma)
HMFfit_cor2 <- setNames(HMFfit_cor2, c(expression(paste("c"[0])), expression(k[r]), expression(sigma[u]),expression(sigma[v]),expression(rho[uv]), expression(sigma[e])))


(HMFfit_corplot2 <-HMFfit_cor2  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))+
    theme(axis.text.x = element_text(angle=70, hjust=1))
) 

# Code for Figure 7:

re_plot1 <- mcmc_plot(HMFfit_multi, pars=c("^r_serial"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("experiment 1","experiment 2","experiment 3","experiment 4","experiment 5","experiment 6"), limits=c("r_serial[1,Intercept]","r_serial[2,Intercept]","r_serial[3,Intercept]","r_serial[4,Intercept]","r_serial[5,Intercept]", "r_serial[6,Intercept]")) +
  labs(x=expression(paste(Delta,italic("(c")[0],") in ",mu,"mol L"^-1)))

re_plot2 <- mcmc_plot(HMFfit_multi, pars=c("^r_serial"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("experiment 1","experiment 2","experiment 3","experiment 4","experiment 5","experiment 6"),limits=c("r_serial[1,time]","r_serial[2,time]","r_serial[3,time]","r_serial[4,time]","r_serial[5,time]", "r_serial[6,time]")) + labs(x=expression(paste(Delta,italic("(k"[r]),") in ",mu,"mol L"^-1,"day"^-1)))+
  xlim(-0.1,0.1)

re_plot_1_2 = re_plot1 + re_plot2
re_plot_1_2[[2]]=re_plot_1_2[[2]]+theme(axis.text.y = element_blank(),
                                        axis.ticks.y= element_blank(),
                                        axis.title.y = element_blank() )
re_plot_1_2

# Code for Table 2:

HMFfit_multi_cor <- summary(HMFfit_multi)

HMFfit_multi_cor2 <- rbind(data.frame(HMFfit_multi_cor$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(do.call(rbind,HMFfit_multi_cor$random))%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(HMFfit_multi_cor$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))


rownames(HMFfit_multi_cor2) <- c("$c_0  ( \\mu \\text {mol L}^{-1})$", "$k_r ( \\mu \\text {mol L}^{-1} \\text {day}^{-1})$", "$\\sigma_u$", "$\\sigma_v$", "$\\rho_{uv}$","$\\sigma_{e}$")

colnames(HMFfit_multi_cor2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  HMFfit_multi_cor2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:HMF-summary-multi)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,3,2,2)
  ),
  escape = FALSE
)

# Code for Figure 8:

newavg <- data.frame(time = seq(from = 0, to = 310, by = 10))

# credible interval
fitavg2 <- cbind(newavg, fitted(HMFfit_multi, newdata = newavg, re_formula = NA)[,-2])
names(fitavg2) <- c("time", "conc", "lower", "upper")

# prediction interval
predavg2 <- cbind(newavg, predict(HMFfit_multi, newdata = newavg, re_formula = NA)[,-2])
names(predavg2) <- c("time", "conc", "lower", "upper")

# CI from complete pooling, adding the CI interval from partial pooling:
HMF_plot_CI_combined <- ggplot(HMF,aes(x=time, y=conc)) +
  geom_point(color="red")+  # data
  geom_line(data=fitavg1, (aes(y=conc)), size=1.5, lty=2)+ # regression line completely pooled
  geom_line(data=fitavg2, (aes(y=conc)), lty=1, size=1.5)+ #regression line partially pooled
  geom_ribbon(data = fitavg1, aes(ymin = lower, ymax = upper), fill = "blue", alpha =0.5)+ # CI from completely pooled added
  geom_ribbon(data = fitavg2, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha =0.7)+ # CI from partially pooled added
  labs(x="time in days", y=expression(paste("HMF in  ", mu,"mol L"^-1)), subtitle = "A")  

# adding the PI interval
HMF_plot_PI_combined <- ggplot(HMF,aes(x=time, y=conc)) +
  geom_point(color="red")+  # data
  geom_line(data=fitavg1, (aes(y=conc)), size=1.5, lty=2)+ # regression line completely pooled
  geom_line(data=fitavg2, (aes(y=conc)), lty=1, size=1.5)+ #regression line partially pooled
  geom_ribbon(data = predavg1, aes(ymin = lower, ymax = upper), fill = "blue", alpha =0.5)+ # PI from completely pooled added
  geom_ribbon(data = predavg2, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha =0.7)+ # PI from partially pooled added
  labs(x="time in days", y=expression(paste("HMF in  ", mu,"mol L"^-1)), subtitle = "B")  

HMF_plot_CI_combined + HMF_plot_PI_combined 

# Code for Figure 9:

# building df with no-pooling results, ic_np=intercepts no pooling, sl_np=slope no pooling
serial <- c(1,2,3,4,5,6)
ic_np <- c(fixef(HMFfit1)[1,1], fixef(HMFfit2)[1,1],fixef(HMFfit3)[1,1],fixef(HMFfit4)[1,1],fixef(HMFfit5)[1,1],fixef(HMFfit6)[1,1])
sl_np <- c(fixef(HMFfit1)[2,1],fixef(HMFfit2)[2,1],fixef(HMFfit3)[2,1],fixef(HMFfit4)[2,1],fixef(HMFfit5)[2,1],fixef(HMFfit6)[2,1])
df_no_pooling <- data.frame(serial, ic_np,sl_np)
df_no_pooling$serial <- as.factor(df_no_pooling$serial)

# building df with completely pooled results, cp=completely pooled
ic_cp <- c(fixef(HMFfit_all)[1,1], fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1])
sl_cp <- c(fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1])
serial <- c(1,2,3,4,5,6)
df_cp <- data.frame(serial,ic_cp, sl_cp)
df_cp$serial <- as.factor(df_cp$serial)

# building df with partially pooled results population level, pp=partially pooled
ic_pp <- c(fixef(HMFfit_multi)[1,1], fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1])
sl_pp <- c(fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1])
serial <- c(1,2,3,4,5,6)
df_pp <- data.frame(serial,ic_pp, sl_pp)
df_pp$serial <- as.factor(df_pp$serial)

#building df with partially pooled results group level, ppg=partially pooled group level
# the regression coefficients for the group level can be extracted using the following command, it gives the coefficients directly so not as deviations

df_ppg <- coef(HMFfit_multi)$serial[, ,] %>% as_tibble() %>% select(Estimate.Intercept, Estimate.time) %>% bind_cols(serial) %>% rename(serial="...3") 
df_ppg$serial <- as.factor(df_ppg$serial)

# plot showing the fits for serial 1:
ggplot(data=HMF_1, aes(x=time, y=conc))+
  geom_point(color="red", size=2)+
  geom_abline(data=df_no_pooling %>% filter(serial %in% c(1)), aes(intercept=ic_np, slope=sl_np),lty=1, size=1.2)+ #black: no-pooling
  #  geom_abline(data=df_cp %>% filter(serial %in% c(1)), aes(intercept=ic_cp, slope=sl_cp), lty=2, color="blue", size=1.2)+ #blue: complete pooling
  geom_abline(data=df_pp %>% filter(serial %in% c(1)), aes(intercept=ic_pp, slope=sl_pp), lty=1, color="red", size=1.2)+ # red: partial at population level (grand mean)
  geom_abline(data=df_ppg %>% filter(serial %in% c(1)), aes(intercept=Estimate.Intercept, slope=Estimate.time), lty=4, color="darkgreen", size=1.2)+
  labs(x="time in days", y=expression(paste("HMF in  ", mu,"mol L"^-1)))

# Code for Table 3:

HMFfit_all <- add_criterion(HMFfit_all, c("loo"), file=here("fits", "HMFfit_all"))
HMFfit_multi <- add_criterion(HMFfit_multi, c("loo"), file=here("fits", "HMFfit_multi"))

loo_all_HMF <- loo(HMFfit_all)
loo_multi_HMF <- loo(HMFfit_multi)

loo_result_HMF <- loo_compare(loo_all_HMF,loo_multi_HMF)

elpd_diff_HMF <- c(loo_result_HMF[1,1], loo_result_HMF[2,1])
se_diff_HMF <- c(loo_result_HMF[1,2], loo_result_HMF[2,2])
loo_df_HMF <- data.frame(elpd_diff_HMF, se_diff_HMF)
colnames(loo_df_HMF) <- c("elpd-difference","SE" )
rownames(loo_df_HMF) <- c("partial pooling", "complete pooling")

apa_table(
  loo_df_HMF,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:HMF-loo)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Code for Figure S10
# building df with no-pooling results, ic_np=intercepts no pooling, sl_np=slope no pooling

# building df with completely pooled results, cp=completely pooled
ic_cp <- c(fixef(HMFfit_all)[1,1], fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1],fixef(HMFfit_all)[1,1])
sl_cp <- c(fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1],fixef(HMFfit_all)[2,1])
serial <- c(1,2,3,4,5,6)
df_cp <- data.frame(serial,ic_cp, sl_cp)
df_cp$serial <- as.factor(df_cp$serial)

# building df with partially pooled results population level, pp=partially pooled
ic_pp <- c(fixef(HMFfit_multi)[1,1], fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1],fixef(HMFfit_multi)[1,1])
sl_pp <- c(fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1],fixef(HMFfit_multi)[2,1])
serial <- c(1,2,3,4,5,6)
df_pp <- data.frame(serial,ic_pp, sl_pp) %>% mutate(serial=str_c("Experiment ", serial))

#building df with partially pooled results group level, ppg=partially pooled group level
# the regression coefficients for the group level can be extracted using the following command, it gives the coefficients directly so not as deviations

df_ppg <- coef(HMFfit_multi)$serial[, ,] %>% as_tibble() %>% select(Estimate.Intercept, Estimate.time) %>% bind_cols(serial) %>% rename(serial="...3") %>% mutate(serial=str_c("Experiment ", serial))

#df_ppg$serial <- as.factor(df_ppg$serial)

ggplot(data=HMF_plot, aes(x=time, y=conc))+
  geom_point(color="red", size=2)+
  geom_abline(data=df_no_pooling, aes(intercept=ic_np, slope=sl_np), size=1.3)+
  #   geom_abline(data=df_cp, aes(intercept=ic_cp, slope=sl_cp), lty=2, size=1.2, color="blue")+
  geom_abline(data=df_pp, aes(intercept=ic_pp, slope=sl_pp), size=1.2, color="red")+ # grand mean
  geom_abline(data=df_ppg, aes(intercept=Estimate.Intercept, slope=Estimate.time), lty=4, size=1.3, color="darkgreen")+
  facet_wrap(~serial)+
  labs(x="time in days", y=expression(paste("HMF in  ",mu,"mol L"^-1)))

# THIAMIN CASE STUDY 2

# Code for thiamin data + Figure S11:

# the following data contain thiamin data but excluding the extra storage data at 35 and 50 C; time in days AND min
thiamin <- read.csv(here("data", "thiamin_horak_fink.csv"),header=TRUE, sep=";")

#The following files contain data with time in min:
thiamin_120_min <- subset(thiamin, temp ==120, select=time:conc)
thiamin_130_min <- subset(thiamin, temp ==130, select=time:conc)
thiamin_140_min <- subset(thiamin, temp ==140, select=time:conc)
thiamin_150_min <- subset(thiamin, temp ==150, select=time:conc)

# all data, time in DAYS and extra storage data at 35 and 50 added (154 data):
thiamin_all_days <- read.csv(here("data", "thiamin_all_days.csv"),header=TRUE, sep=";")

thiamin_35 <- subset(thiamin_all_days, temp ==35, select=time:conc)
thiamin_50 <- subset(thiamin_all_days, temp ==50, select=time:conc)
thiamin_72 <- subset(thiamin_all_days, temp ==72, select=time:conc)
thiamin_85 <- subset(thiamin_all_days, temp ==85, select=time:conc)
thiamin_120 <- subset(thiamin_all_days, temp ==120, select=time:conc)

thiamin_all_days %>% mutate(Temp=str_c("T=",temp, " °C")) %>% ggplot(aes(x=time, y=conc, group=temp))+
  geom_point(aes(x=time,y=conc))+
  facet_wrap(~Temp, ncol=4, scales="free")+
  labs(x = "time (days)", y = "thiamin (mg/L)")+
  theme(axis.text.x = element_text(angle=70, hjust=1))

# Code for individual regression of the thiamin data:

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0", lb=0),
           prior(normal(0.018,0.002), nlpar="k", lb=0),
           prior(normal(2,0.5), nlpar="nt", lb=0),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_35 <- brm(formula=nlform, data=thiamin_35, family = gaussian(), 
                        prior = nlprior, iter = 16000, warmup=8000, 
                        control = list(adapt_delta = 0.999),
                        file=here("fits", "thiamin_model_35"))

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0", lb=0),
           prior(normal(0.06,0.005), nlpar="k", lb=0),
           prior(normal(2,0.5), nlpar="nt"),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_50 <- brm(formula=nlform, data=thiamin_50, family = gaussian(), 
                        prior = nlprior, iter = 16000, warmup=8000, 
                        control = list(adapt_delta = 0.999),
                        file=here("fits","thiamin_model_50"))

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0", lb=0),
           prior(normal(0.9,0.1), nlpar="k", lb=0),
           prior(normal(2,0.5), nlpar="nt", lb=0),
           prior(cauchy(0,10), class="sigma")
)           


thiamin_model_72 <- brm(formula=nlform, data=thiamin_72, family = gaussian(), 
                        prior = nlprior, iter = 16000, warmup=8000, 
                        control = list(adapt_delta = 0.999),
                        file=here("fits", "thiamin_model_72"))


nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0", lb=0),
           prior(normal(2.6,0.5), nlpar="k", lb=0),
           prior(normal(2,0.5), nlpar="nt",lb=0),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_85 <- brm(formula=nlform, data=thiamin_85, family = gaussian(), 
                        prior = nlprior, iter = 16000, warmup=8000, 
                        control = list(adapt_delta = 0.999),
                        file=here("fits", "thiamin_model_85"))

# Take care: the above results are for time in days, below is for time in min. Time in mins gives more stable results, used to get the best estimate of the order nt for this range of temperatures

#thiamin_120 gives however a better fit when days are used instead of mins:
nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.32,0.1), nlpar = "c0",lb=0),
           prior(normal(300,30), nlpar="k",lb=0),
           prior(normal(2,0.5), nlpar="nt",lb=0),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_120 <- brm(formula=nlform, data=thiamin_120, family = gaussian(), 
                         prior = nlprior, iter = 16000, warmup=8000, 
                         control = list(adapt_delta = 0.999),
                         file=here("fits", "thiamin_model_120"))


nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0",lb=0),
           prior(normal(0.05,0.01), nlpar="k",lb=0),
           prior(normal(2,0.1), nlpar="nt",lb=0),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_130_min <- brm(formula=nlform, data=thiamin_130_min, family = gaussian(), 
                             prior = nlprior, iter = 16000, warmup=8000, 
                             control = list(adapt_delta = 0.999),
                             file=here("fits", "thiamin_model_130_min"))


nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0",lb=0),
           prior(normal(0.16,0.01), nlpar="k",lb=0),
           prior(normal(2,0.1), nlpar="nt",lb=0),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_140_min <- brm(formula=nlform, data=thiamin_140_min, family = gaussian(), 
                             prior = nlprior, iter = 16000, warmup=8000, 
                             control = list(adapt_delta = 0.999),
                             file=here("fits", "thiamin_model_140_min"))


nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(0.4,0.1), nlpar = "c0",lb=0),
           prior(normal(0.5,0.05), nlpar="k",lb=0),
           prior(normal(2,0.1), nlpar="nt",lb=0),
           prior(cauchy(0,10), class="sigma")
)           

thiamin_model_150_min <- brm(formula=nlform, data=thiamin_150_min, family = gaussian(), 
                             prior = nlprior, iter = 16000, warmup=8000, 
                             control = list(adapt_delta = 0.999),
                             file=here("fits", "thiamin_model_150_min"))

#The next code calculates individual fits from the individual regressions and puts them in a dataframe fits_unpooled; since time varies between days and mins it needs to be done this way
newvary1 <- expand.grid(time=seq(from = 0, to =210, by=5),temp=c(35))
newvary2 <- expand.grid(time=seq(from = 0, to =140, by=4),temp=c(50))
newvary3 <- expand.grid(time=seq(from = 0, to =8, by=0.2),temp=c(72))
newvary4 <- expand.grid(time=seq(from = 0, to =6, by=0.2),temp=c(85))
newvary5 <- expand.grid(time=seq(from = 0, to =0.06, by=0.002),temp=c(120))
newvary6 <- expand.grid(time=seq(from = 0, to =0.06, by=0.002),temp=c(130))
newvary7 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(140))
newvary8 <- expand.grid(time=seq(from = 0, to =0.012, by=0.0004),temp=c(150))
#newvary <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6, newvary7, newvary8)

c0_35 <- rep(c(fixef(thiamin_model_35)[1,1]), times=nrow(newvary1))
nt_35<- rep(c(fixef(thiamin_model_35)[3,1]), times=nrow(newvary1))
kr_35<- rep(c(fixef(thiamin_model_35)[2,1]), times=nrow(newvary1))
fit_35 <- cbind(newvary1,c0_35,nt_35,kr_35) %>% mutate(fit=(c0_35^(1-nt_35)+(nt_35-1)*kr_35*time)^(1/(1-nt_35))) %>% select(-c0_35,-nt_35,-kr_35)

c0_50 <- rep(c(fixef(thiamin_model_50)[1,1]), times=nrow(newvary2))
nt_50<- rep(c(fixef(thiamin_model_50)[3,1]), times=nrow(newvary2))
kr_50<- rep(c(fixef(thiamin_model_50)[2,1]), times=nrow(newvary2))
fit_50 <- cbind(newvary2,c0_50,nt_50,kr_50) %>% mutate(fit=(c0_50^(1-nt_50)+(nt_50-1)*kr_50*time)^(1/(1-nt_50))) %>% select(-c0_50,-nt_50,-kr_50)

c0_72 <- rep(c(fixef(thiamin_model_72)[1,1]), times=nrow(newvary3))
nt_72<- rep(c(fixef(thiamin_model_72)[3,1]), times=nrow(newvary3))
kr_72<- rep(c(fixef(thiamin_model_72)[2,1]), times=nrow(newvary3))
fit_72 <- cbind(newvary3,c0_72,nt_72,kr_72) %>% mutate(fit=(c0_72^(1-nt_72)+(nt_72-1)*kr_72*time)^(1/(1-nt_72))) %>% select(-c0_72,-nt_72,-kr_72)

c0_85 <- rep(c(fixef(thiamin_model_85)[1,1]), times=nrow(newvary4))
nt_85<- rep(c(fixef(thiamin_model_85)[3,1]), times=nrow(newvary4))
kr_85<- rep(c(fixef(thiamin_model_85)[2,1]), times=nrow(newvary4))
fit_85 <- cbind(newvary4,c0_85,nt_85,kr_85) %>% mutate(fit=(c0_85^(1-nt_85)+(nt_85-1)*kr_85*time)^(1/(1-nt_85))) %>% select(-c0_85,-nt_85,-kr_85)

c0_120 <- rep(c(fixef(thiamin_model_120)[1,1]), times=nrow(newvary5))
nt_120<- rep(c(fixef(thiamin_model_120)[3,1]), times=nrow(newvary5))
kr_120<- rep(c(fixef(thiamin_model_120)[2,1]), times=nrow(newvary5))
fit_120 <- cbind(newvary5,c0_120,nt_120,kr_120) %>% mutate(fit=(c0_120^(1-nt_120)+(nt_120-1)*kr_120*time)^(1/(1-nt_120))) %>% select(-c0_120,-nt_120,-kr_120)

c0_130 <- rep(c(fixef(thiamin_model_130_min)[1,1]), times=nrow(newvary6))
nt_130<- rep(c(fixef(thiamin_model_130_min)[3,1]), times=nrow(newvary6))
kr_130<- rep(c(fixef(thiamin_model_130_min)[2,1]*1440), times=nrow(newvary6))
fit_130 <- cbind(newvary6,c0_130,nt_130,kr_130) %>% mutate(fit=(c0_130^(1-nt_130)+(nt_130-1)*kr_130*time)^(1/(1-nt_130))) %>% select(-c0_130,-nt_130,-kr_130)

c0_140 <- rep(c(fixef(thiamin_model_140_min)[1,1]), times=nrow(newvary7))
nt_140<- rep(c(fixef(thiamin_model_140_min)[3,1]), times=nrow(newvary7))
kr_140<- rep(c(fixef(thiamin_model_140_min)[2,1]*1440), times=nrow(newvary7))
fit_140 <- cbind(newvary7,c0_140,nt_140,kr_140) %>% mutate(fit=(c0_140^(1-nt_140)+(nt_140-1)*kr_140*time)^(1/(1-nt_140))) %>% select(-c0_140,-nt_140,-kr_140)

c0_150 <- rep(c(fixef(thiamin_model_150_min)[1,1]), times=nrow(newvary8))
nt_150<- rep(c(fixef(thiamin_model_150_min)[3,1]), times=nrow(newvary8))
kr_150<- rep(c(fixef(thiamin_model_150_min)[2,1]*1440), times=nrow(newvary8))
fit_150 <- cbind(newvary8,c0_150,nt_150,kr_150) %>% mutate(fit=(c0_150^(1-nt_150)+(nt_150-1)*kr_150*time)^(1/(1-nt_150))) %>% select(-c0_150,-nt_150,-kr_150)

fits_unpooled <- rbind(fit_35,fit_50,fit_72,fit_85,fit_120,fit_130,fit_140,fit_150)

# Code for Figure 10:

post1 <- posterior_samples(thiamin_model_35)
post2 <- posterior_samples(thiamin_model_50)
post3 <- posterior_samples(thiamin_model_72)
post4 <- posterior_samples(thiamin_model_85)
post5 <- posterior_samples(thiamin_model_120)
post6 <- posterior_samples(thiamin_model_130_min)
post7 <- posterior_samples(thiamin_model_140_min)
post8 <- posterior_samples(thiamin_model_150_min)


p_thiamin_nt <-
  bind_rows(
    post1,
    post2,
    post3,
    post4,
    post5,
    post6,
    post7,
    post8
  )
iter <- 32000

p_thiamin_nt <- 
  p_thiamin_nt %>% 
  mutate(temperature = rep(c("35 °C","50 °C","72 °C","85 °C", "120 °C", "130 °C", "140 °C", "150 °C"),each = iter)) %>%  mutate(temperature=fct_relevel(temperature, "35 °C","50 °C","72 °C","85 °C", "120 °C", "130 °C", "140 °C", "150 °C"))
p_thiamin_nt %>% 
  ggplot(aes(x = b_nt_Intercept, y = temperature)) +
  stat_halfeye(fill = "green4", alpha=0.5,
               point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste(n[t], "-value")), y=expression(paste("n"[t]," density for experiment at T =")))


# Code for global regression of the thiamin data with the nth-order model, times are in days for all temperatures,  Tb = 367.3 K

# equation for the model:
nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*kb*exp(Ea*(1-367.3/(temp+273)))*time)^(1/(1-nt)), 
           c0 + nt+kb+Ea ~1, 
           nl=TRUE) 

# set priors for the parameters:
nlprior<-c(prior(normal(0.37,0.05), nlpar = "c0", lb=0),
           prior(normal(9,1), nlpar="kb", lb=0),
           prior(normal(2,0.5), nlpar="nt", lb=0),
           prior(normal(32.6, 1), nlpar="Ea", lb=0),
           prior(cauchy(0,25), class="sigma")
)

thiamin_all_model_nth<-brm(formula=nlform, data=thiamin_all_days, family = gaussian(), prior = nlprior, warmup=4000, iter=8000, control = list(adapt_delta = 0.9, max_treedepth=12), file=here("fits", "thiamin_all_model_nth"))

# Code for Figure S12:

preds_thiamin <- broom.mixed::augment(thiamin_all_model_nth)
preds_thiamin <- preds_thiamin %>% mutate(temp=str_c("T=", temp,"°C"))
thiamin_all_days <- thiamin_all_days %>% mutate(temp=str_c("T=", temp,"°C"))
ggplot(thiamin_all_days, aes(time, conc)) +
  geom_point(col = 'red') +
  geom_line(aes(time, .fitted), preds_thiamin) +
  facet_wrap(~temp, ncol=4, scales = "free")+
  ylab(expression(paste('thiamin concentration in mg L'^-1))) +
  xlab('time in days') +
  theme(axis.text.x = element_text(angle=70, hjust=1))

# Code for Figure S13:

thiamin_all_model_nth_post <- posterior_samples(thiamin_all_model_nth)
thiaminfit_cor2 <- dplyr::select(thiamin_all_model_nth_post,b_c0_Intercept:sigma)
thiaminfit_cor2 <- setNames(thiaminfit_cor2, c(expression(paste(c[0])), expression(n[t]),expression(k[ref]), expression(paste("E"[a],"RT"[ref])), expression(sigma[e])))


(thiaminfit_corplot2 <-thiaminfit_cor2  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))+
    theme(axis.text.x = element_text(angle=70, hjust=1))
) 

# Code for global regression with second order model, recalculated time in days, Tb = 367.3

# equation for the model:
nlform<-bf(conc ~ c0/(1+c0*10^kb*exp(Ea*(1-367.3/(temp+273.0)))*time), 
           c0 ~1+(1|temp),
           kb+Ea ~1, 
           nl=TRUE) 

# set priors for the parameters:
nlprior<-c(prior(normal(0.37,0.05), nlpar = "c0", lb=0),
           prior(normal(9,1), nlpar="kb"),
           prior(normal(32.6, 0.1), nlpar="Ea", lb=0),
           prior(cauchy(0,25), class="sigma")
)

thiamin_all_model_days2<-brm(formula=nlform, data=thiamin_all_days, family = gaussian(), prior = nlprior, warmup=4000, iter=8000, control = list(adapt_delta = 0.9, max_treedepth=12), file=here("fits", "thiamin_all_model_days2"))

thiamin_all_model_days2_post <- posterior_samples(thiamin_all_model_days2)

# Code for Table 4:

thiaminfit_multi_cor <- summary(thiamin_all_model_days2)

thiaminfit_multi_cor2 <- rbind(data.frame(thiaminfit_multi_cor$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(do.call(rbind,thiaminfit_multi_cor$random))%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(thiaminfit_multi_cor$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

# back transform Ea and kref:
thiaminfit_multi_cor2[3,] <- (thiaminfit_multi_cor2[3,]*8.314*367.3)/1000

thiaminfit_multi_cor2[2,] <- 10^thiaminfit_multi_cor2[2,]


rownames(thiaminfit_multi_cor2) <- c("$c_0 \\text { (mg/dm}^3)$", "$k_{ref} \\text { (dm}^3 \\text {mol}^{-1} \\text {day}^{-1})$", "$E_a \\text {  (kJ/mol)}$", "$\\sigma_{c_0}$","$\\sigma_{e}$")

colnames(thiaminfit_multi_cor2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  thiaminfit_multi_cor2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:thiamin-summary-multi)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# Code for Figure 11:

# combines the data with predictions using the random effects with re_formula = NULL
newvary1 <- expand.grid(time=seq(from = 0, to =210, by=5),temp=c(35))
newvary2 <- expand.grid(time=seq(from = 0, to =140, by=4),temp=c(50))
newvary3 <- expand.grid(time=seq(from = 0, to =8, by=0.2),temp=c(72))
newvary4 <- expand.grid(time=seq(from = 0, to =6, by=0.2),temp=c(85))
newvary5 <- expand.grid(time=seq(from = 0, to =0.06, by=0.002),temp=c(120))
newvary6 <- expand.grid(time=seq(from = 0, to =0.06, by=0.002),temp=c(130))
newvary7 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(140))
newvary8 <- expand.grid(time=seq(from = 0, to =0.012, by=0.0004),temp=c(150))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6, newvary7, newvary8)

thiamin_ind <- cbind(newvary, predict(thiamin_all_model_days2, newdata=newvary, re_formula = NA)[,-2]) # fit for the population level
thiamin_ind2 <- cbind(newvary, predict(thiamin_all_model_days2, newdata=newvary, re_formula = NULL)[,-2]) #fit for the group level

names(thiamin_ind) <- c("time", "temp", "conc", "lower", "upper")
thiamin_ind$temp=as.factor(thiamin_ind$temp)
names(thiamin_ind2) <- c("time", "temp", "conc", "lower", "upper")
thiamin_ind2$temp=as.factor(thiamin_ind2$temp)
# labeling the facets:
thiamin_all_days <- thiamin_all_days %>% mutate(temp=str_c("T=", temp," °C"))
thiamin_ind <- thiamin_ind %>% mutate(temp=str_c("T=", temp," °C"))
thiamin_ind2 <- thiamin_ind2 %>% mutate(temp=str_c("T=", temp," °C"))
fits_unpooled <- fits_unpooled %>% mutate(temp=str_c("T=", temp," °C"))

(ind_fit <-  ggplot(thiamin_all_days, aes(x=time, y=conc))+
    geom_point(colour = "#2c3e50",fill = "#2c3e50")+
    facet_wrap(~temp, ncol=4, scales="free_x")+
    geom_line(data = thiamin_ind, aes(y = conc), size = 1, colour="blue") +
    geom_line(data = thiamin_ind2, aes(y = conc), size = 1, colour="red", lty=2) + #fit for the group level
    geom_line(data=fits_unpooled, aes(y=fit))+ # fit for the individual level
    labs(x="time in days", y=expression(paste("thiamin in mg L"^-1)))
)

# Code for Table 5:

thiamin_all_model_days2 <- add_criterion(thiamin_all_model_days2, c("loo"), file=here("fits", "thiamin_all_model_days2"))
thiamin_all_model_nth <- add_criterion(thiamin_all_model_nth, c("loo"), file=here("fits", "thiamin_all_model_nth"))

loo_thiamin_nth <- loo(thiamin_all_model_nth)
loo_thiamin_multi <- loo(thiamin_all_model_days2)

loo_result_thiamin <- loo_compare(loo_thiamin_nth,loo_thiamin_multi)

elpd_diff_thiamin <- c(loo_result_thiamin[1,1], loo_result_thiamin[2,1])
se_diff_thiamin <- c(loo_result_thiamin[1,2], loo_result_thiamin[2,2])
loo_df_thiamin <- data.frame(elpd_diff_thiamin, se_diff_thiamin)
colnames(loo_df_thiamin) <- c("elpd-difference","SE" )
rownames(loo_df_thiamin) <- c("partial pooling", "complete pooling")

apa_table(
  loo_df_thiamin,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:thiamin-loo)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Code for Table S4:

p1 <- summary(thiamin_all_model_nth)
summary_p1 <- rbind(data.frame(p1$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(p1$spec_pars)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

# back transform Ea:
summary_p1[4,] <- summary_p1[4,]*8.314*373.5747/1000

rownames(summary_p1) <- c("$c_0 \\text { in mg L}^{-1}$", "$n_t \\text { (-)}$", "$k_{\\text {ref}} \\text { in (L mg}^{-1})^{n_t-1} \\text {day}^{-1}$","$E_a \\text{ in kJ mol}^{-1}$", "$\\sigma_e \\text { in mg L}^{-1}$")
colnames(summary_p1) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  summary_p1,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:thiamin-global-summary)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)
# Code for two-step regression of temperature effect on thiamin:

#determination of the rate constant for a second-order (so) reaction after linearization of the rate equation

thiamin_so_35_lr <- brm(data=thiamin_35, 
                        formula= 1/conc ~1+time, 
                        family = gaussian(), 
                        prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                  set_prior("normal(0.01, 0.01)", class = "b"),
                                  set_prior("cauchy(0,10)", class = "sigma")),
                        iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                        file=here("fits", "thiamin_so_35_lr"))

thiamin_so_50_lr <- brm(data=thiamin_50, 
                        formula= 1/conc ~1+time, 
                        family = gaussian(), 
                        prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                  set_prior("normal(0.04, 0.01)", class = "b"),
                                  set_prior("cauchy(0,10)", class = "sigma")),
                        iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                        file=here("fits", "thiamin_so_50_lr"))

thiamin_so_72_lr <- brm(data=thiamin_72, 
                        formula= 1/conc ~1+time, 
                        family = gaussian(), 
                        prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                  set_prior("normal(1, 0.1)", class = "b"),
                                  set_prior("cauchy(0,10)", class = "sigma")),
                        iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                        file=here("fits", "thiamin_so_72_lr"))

thiamin_so_85_lr <- brm(data=thiamin_85, 
                        formula=(1/conc) ~1+time, 
                        family = gaussian(), 
                        prior = c(set_prior("normal(2.5,0.001)", class = "Intercept"),
                                  set_prior("normal(5, 0.5)", class = "b"),
                                  set_prior("cauchy(0,25)", class = "sigma")),
                        iter = 2000, warmup=1000, inits=0, control = list(adapt_delta = 0.9),
                        file=here("fits", "thiamin_so_85_lr"))

#time in mins from here but also checked with days:
thiamin_so_120_lr <- brm(data=thiamin_120_min, 
                         formula= 1/conc ~1+time, 
                         family = gaussian(), 
                         prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                   set_prior("normal(0.06, 0.01)", class = "b"),
                                   set_prior("cauchy(0,10)", class = "sigma")),
                         iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                         file=here("fits", "thiamin_so_120_lr"))

thiamin_so_120_lr_days <- brm(data=thiamin_120, 
                              formula= 1/conc ~1+time, 
                              family = gaussian(), 
                              prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                        set_prior("normal(82, 5)", class = "b"),
                                        set_prior("cauchy(0,10)", class = "sigma")),
                              iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                              file=here("fits", "thiamin_so_120_lr_days"))

thiamin_so_130_lr <- brm(data=thiamin_130_min, 
                         formula= 1/conc ~1+time, 
                         family = gaussian(), 
                         prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                   set_prior("normal(0.14, 0.01)", class = "b"),
                                   set_prior("cauchy(0,10)", class = "sigma")),
                         iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                         file=here("fits", "thiamin_so_130_lr"))

thiamin_so_130_lr_days <- brm(data=thiamin_130, 
                              formula= 1/conc ~1+time, 
                              family = gaussian(), 
                              prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                        set_prior("normal(222, 10)", class = "b"),
                                        set_prior("cauchy(0,10)", class = "sigma")),
                              iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                              file=here("fits", "thiamin_so_130_lr_days"))

thiamin_so_140_lr <- brm(data=thiamin_140_min, 
                         formula= 1/conc ~1+time, 
                         family = gaussian(), 
                         prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                   set_prior("normal(0.24, 0.01)", class = "b"),
                                   set_prior("cauchy(0,10)", class = "sigma")),
                         iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                         file=here("fits", "thiamin_so_140_lr"))

thiamin_so_140_lr_days <- brm(data=thiamin_140, 
                              formula= 1/conc ~1+time, 
                              family = gaussian(), 
                              prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                        set_prior("normal(340, 10)", class = "b"),
                                        set_prior("cauchy(0,10)", class = "sigma")),
                              iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                              file=here("fits", "thiamin_so_140_lr_days"))

thiamin_so_150_lr <- brm(data=thiamin_150_min, 
                         formula= 1/conc ~1+time, 
                         family = gaussian(), 
                         prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                   set_prior("normal(0.54, 0.01)", class = "b"),
                                   set_prior("cauchy(0,10)", class = "sigma")),
                         iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                         file=here("fits", "thiamin_so_150_lr"))

thiamin_so_150_lr_days <- brm(data=thiamin_150, 
                              formula= 1/conc ~1+time, 
                              family = gaussian(), 
                              prior = c(set_prior("normal(2.5,1)", class = "Intercept"),
                                        set_prior("normal(650, 60)", class = "b"),
                                        set_prior("cauchy(0,10)", class = "sigma")),
                              iter = 4000, warmup=2000, control = list(adapt_delta = 0.9),
                              file=here("fits", "thiamin_so_150_lr_days"))

# Code for Figure 12:

Tinv <- c(1/(8.314*308),1/(8.314*323),1/(8.314*345),1/(8.314*358),1/(8.314*393),1/(8.314*403),1/(8.314*413),1/(8.314*423))
# data from 35-85 are in days, 120-150 are in minutes recalculated to days (multiplied by 24*60=1440), 
ln_k <- c(log(fixef(thiamin_so_35_lr)[2,1]),log(fixef(thiamin_so_50_lr)[2,1]),log(fixef(thiamin_so_72_lr)[2,1]), log(fixef(thiamin_so_85_lr)[2,1]),log(fixef(thiamin_so_120_lr_days)[2,1]), log(fixef(thiamin_so_130_lr_days)[2,1]),log(fixef(thiamin_so_140_lr_days)[2,1]),log(fixef(thiamin_so_150_lr_days)[2,1]))
df_linArrh <- data.frame(Tinv,ln_k)

linArrhfit <- 
  brm(data = df_linArrh, family = gaussian,
      formula = ln_k ~ 1 + Tinv,
      prior = c(set_prior("normal(25, 10)", class = "Intercept"),
                set_prior("normal(-100000, 10000)", class = "b"),
                set_prior("cauchy(0,10)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits", "linArrhfit"))

# calculate k0 from its ln, and recalculate Ea in kJ/mol:
linArrh_post <- posterior_samples(linArrhfit) %>% mutate(b_Tinv=-b_Tinv) %>% mutate(k0=exp(b_Intercept)) %>% mutate(Ea=b_Tinv/1000)

# range to plot:
T.seq <- data.frame(Tinv = seq(from = 1/(8.314*430), to = 1/(305*8.314), by = 0.00001/8.314))

#fitted is about mu:
muSummary <-
  fitted(linArrhfit, 
         newdata = T.seq) %>%
  as_tibble() %>%
  bind_cols(T.seq)

#predict is about future individual values:
pred.Arrh <-
  predict(linArrhfit,
          newdata = T.seq) %>%
  as_tibble() %>%
  bind_cols(T.seq)

#plot of fitted and predicted values
plot_Arrh <- df_linArrh %>%
  ggplot(aes(x = Tinv, y = ln_k)) +
  geom_ribbon(data = pred.Arrh, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = muSummary, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = muSummary, aes(y = Estimate)) +
  geom_point(color = "navyblue", shape = 16, size = 1.5, alpha = 2/3)+
  labs(x=expression(paste("1/R",italic("T"))), y=expression(paste("ln(",italic("k")[r],")")))

plot_Arrh

# Code for Figure 13:

thiamin_all_model_days2_post <- thiamin_all_model_days2_post %>% mutate(Ea=b_Ea_Intercept*8.314*367.3/1000) %>% mutate(lnk0=log(10^b_kb_Intercept/(exp(-Ea*1000/(8.314*367.3))))) 

Ea_density <- ggplot(data=linArrh_post, aes(x=Ea))+
  geom_density(fill="red", alpha=0.4)+ 
  geom_density(data=thiamin_all_model_days2_post, aes(x=Ea), fill="blue", alpha=0.3)+
  labs(x=expression(paste(italic("E")[a], " in kJ mol"^-1)), subtitle = "A")  

lnk0_density <- ggplot(data=linArrh_post, aes(x=b_Intercept))+
  geom_density(fill="red", alpha=0.4)+
  geom_density(data=thiamin_all_model_days2_post, aes(x=lnk0), fill="blue",alpha=0.3)+
  labs(x=expression(paste("ln",italic("k")[0])), subtitle = "B")

Ea_density + lnk0_density


# Code for Figure 14:

# this part calculates what the thiamin content would be for heating milk with c0=0.4 mg/L for 20 min at 120 C followed by storage at 20 C for 50, 75 and 100 days::

mu_120 <- data.frame(ct=0.4/(1+0.4*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/393))*20/(60*24)))%>% 
  mutate(ct2=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*10)) %>% 
  mutate(ct3=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*50))%>% 
  mutate(ct4=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*75))%>% 
  mutate(ct5=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*100))

plotA <- ggplot(mu_120, aes(x=ct))+
  geom_density(fill="turquoise", alpha=0.5)+
  geom_density(aes(x=ct2), fill="red", alpha=0.5)+
  geom_density(aes(x=ct3), fill="blue",alpha=0.5)+
  geom_density(aes(x=ct4), fill="yellow", alpha=0.5)+
  geom_vline(xintercept = 0.4, size=3)+
  geom_density(aes(x=ct5), fill="purple", alpha=0.5)+
  labs(x=expression(paste(c[t]," in mg L"^-1)), subtitle="A")

# this part calculates what the thiamin content would be for heating milk with c0=0.4 mg/L for 5 s at 150 C followed by storage at 20 C for 50, 75 and 100 days:

mu_150 <- data.frame(ct=0.4/(1+0.4*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/423))*5/(3600*24)))%>% 
  mutate(ct2=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*10))%>% 
  mutate(ct3=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*50))%>% 
  mutate(ct4=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*75))%>% 
  mutate(ct5=ct/(1+ct*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*100))

plotB <- ggplot(mu_150, aes(x=ct))+
  geom_density(fill="turquoise", alpha=0.5)+
  geom_density(aes(x=ct2), fill="red", alpha=0.5)+
  geom_density(aes(x=ct3), fill="blue",alpha=0.5)+
  geom_density(aes(x=ct4), fill="yellow", alpha=0.5)+
  geom_density(aes(x=ct5), fill="purple", alpha=0.5)+
  geom_vline(xintercept = 0.4, size=3)+
  labs(x=expression(paste(c[t]," in mg L"^-1)), subtitle="B")

# this part calculates what the thiamin content woudl be for heating milk with c0=0.4 mg/L for 20 min at 120 C followed by storage at 20 C for 50, 75 and 100 days::

thiamin_calc <- data.frame(ct120=0.4/(1+0.4*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/393))*20/(60*24)))%>% 
  mutate(ct120_2=ct120/(1+ct120*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*10)) %>% 
  mutate(ct120_3=ct120/(1+ct120*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*50))%>% 
  mutate(ct120_4=ct120/(1+ct120*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*75))%>% 
  mutate(ct120_5=ct120/(1+ct120*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*100)) %>% 
  mutate(ct150=0.4/(1+0.4*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/423))*5/(3600*24)))%>% 
  mutate(ct150_2=ct150/(1+ct150*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*10))%>% 
  mutate(ct150_3=ct150/(1+ct150*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*50))%>% 
  mutate(ct150_4=ct150/(1+ct150*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*75))%>% 
  mutate(ct150_5=ct150/(1+ct150*10^thiamin_all_model_days2_post$b_kb_Intercept*exp(thiamin_all_model_days2_post$b_Ea_Intercept*(1-367.3/293))*100))

plotC <- ggplot(thiamin_calc)+
  geom_density(aes(x=ct120),fill="turquoise", alpha=0.5)+
  geom_density(aes(x=ct120_2), fill="red", alpha=0.5)+
  geom_density(aes(x=ct120_3), fill="blue",alpha=0.5)+
  geom_density(aes(x=ct120_5), fill="purple", alpha=0.5)+
  annotate("text",x=0.287, y=200, label="E", colour="green")+
  annotate("text",x=0.284, y=200, label="F", color="red")+
  annotate("text",x=0.277, y=200, label="G", color="blue")+
  annotate("text",x=0.269, y=200, label="H", color="purple")+
  
  geom_density(aes(x=ct150),fill="turquoise", alpha=0.5)+
  geom_density(aes(x=ct150_2), fill="red", alpha=0.5)+
  geom_density(aes(x=ct150_3), fill="yellow",alpha=0.5)+
   geom_density(aes(x=ct150_5), fill="purple", alpha=0.5)+
  geom_vline(xintercept = 0.4, size=3)+
  labs(x=expression(paste(c[t]," in mg L"^-1)), y = "density")+
  annotate("text",x=0.392, y=3700, label="A", colour="darkblue")+
  annotate("text",x=0.388, y=2000, label="B", color="red")+
  annotate("text",x=0.378, y=600, label="C", color="black")+
  annotate("text",x=0.362, y=400, label="D", color="purple")+
  xlim(c(0.4,0.25))
plotC

# MULTIRESPONSe CASE STUDY #

# Data from PhD thesis Carline Brands:

df_carline <- read.csv(file=here("data_carline.csv"), header=TRUE, sep=",")
fits_carline <- read.csv(file=here("fits_carline.csv"), header=TRUE, sep=",")

# Code for Figure 17:

plot_la_lu <- ggplot(df_carline, aes(x=time))+
  geom_point(aes(y=la), shape=21, size=2, stroke=1, fill="blue")+
  geom_point(aes(y=lu),shape=15, size=2, stroke=1, fill="green")+
  geom_line(data=fits_carline, aes(y=la_mod))+
  geom_line(data=fits_carline, aes(y=lu_mod))+
  annotate("text", x=20, y=140, label="lactose")+
  annotate("text", x=20, y=40, label="lactulose")+
  labs(x="time, min", y="mmol/L")+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
plot_gal_tag <- ggplot(df_carline, aes(x=time))+
  geom_point(aes(y=gal), shape=24, size=2, stroke=1, fill="red")+
  geom_point(aes(y=tag), shape=23, size=2, stroke=1, fill="brown")+
  geom_point(aes(y=formic), shape=21, size=2, stroke=1, fill="turquoise")+
  geom_line(data=fits_carline, aes(y=gal_mod))+
  geom_line(data=fits_carline, aes(y=tag_mod))+
  geom_line(data=fits_carline, aes(y=formic_mod))+
  annotate("text", x=45, y=18.5, label="galactose", color="red")+
  annotate("text", x=50, y=9, label="formic acid", color="purple")+
  annotate("text", x=50, y=4, label="tagatose", color="brown")+
  labs(x="time, min", y="mmol/L")+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
plot_lys_amad <- ggplot(df_carline, aes(x=time))+
  geom_point(aes(y=amadori))+
  geom_point(aes(y=lys), shape=18, size=2, stroke=1, fill="red")+
  geom_line(data=fits_carline, aes(y=lys_mod))+
  geom_line(data=fits_carline, aes(y=amad_mod))+
  annotate("text", x=20, y=14, label="lysine", color="brown")+
  annotate("text", x=20, y=3, label="Amadori")+
  labs(x="time, min", y="mmol/L")+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
plot_lys_amad <- ggplot(df_carline, aes(x=time))+
  geom_point(aes(y=amadori))+
  geom_point(aes(y=lys), shape=18, size=2, stroke=1, fill="red")+
  geom_line(data=fits_carline, aes(y=lys_mod))+
  geom_line(data=fits_carline, aes(y=amad_mod))+
  annotate("text", x=20, y=14, label="lysine", color="brown")+
  annotate("text", x=20, y=3, label="Amadori")+
  labs(x="time, min", y="mmol/L")+
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
plot_la_lu+plot_gal_tag+plot_lys_amad+plot_mel


# CASE STUY 4 ON ALPHA-LACTALBUMIN

# Read data:
aLa <- read.csv(file=here("data","a_LA_control.csv"), header=TRUE, sep=",")
aLa$conc <-aLa$conc/1000
aLa$time <- aLa$time/60
# Code for Figure S14:
aLa_plot <- aLa %>% mutate(temp=str_c("T=", temp, " °C"))
(aLaplot <- ggplot(data=aLa_plot, aes(x=time, y=conc))+
    geom_point(pch=21, size=1, stroke=1.4, fill="#41b6c4")+
    facet_wrap(~temp, scales = "free")+
    labs(x="time in h", y=expression(paste(alpha,"-LA in mg L"^-1)))
)

# Code for individual regressions:

aLa67_3 <- subset(aLa, temp ==67.3, select=time:conc)


nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(1.3,0.1), nlpar = "c0"),
           prior(normal(0.5,0.05), nlpar="k"),
           prior(normal(1.5,0.3), nlpar="nt")
)           

aLa67_3_model <- brm(formula=nlform, data=aLa67_3, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.999, max_treedepth=15), file=here("fits","aLa67_3"))

aLa69_8 <- subset(aLa, temp ==69.8, select=time:conc)


nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(1.3,0.1), nlpar = "c0"),
           prior(normal(0.8,0.05), nlpar="k"),
           prior(normal(1.5,0.3), nlpar="nt")
)           


aLa69_8_model <- brm(formula=nlform, data=aLa69_8, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","aLa69_8"))

aLa72_5 <- subset(aLa, temp ==72.5, select=time:conc)

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(1.4,0.1), nlpar = "c0"),
           prior(normal(1.8,0.1), nlpar="k"),
           prior(normal(1.4,0.1), nlpar="nt")
)           

aLa72_5_model <- brm(formula=nlform, data=aLa72_5, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","aLa72_5"))

aLa75 <- subset(aLa, temp ==75.0, select=time:conc)

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(1.4,0.1), nlpar = "c0"),
           prior(normal(2.5,0.1), nlpar="k"),
           prior(normal(1.3,0.1), nlpar="nt")
)           

aLa75_model <- brm(formula=nlform, data=aLa75, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","aLa75"))

aLa77_5 <- subset(aLa, temp ==77.5, select=time:conc)

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(1.3,0.1), nlpar = "c0"),
           prior(normal(5,0.5), nlpar="k"),
           prior(normal(1,0.1), nlpar="nt")
)           

aLa77_5_model <- brm(formula=nlform, data=aLa77_5, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","aLa77_5"))

aLa79_6 <- subset(aLa, temp ==79.6, select=time:conc)

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(1.4,0.1), nlpar = "c0"),
           prior(normal(8,0.8), nlpar="k"),
           prior(normal(1,0.1), nlpar="nt")
)           

aLa79_6_model <- brm(formula=nlform, data=aLa79_6, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","aLa79_6"))

# Coe for Figure S15
# calculating the regression lines for the individual models: 

# 67.3 C

time.seq <- data.frame(time = seq(from = 0, to = 6, by = 0.1))

fitlin67_3 <-
  fitted(aLa67_3_model, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predlin67_3 <-
  predict(aLa67_3_model, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline67_3 <- ggplot(data = aLa67_3, aes(x = time, y = conc)) +
  geom_point(color = "red", shape = 19, size = 2, alpha = 1) +  
  geom_ribbon(data = predlin67_3, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`), fill = "blue", alpha=0.4) +
  geom_ribbon(data = fitlin67_3, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`), fill = "lightblue") +
  geom_line(data = fitlin67_3, aes(y = Estimate), size = 1/4) +
  labs(x="time in h", y=expression(paste(alpha,"-LA in mg L"^-1)), subtitle = expression(paste("T= 67.3 ", degree*C)))

#69.8 C

time.seq <- data.frame(time = seq(from = 0, to = 6, by = 0.1))

fitlin69_8 <-
  fitted(aLa69_8_model, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predlin69_8 <-
  predict(aLa69_8_model, 
        newdata = time.seq) %>%
  as_tibble() %>%
 bind_cols(time.seq)

regrline69_8 <- ggplot(data = aLa69_8, aes(x = time, y = conc)) +
  geom_point(color = "red", shape = 19, size = 1.5, alpha = 1) +
  geom_ribbon(data = predlin69_8, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.4) +
  geom_ribbon(data = fitlin69_8, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fitlin69_8, aes(y = Estimate), size = 1/4) +
  labs(x="time in h", y=expression(paste(alpha,"-LA in mg L"^-1)), subtitle = expression(paste("T= 69.8 ", degree*C)))

#72.5 C

time.seq <- data.frame(time = seq(from = 0, to = 4, by = 0.1))

fitlin72_5 <-
  fitted(aLa72_5_model, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predlin72_5 <-
  predict(aLa72_5_model, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline72_5 <- ggplot(data = aLa72_5, aes(x = time, y = conc)) +
  geom_point(color = "red", shape = 19, size = 1.5, alpha = 1) +  
  geom_ribbon(data = predlin72_5, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),fill = "blue", alpha=0.4) +
  geom_ribbon(data = fitlin72_5, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fitlin72_5, aes(y = Estimate), size = 1/4) +
  labs(x="time in h", y=expression(paste(alpha,"-LA in mg L"^-1)), subtitle = expression(paste("T= 72.5 ", degree*C)))

#75 C

time.seq <- data.frame(time = seq(from = 0, to = 1, by = 0.04))

fitlin75 <-
  fitted(aLa75_model, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predlin75 <-
  predict(aLa75_model, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline75 <- ggplot(data = aLa75, aes(x = time, y = conc)) +
  geom_point(color = "red", shape = 19, size = 2, alpha = 1) +
  geom_ribbon(data = predlin75, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`), fill = "blue", alpha=0.4) +
  geom_ribbon(data = fitlin75, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),fill = "lightblue") +
  geom_line(data = fitlin75, aes(y = Estimate), size = 1/4) +
  labs(x="time in h", y=expression(paste(alpha,"-LA in mg L"^-1)), subtitle = expression(paste("T= 75.0 ", degree*C)))

# 77.5 C

time.seq <- data.frame(time = seq(from = 0, to = 1, by =0.04))

fitlin77_5 <-
  fitted(aLa77_5_model, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predlin77_5 <-
  predict(aLa77_5_model, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline77_5 <- ggplot(data = aLa77_5,  aes(x = time, y = conc)) +
  geom_point(color = "red", shape = 19, size = 2, alpha = 1) +
  geom_ribbon(data = predlin77_5, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.4) +
  geom_ribbon(data = fitlin77_5, aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),fill = "lightblue") +
  geom_line(data = fitlin77_5, aes(y = Estimate), size = 1/4) +
  labs(x="time in h", y=expression(paste(alpha,"-LA in mg L"^-1)), subtitle = expression(paste("T= 77.5 ", degree*C)))

#79.6 C

time.seq <- data.frame(time = seq(from = 0, to = 0.8, by = 0.02))

fitlin79_6 <-
  fitted(aLa79_6_model, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predlin79_6 <-
  predict(aLa79_6_model, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline79_6 <- ggplot(data = aLa79_6, aes(x = time, y = conc)) +
  geom_point(color = "red", shape = 19, size = 2, alpha = 1) +
  geom_ribbon(data = predlin79_6, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.4) +
  geom_ribbon(data = fitlin79_6, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fitlin79_6, 
            aes(y = Estimate), size = 1/4) +
  
  theme_bw()+
  labs(x="time in h", y=expression(paste(alpha,"-LA, mg/l"^-1)), subtitle = expression(paste("T= 79.6 ", degree*C)))

regrline67_3 + regrline69_8 + regrline72_5 + regrline75 + regrline77_5 + regrline79_6

# Code for Figure S16:

# Figure S16
post_aLa72_5 <- posterior_samples(aLa72_5_model)
cor_725 <- dplyr::select(post_aLa72_5,b_c0_Intercept:sigma)
# change the names of the columns to be displayed in the panels
cor_725<- setNames(cor_725, c(expression(paste("c"[0])), expression(paste("k"[r])), expression(paste("n"[t])),expression(sigma[e])))

# use ggally for a pairs plot
(cor_plot725 <-cor_725  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")),
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")+
            theme(axis.text.x = element_text(angle=70, hjust=1))
    ))

# Code for Figure S17:

aLa_post1 <- posterior_samples(aLa67_3_model)
aLa_post2 <- posterior_samples(aLa69_8_model)
aLa_post3 <- posterior_samples(aLa72_5_model)
aLa_post4 <- posterior_samples(aLa75_model)
aLa_post5 <- posterior_samples(aLa77_5_model)
aLa_post6 <- posterior_samples(aLa79_6_model)

aLa_post_all <-
  bind_rows(
    aLa_post1,
    aLa_post2,
    aLa_post3,
    aLa_post4,
    aLa_post5,
    aLa_post6
  )
iter <- 8000

aLa_post_all <- 
  aLa_post_all %>% 
  mutate(temperature = rep(c("67.3 °C","69.8 °C","72.5 °C","75.0 °C","77.5 °C","79.6 °C"),each = iter))


aLa_post_all <- 
  aLa_post_all %>% 
  mutate(temperature = rep(c("67.3 °C","69.8 °C","72.5 °C","75.0 °C","77.5 °C","79.6 °C"),each = iter))
aLa_post_all %>% 
  ggplot(aes(x = b_c0_Intercept, y = temperature)) +
  stat_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("c"[0]," in mg L"^-1)), y=expression(paste("c"[0]," density for experiment at T =")))


# Code for Figure 18:

aLa_post1 <- posterior_samples(aLa67_3_model)
aLa_post2 <- posterior_samples(aLa69_8_model)
aLa_post3 <- posterior_samples(aLa72_5_model)
aLa_post4 <- posterior_samples(aLa75_model)
aLa_post5 <- posterior_samples(aLa77_5_model)
aLa_post6 <- posterior_samples(aLa79_6_model)

aLa_post_all <-
  bind_rows(
    aLa_post1,
    aLa_post2,
    aLa_post3,
    aLa_post4,
    aLa_post5,
    aLa_post6
  )
iter <- 8000

aLa_post_all <- 
  aLa_post_all %>% 
  mutate(temperature = rep(c("67.3 °C","69.8 °C","72.5 °C","75.0 °C","77.5 °C","79.6 °C"),each = iter))
aLa_post_all %>% 
  ggplot(aes(x = b_nt_Intercept, y = temperature)) +
  stat_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(n[t]), y=expression(paste("n"[t]," density for experiment at T (in °C)")))

# Code for pooled regression, Tref=73.6+273=346.6, pooled model for nt:

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*(temp+273)*(10^kref)*exp(DH*(1-346.6/(temp+273)))*time)^(1/(1-nt)), 
           c0~1,
           nt~1,
           kref~1,
           DH ~1, 
           nl=TRUE)

# set priors for the parameters:
nlprior<-c(prior(normal(1.3,1), nlpar = "c0", lb=0),
           prior(normal(-3.5,1), nlpar="kref"),
           prior(normal(1.3,0.3), nlpar="nt"),
           prior(normal(78, 10), nlpar="DH"),
           prior(cauchy(0,25), class="sigma")
)

aLa_eyring_pooled<-brm(formula=nlform, data=aLa, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, control = list(adapt_delta = 0.95, max_treedepth=12), file=here("fits", "aLa_eyring_pooled"))

aLa_eyring_pooled_post <- posterior_samples(aLa_eyring_pooled) %>% select(-lp__) %>% mutate(DeltaH=b_DH_Intercept*8.314*346.6) %>%mutate(X=10^b_kref_Intercept) %>%  mutate(Y=(((1.381*10^-23)/(6.626*10^-34))*exp(-b_DH_Intercept))) %>% mutate(DeltaS=8.314*log(X/Y)) %>% mutate(DeltaG=DeltaH-347*DeltaS)

# Code for Figure S18:

pairs_pooled <- dplyr::select(aLa_eyring_pooled_post,b_c0_Intercept:X) %>% mutate(DH=DeltaH/1000) %>% dplyr::relocate(X,.after = b_nt_Intercept) %>% dplyr::relocate(DH, .after=b_kref_Intercept) %>% dplyr::select(-b_kref_Intercept) %>% dplyr::select(-b_DH_Intercept) %>% dplyr::select(-DeltaH)
# change the names of the columns to be displayed in the panels
pairs_pooled<- setNames(pairs_pooled, c(expression(paste("c"[0])), expression(paste("n"[t])), expression(paste("k"[ref])), expression(paste(Delta, "H")), expression(sigma[e])))

# use ggally for a pairs plot
(pairs_pooled_plot <-pairs_pooled  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")),
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")+
            theme(axis.text.x = element_text(angle=45, hjust = 1,vjust=0.5))
    ))

# Code for Table S5:

# collect the results from the posterior:
df_c0_pooled <- c(mean(aLa_eyring_pooled_post$b_c0_Intercept), sd(aLa_eyring_pooled_post$b_c0_Intercept),quantile(aLa_eyring_pooled_post$b_c0_Intercept, probs=c(0.025)), quantile(aLa_eyring_pooled_post$b_c0_Intercept, probs=0.975))

df_nt_pooled <- c(mean(aLa_eyring_pooled_post$b_nt_Intercept), sd(aLa_eyring_pooled_post$b_nt_Intercept),quantile(aLa_eyring_pooled_post$b_nt_Intercept, probs=c(0.025)), quantile(aLa_eyring_pooled_post$b_nt_Intercept, probs=0.975))

df_sigma_pooled <- c(mean(aLa_eyring_pooled_post$sigma), sd(aLa_eyring_pooled_post$sigma),quantile(aLa_eyring_pooled_post$sigma, probs=c(0.025)), quantile(aLa_eyring_pooled_post$sigma, probs=0.975))

# recalculate in kJ/mol:
df_DeltaH_pooled <- c(mean(aLa_eyring_pooled_post$DeltaH/1000), sd(aLa_eyring_pooled_post$DeltaH/1000),quantile(aLa_eyring_pooled_post$DeltaH/1000, probs=c(0.025)), quantile(aLa_eyring_pooled_post$DeltaH/1000, probs=0.975))

df_DeltaS_pooled <- c(mean(aLa_eyring_pooled_post$DeltaS), sd(aLa_eyring_pooled_post$DeltaS),quantile(aLa_eyring_pooled_post$DeltaS, probs=c(0.025)), quantile(aLa_eyring_pooled_post$DeltaS, probs=0.975))

#recalculate in kJ/mol:
df_DeltaG_pooled <- c(mean(aLa_eyring_pooled_post$DeltaG/1000), sd(aLa_eyring_pooled_post$DeltaG/1000),quantile(aLa_eyring_pooled_post$DeltaG/1000, probs=c(0.025)), quantile(aLa_eyring_pooled_post$DeltaG/1000, probs=0.975))

#build a dataframe with these results and then interchange rows and columns:
df_aLa_pooled <- as.data.frame(matrix(data=c(df_c0_pooled,df_nt_pooled,df_DeltaH_pooled,df_DeltaS_pooled,df_DeltaG_pooled,df_sigma_pooled), nrow=6, ncol=4, byrow = TRUE))

rownames(df_aLa_pooled) <- c("$c_0 \\text { in mg L}^{-1}$", "$n_t \\text { (-)}$", "$\\Delta H^{o \\ddagger} \\text { in kJ mol}^{-1}$", "$\\Delta S^{o \\ddagger} \\text { in J mol}^{-1} \\text {K}^{-1}$", "$\\Delta G^{o \\ddagger} \\text { in kJ mol}^{-1}$", "$\\sigma_e$")
colnames(df_aLa_pooled) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  df_aLa_pooled,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:aLa-pooled-summary)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# Code for Figure S19

newvary1 <- expand.grid(time=seq(from = 0, to =6, by=0.1),temp=c(67.3))
newvary2 <- expand.grid(time=seq(from = 0, to =5, by=0.1),temp=c(69.8))
newvary3 <- expand.grid(time=seq(from = 0, to =4, by=0.1),temp=c(72.5))
newvary4 <- expand.grid(time=seq(from = 0, to =1, by=0.05),temp=c(75.0))
newvary5 <- expand.grid(time=seq(from = 0, to =1, by=0.05),temp=c(77.5))
newvary6 <- expand.grid(time=seq(from = 0, to =0.8, by=0.04),temp=c(79.6))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6)

aLa_eyring_pooled_NA <- cbind(newvary, predict(aLa_eyring_pooled, newdata=newvary,re_formula = NA)[,-2])
names(aLa_eyring_pooled_NA) <- c("time", "temp", "conc", "lower", "upper")

(multi_fit <-  ggplot(aLa_plot, aes(x=time, y=conc))+
    geom_point(shape=21, fill="blue", color="black")+
    geom_line(data = aLa_eyring_pooled_NA, aes(y = conc), size = 1, colour="blue") + #population level
    facet_wrap(~temp, ncol=3, scales="free_x")+
    labs(x="time in min", y=expression(paste(alpha, "-LA in g L"^-1)))
)

# Code for partially pooled regression, Eyring partially pooled (pp), varying nt, Tref=73.6+273=346.6

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*(temp+273)*(10^kref)*exp(DH*(1-346.6/(temp+273)))*time)^(1/(1-nt)), 
           c0~1,
           nt~1+(1|temp),
           kref~1,
           DH ~1, 
           nl=TRUE)

# set priors for the parameters:
nlprior<-c(prior(normal(1.37,0.1), nlpar = "c0", lb=0),
           prior(normal(-2.3,0.1), nlpar="kref"),
           prior(normal(1.27,0.015), nlpar="nt"),
           prior(normal(73, 10), nlpar="DH"),
           prior(cauchy(0,10), class="sigma")
)

aLa_eyring_pp<-brm(formula=nlform, data=aLa, family = gaussian(), prior = nlprior, warmup=8000, iter=16000, control = list(adapt_delta = 0.999, max_treedepth=12), file=here("fits", "aLa_eyring_pp"))

aLa_eyring_pp_post <- posterior_samples(aLa_eyring_pp) %>% select(-lp__) %>% mutate(DeltaH=b_DH_Intercept*8.314*346.6) %>%mutate(X=10^b_kref_Intercept) %>%  mutate(Y=(((1.381*10^-23)/(6.626*10^-34))*exp(-b_DH_Intercept))) %>% mutate(DeltaS=8.314*log(X/Y)) %>% mutate(DeltaG=DeltaH-347*DeltaS)

# Code for Figure 20:

re_plot_ala <- mcmc_plot(aLa_eyring_pp, pars=c("^r_temp__"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("67.3 °C","69.8 °C","72.5 °C","75 °C","77.5 °C","79.6 °C"))+
  labs(x=expression(paste(Delta,"(n"[t],")")), y=expression(paste("n"[t]," for experiment at T=")))
re_plot_ala

# Code for Figure S20
pairs_pp <- dplyr::select(aLa_eyring_pp_post,b_c0_Intercept:sigma) %>%
  mutate(DeltaH=b_DH_Intercept*8.314*346.6/1000) %>%mutate(X=10^b_kref_Intercept)  %>% dplyr::select(-b_kref_Intercept) %>% select(-b_DH_Intercept) %>% relocate(X,.after = b_nt_Intercept) %>% relocate(DeltaH,.after = X)

# change the names of the columns to be displayed in the panels
pairs_pp<- setNames(pairs_pp, c(expression(paste("c"[0])), expression(paste("n"[t])), expression(paste("k"[ref])), expression(paste(Delta, "H")), expression(sigma[c][0]), expression(sigma[e])))

# use ggally for a pairs plot
(pairs_pp_plot <-pairs_pp  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")),
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")+
            theme(axis.text.x = element_text(angle=90))
    ))

# Code for Table 6:

# collect the results from the posterior:
df_c0_pp <- c(mean(aLa_eyring_pp_post$b_c0_Intercept), sd(aLa_eyring_pp_post$b_c0_Intercept),quantile(aLa_eyring_pp_post$b_c0_Intercept, probs=c(0.025)), quantile(aLa_eyring_pp_post$b_c0_Intercept, probs=0.975))

df_nt_pp <- c(mean(aLa_eyring_pp_post$b_nt_Intercept), sd(aLa_eyring_pp_post$b_nt_Intercept),quantile(aLa_eyring_pp_post$b_nt_Intercept, probs=c(0.025)), quantile(aLa_eyring_pp_post$b_nt_Intercept, probs=0.975))

df_sigma_pp <- c(mean(aLa_eyring_pp_post$sigma), sd(aLa_eyring_pp_post$sigma),quantile(aLa_eyring_pp_post$sigma, probs=c(0.025)), quantile(aLa_eyring_pp_post$sigma, probs=0.975))

# recalculate in kJ/mol:
df_DeltaH_pp <- c(mean(aLa_eyring_pp_post$DeltaH/1000), sd(aLa_eyring_pp_post$DeltaH/1000),quantile(aLa_eyring_pp_post$DeltaH/1000, probs=c(0.025)), quantile(aLa_eyring_pp_post$DeltaH/1000, probs=0.975))

df_DeltaS_pp <- c(mean(aLa_eyring_pp_post$DeltaS), sd(aLa_eyring_pp_post$DeltaS),quantile(aLa_eyring_pp_post$DeltaS, probs=c(0.025)), quantile(aLa_eyring_pp_post$DeltaS, probs=0.975))

#recalculate in kJ/mol:
df_DeltaG_pp <- c(mean(aLa_eyring_pp_post$DeltaG/1000), sd(aLa_eyring_pp_post$DeltaG/1000),quantile(aLa_eyring_pp_post$DeltaG/1000, probs=c(0.025)), quantile(aLa_eyring_pp_post$DeltaG/1000, probs=0.975))

#build a dataframe with these results:
df_aLa_pp <- as.data.frame(matrix(data=c(df_c0_pp,df_nt_pp,df_DeltaH_pp,df_DeltaS_pp,df_DeltaG_pp,df_sigma_pp), nrow=6, ncol=4, byrow = TRUE))

rownames(df_aLa_pp) <- c("$c_0 \\text { (mg/dm}^3)$", "$n_t \\text { (-)}$", "$\\Delta H^{o \\ddagger} \\text { (kJ/mol)}$", "$\\Delta S^{o \\ddagger} \\text { (J/mol/K)}$", "$\\Delta G^{o \\ddagger} \\text { (kJ/mol)}$", "$\\sigma_e$")
colnames(df_aLa_pp) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  df_aLa_pp,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:aLa-pp-summary)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# Code for Figure 21, fits resulting from the multilevel regression:

newvary1 <- expand.grid(time=seq(from = 0, to =6, by=0.1),temp=c(67.3))
newvary2 <- expand.grid(time=seq(from = 0, to =5, by=0.1),temp=c(69.8))
newvary3 <- expand.grid(time=seq(from = 0, to =4, by=0.1),temp=c(72.5))
newvary4 <- expand.grid(time=seq(from = 0, to =1, by=0.02),temp=c(75.0))
newvary5 <- expand.grid(time=seq(from = 0, to =1, by=0.02),temp=c(77.5))
newvary6 <- expand.grid(time=seq(from = 0, to =0.7, by=0.01),temp=c(79.6))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6)

aLa_eyring_NA <- cbind(newvary, predict(aLa_eyring_pp, newdata=newvary,re_formula = NA)[,-2])
names(aLa_eyring_NA) <- c("time", "temp", "conc", "lower", "upper")
aLa_eyring_NA$temp=as.factor(aLa_eyring_NA$temp)

aLa_eyring_NULL <- cbind(newvary, predict(aLa_eyring_pp, newdata=newvary,re_formula = NULL)[,-2])
names(aLa_eyring_NULL) <- c("time", "temp", "conc", "lower", "upper")
aLa_eyring_NULL$temp=as.factor(aLa_eyring_NULL$temp)

aLa <- aLa %>% mutate(temp=str_c("T=", temp, " °C"))
aLa_eyring_NA <- aLa_eyring_NA %>% mutate(temp=str_c("T=", temp, " °C"))
aLa_eyring_NULL <- aLa_eyring_NULL %>% mutate(temp=str_c("T=", temp, " °C"))
(
  multi_fit <-  ggplot(aLa, aes(x=time, y=conc))+
    geom_point(shape=21, fill="blue", color="black")+
    geom_line(data = aLa_eyring_NA, aes(y = conc), size = 1, colour="blue") + #population level
    geom_line(data=aLa_eyring_NULL, aes(x=time, y=conc), color="red")+ # group level
    facet_wrap(~temp, ncol=3, scales="free_x")+
    labs(x="time in h", y=expression(paste(alpha,"-lactalbumin in mg L"^-1)))
)

# Code for Table 7:

aLa_eyring_pooled <- add_criterion(aLa_eyring_pooled, c("loo"), file=here("fits", "aLa_eyring_pooled"))
aLa_eyring_pp <- add_criterion(aLa_eyring_pp, c("loo"), file=here("fits", "aLa_eyring_pp"))

loo_aLa_pooled <- loo(aLa_eyring_pooled)
loo_aLa_pp <- loo(aLa_eyring_pp)

loo_result_aLa <- loo_compare(loo_aLa_pooled,loo_aLa_pp)

elpd_diff_aLa <- c(loo_result_aLa[1,1], loo_result_aLa[2,1])
se_diff_aLa <- c(loo_result_aLa[1,2], loo_result_aLa[2,2])
loo_df_aLa <- data.frame(elpd_diff_aLa, se_diff_aLa)
colnames(loo_df_aLa) <- c("elpd-difference","SE" )
rownames(loo_df_aLa) <- c("partial pooling", "complete pooling")

apa_table(
  loo_df_aLa,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:aLa-loo)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Code for Table S6:

aLa_all_Arrh_post <- aLa_all_Arrh_post %>% dplyr::select(-b_kref_Intercept, -b_Ea_Intercept, -lp__,-kref,-k0) %>% relocate(sigma,.after = lnk0)

p_Arrh1 <- posterior_summary(aLa_all_Arrh_post)


rownames(p_Arrh1) <- c("$n_t (-)$", "$E_a \\text {  in kJ mol}^{-1}$", "$\\text {ln}(k_0) (-)$","$\\sigma_e \\text {  in g L}^{-1}$")
colnames(p_Arrh1) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  p_Arrh1,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:Arrh1-table)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# Code for Table S7:

aLa_all_Arrh2_post <- aLa_all_Arrh2_post %>% dplyr::select(-b_kref_Intercept, -b_Ea_Intercept, -lp__,-kref,-k0) %>% relocate(sigma,.after = lnk0)
p_Arrh2 <- posterior_summary(aLa_all_Arrh2_post)

rownames(p_Arrh2) <- c("$E_a \\text {  in kJ mol}^{-1}$", "$\\text {ln}(k_0) (-)$","$\\sigma_e \\text {  in g L}^{-1}$")
colnames(p_Arrh2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  p_Arrh2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:Arrh2-table)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)


# CASE STUDY 4 ON MICROBIAL INACTIVATION

# Code to retrieve the data and for Figure 23::

cava <- read.csv(file=here("data", "cava.csv"), header=TRUE, sep=";")
cava_55 <- subset(cava, temp ==55, select=time:logN)
cava_58 <- subset(cava, temp ==58, select=time:logN)
cava_60 <- subset(cava, temp ==60, select=time:logN)
cava_62 <- subset(cava, temp ==62, select=time:logN)

cava_plot <- cava %>% mutate(temp=str_c("T=",temp," °C"))

cava_plot %>% ggplot(aes(x=time, y=logN))+
  geom_point(shape=21, size=1.2, stroke=2, color="black", fill="red")+
  facet_wrap(~temp, scales="free")+
  labs(x="time in min", y=expression(paste(log[10],N[0])))+
  theme(strip.background = element_rect(color="black", fill="lightblue", size=1.5, linetype="solid"))

# Code for Figure S21:

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha", lb=0),
           prior(normal(0.01,0.01), nlpar="beta",lb=0),
           prior(cauchy(0,10), class="sigma")
)

cava_ppc <- brm(formula=nlform, data=cava_55, family = gaussian(), sample_prior="only", prior = nlprior, chains= 4, warmup=1000, iter = 2000, control = list(adapt_delta = 0.9), file=here("fits", "cava_ppc"))

fe_only <- tibble(time=seq(min(cava_55$time),max(cava_55$time), length.out = 100)) %>% add_fitted_draws(cava_ppc, re_formula = NA, n=200)

ggplot(fe_only, aes(x=time, y=.value, group=.draw))+
  geom_line(colour="blue")+
  ylim(0,8)


# Code for individual regressions with the full Weibull model

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,0.1), nlpar = "logN0"),
           prior(normal(1.5,0.1), nlpar="alpha"),
           prior(normal(15,2), nlpar="beta"),
           prior(cauchy(0,10), class="sigma")
)

cava_55_model <- brm(formula=nlform, data=cava_55, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.95), file=here("fits", "cava_55_model"))

post_cava_55 <- posterior_samples(cava_55_model)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1.5,0.1), nlpar="alpha"),
           prior(normal(7,2), nlpar="beta"),
           prior(cauchy(0,10), class="sigma")
)

cava_58_model <- brm(formula=nlform, data=cava_58, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.9), file=here("fits", "cava_58_model"))

post_cava_58 <- posterior_samples(cava_58_model)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,0.1), nlpar = "logN0"),
           prior(normal(1.2,0.1), nlpar="alpha"),
           prior(normal(1.7,0.1), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

cava_60_model <- brm(formula=nlform, data=cava_60, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.95), file=here("fits", "cava_60_model"))

post_cava_60 <- posterior_samples(cava_60_model)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,0.1), nlpar = "logN0"),
           prior(normal(1.3,0.1), nlpar="alpha"),
           prior(normal(1.0,0.1), nlpar="beta"),
           prior(cauchy(0,10), class="sigma")
)

cava_62_model <- brm(formula=nlform, data=cava_62, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.9), file=here("fits", "cava_62_model"))

post_cava_62 <- posterior_samples(cava_62_model)

# Code for Figure 23:

p_cava <-
  bind_rows(
    post_cava_55,
    post_cava_58,
    post_cava_60,
    post_cava_62
  )
iter <- 8000

p_cava <- 
  p_cava %>% 
  mutate(Temperature = rep(c("55 °C","58 °C","60 °C","62 °C"),each = iter)) 

plot_logN0 <- p_cava %>% 
  ggplot(aes(x = b_logN0_Intercept, y = Temperature)) +
  geom_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste(log[10],N[0])), y="densities for the experiment at T =")


plot_alpha <- p_cava %>% 
  ggplot(aes(x = b_alpha_Intercept, y = Temperature)) +
  geom_halfeyeh(fill = "lightblue", point_interval = mean_qi, .width = .95) +
  labs(x=expression(alpha))

plot_beta <- p_cava %>% 
  ggplot(aes(x = log10(b_beta_Intercept), y = Temperature)) +
  geom_halfeyeh(fill = "lightblue", point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste(log[10],beta)))

plot_all=plot_logN0+plot_alpha+plot_beta

plot_all[[2]]=plot_all[[2]]+theme(axis.text.y=element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank() )

plot_all[[3]]=plot_all[[3]]+theme(axis.text.y=element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank() )

plot_all

# Code for pooled regression, Tref chosen as the average temperature: 58.75

nlform<-bf(logN ~ logN0-(1/2.303)*(time/(10^bref*10^(-c1*(temp-58.75))))^alpha,
           logN0~1, 
           c1~1, 
           bref~1,
           alpha~1, 
           nl=TRUE)

nlprior<-c(prior(normal(7,0.01), nlpar = "logN0"),
           prior(normal(1.3,0.1), nlpar="alpha"),
           prior(normal(0.1,0.1), nlpar="c1"),
           prior(normal(1,1), nlpar="bref")
)

cava_all_pooled <- brm(formula=nlform, data=cava, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, control = list(adapt_delta = 0.9), file=here("fits", "cava_all_pooled"))
cava_all_pooled_post <- posterior_samples(cava_all_pooled)

# Code for Table S8

p_cava1 <- posterior_summary(cava_all_pooled_post)

p_cava2 <- p_cava1[c(1,2,5,4,3),]

rownames(p_cava2) <- c("$log_{10} N_0 (-)$", "$\\alpha (-)$", "$\\beta_{\\text{ref}} \\text{  in min} $","$z \\text{  in °C}$", "$\\sigma_e (-)$")
colnames(p_cava2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  p_cava2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:summary-bigelow-single)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# code for Figure S22:

newvary1 <- expand.grid(time=seq(from = 0, to =100, by=2),temp=c(55))
newvary2 <- expand.grid(time=seq(from = 0, to =25, by=0.5),temp=c(58))
newvary3 <- expand.grid(time=seq(from = 0, to =13, by=0.3),temp=c(60))
newvary4 <- expand.grid(time=seq(from = 0, to =7, by=0.2),temp=c(62))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4)

cava_pooled <- cbind(newvary, predict(cava_all_pooled, newdata=newvary)[,-2])
cava_plot <- cava %>% mutate(temp=str_c("T=",temp," °C"))
names(cava_pooled) <- c("time", "temp", "logN", "lower", "upper")
cava_pooled <- cava_pooled %>% mutate(temp=str_c("T=", temp, " °C"))

(pooled_fit <-  ggplot(cava_plot, aes(x=time, y=logN))+
    geom_point(shape=21, fill="blue", color="black")+
    geom_line(data = cava_pooled, aes(y = logN), size = 1, colour="blue") +
    geom_line(data = cava_pooled, aes(y = lower), lty = 2) +
    geom_line(data = cava_pooled, aes(y = upper), lty = 2) +
    facet_wrap(~temp, ncol=5, scales="free_x")+
    labs(x="time in min", y=expression(paste("log"[10],"N"))))


# Code for partially pooled regression:

nlform<-bf(logN ~ logN0-(1/2.303)*(time/(10^bref*10^(-c1*(temp-58.75))))^alpha,
           logN0~1, 
           c1~1, 
           bref~1,
           alpha~1+(1|temp), 
           nl=TRUE)

nlprior<-c(prior(normal(7,0.1), nlpar = "logN0"),
           prior(normal(1.3,0.1), nlpar="alpha"),
           prior(normal(0.18,0.01), nlpar="c1"),
           prior(normal(0.6,0.06), nlpar="bref")
)

cava_all_pooled_multi <- brm(formula=nlform, data=cava, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits", "cava_all_pooled_multi"))

#build dataframe from posterior then remove some columns and recalculate Z and beta_ref
cava_all_pooled_multi_post <- posterior_samples(cava_all_pooled_multi)%>% select(-`r_temp__alpha[55,Intercept]`,-`r_temp__alpha[58,Intercept]`,-`r_temp__alpha[60,Intercept]`,-`r_temp__alpha[62,Intercept]`,-lp__)%>% mutate(Z=1/b_c1_Intercept) %>% mutate(beta_ref=10^b_bref_Intercept) %>% select(-b_c1_Intercept) %>% select(-b_bref_Intercept)

# Code for Figure 24:

re_plot_cava <- mcmc_plot(cava_all_pooled_multi, pars=c("^r_temp__"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("55 °C","58 °C","60 °C","62 °C"))+
  labs(x=expression(paste(Delta,"(",alpha,")")), y=expression(paste(alpha," for experiment done at T =")))

re_plot_cava

# Code for Figure 25:

# combines the data with predictions using the random effects with re_formula = NULL and NA
newvary1 <- expand.grid(time=seq(from = 0, to =100, by=2),temp=c(55))
newvary2 <- expand.grid(time=seq(from = 0, to =25, by=0.5),temp=c(58))
newvary3 <- expand.grid(time=seq(from = 0, to =13, by=0.3),temp=c(60))
newvary4 <- expand.grid(time=seq(from = 0, to =7, by=0.2),temp=c(62))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4)

cava_multi_NA <- cbind(newvary, predict(cava_all_pooled_multi, newdata=newvary, re_formula = NA)[,-2])
cava_multi_NULL <- cbind(newvary, predict(cava_all_pooled_multi, newdata=newvary, re_formula = NULL)[,-2])
cava_multi_NULL$temp=as.character(cava_multi_NULL$temp)
cava_multi_NULL <- cava_multi_NULL %>% mutate(temp=str_c("T=", temp, " °C"))
names(cava_multi_NA) <- c("time", "temp", "logN", "lower", "upper")
names(cava_multi_NULL) <- c("time", "temp", "logN", "lower", "upper")

cava_multi_NA <- cava_multi_NA %>% mutate(temp=str_c("T=", temp, " °C"))


(cava_fit_multi_NULL <-  ggplot(cava_plot, aes(x=time, y=logN))+
    geom_point(shape=21, size=1.2, stroke=2,color="black", fill="red")+
    facet_wrap(~temp,ncol=5, scales="free_x")+
    geom_line(data = cava_multi_NULL, aes(y = logN), size = 1, colour="blue") +
    geom_line(data=cava_multi_NA, aes(y=Estimate), color="red", lty=2, size=1.5)+
    labs(x="time in min", y=expression(paste(log[10],N[0])))+
    theme(strip.background = element_rect(color="black", fill="lightblue", size=1.2, linetype="solid"))
)

# Code for Table 8:

p_cava_multi1 <- posterior_summary(cava_all_pooled_multi_post)
p_cava_multi2 <- p_cava_multi1[c(1,2,6,5,3,4),]

rownames(p_cava_multi2) <- c("$\\log_{10} N_0$", "$\\alpha (-)$", "$\\beta_{ref} ( min) $","$Z (^oC)$","$\\sigma_{\\alpha}$", "$\\sigma_e$")
colnames(p_cava_multi2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  p_cava_multi2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:cava-multi-summary)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# Code for Table 9:

cava_all_pooled_multi <- add_criterion(cava_all_pooled_multi, c("loo"), file=here("fits", "cava_all_pooled_multi"))
cava_all_pooled <- add_criterion(cava_all_pooled, c("loo"), file=here("fits", "cava_all_pooled"))

loo_cava_pooled <- loo(cava_all_pooled)
loo_cava_pp <- loo(cava_all_pooled_multi)

loo_result_cava <- loo_compare(loo_cava_pooled,loo_cava_pp)

elpd_diff_cava <- c(loo_result_cava[1,1], loo_result_cava[2,1])
se_diff_cava <- c(loo_result_cava[1,2], loo_result_cava[2,2])
loo_df_cava <- data.frame(elpd_diff_cava, se_diff_cava)
colnames(loo_df_cava) <- c("elpd-difference","SE" )
rownames(loo_df_cava) <- c("partial pooling", "complete pooling")

apa_table(
  loo_df_cava,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:aLa-loo)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Code for Figure S23: applying log-linear regression to the cava data to obtain D-values at each temperature:
cava55D <- 
  brm(data = cava_55, family = gaussian,
      formula = logN ~ 1 + time,
      prior = c(set_prior("normal(7, 1)", class = "Intercept"),
                set_prior("normal(2, 1)", class = "b"),
                set_prior("cauchy(0,5)", class = "sigma")),
      chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta =0.95), file=here("fits", "cava55D"))

cava55D_plot <- ggplot(cava_55, aes(x=time, y=logN))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="time (min)", y=expression(paste("log"[10],"N")), subtitle = expression(paste("T= 55 ", degree*C)))

cava58D <- 
  brm(data = cava_58, family = gaussian,
      formula = logN ~ 1 + time,
      prior = c(set_prior("normal(7, 1)", class = "Intercept"),
                set_prior("normal(2, 1)", class = "b"),
                set_prior("cauchy(0,5)", class = "sigma")),
      chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta =0.95), file=here("fits", "cava58D"))

cava58D_plot <- ggplot(cava_58, aes(x=time, y=logN))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="time (min)", y=expression(paste("log"[10],"N")), subtitle = expression(paste("T= 58 ", degree*C)))


cava60D <- 
  brm(data = cava_60, family = gaussian,
      formula = logN ~ 1 + time,
      prior = c(set_prior("normal(7, 1)", class = "Intercept"),
                set_prior("normal(2, 1)", class = "b"),
                set_prior("cauchy(0,5)", class = "sigma")),
      chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta =0.95), file=here("fits", "cava60D"))

cava60D_plot <- ggplot(cava_60, aes(x=time, y=logN))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="time (min)", y=expression(paste("log"[10],"N")), subtitle = expression(paste("T= 60 ", degree*C)))


cava62D <- 
  brm(data = cava_62, family = gaussian,
      formula = logN ~ 1 + time,
      prior = c(set_prior("normal(7, 1)", class = "Intercept"),
                set_prior("normal(2, 1)", class = "b"),
                set_prior("cauchy(0,5)", class = "sigma")),
      chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta =0.95), file=here("fits", "cava62D"))

cava62D_plot <- ggplot(cava_62, aes(x=time, y=logN))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="time (min)", y=expression(paste("log"[10],"N")), subtitle = expression(paste("T= 62 ", degree*C)))

cava55D_post <- posterior_samples(cava55D) %>% select(-lp__) %>% mutate(D55=-1/b_time) 
cava58D_post <- posterior_samples(cava58D) %>% select(-lp__) %>% mutate(D58=-1/b_time) 
cava60D_post <- posterior_samples(cava60D) %>% select(-lp__) %>% mutate(D60=-1/b_time) 
cava62D_post <- posterior_samples(cava62D) %>% select(-lp__) %>% mutate(D62=-1/b_time) 

cava55D_plot+cava58D_plot+cava60D_plot+cava62D_plot

# Code for Figure S24:

#Tref=58.75
T_df <- c(55,58,60,62)
D_df <- c(mean(cava55D_post$D55), mean(cava58D_post$D58), mean(cava60D_post$D60), mean(cava62D_post$D62))
TDT_df <- data.frame(T_df,D_df) %>% mutate(logD=log10(D_df)) %>% mutate(T_c=T_df-58.75)

TDT_model <- 
  brm(data = TDT_df, family = gaussian,
      formula = logD ~ 1 + T_c,
      prior = c(set_prior("normal(1, 1)", class = "Intercept"),
                set_prior("normal(-0.16, 0.01)", class = "b"),
                set_prior("cauchy(0,1)", class = "sigma")),
      chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta =0.95), file=here("fits", "TDT_model"))
TDT_model_post <- posterior_samples(TDT_model) %>% select(-lp__) 
#print(TDT_model, digits=4)
#plot(conditional_effects(TDT_model), points=TRUE)

TDT_plot <- plot(conditional_effects(TDT_model), plot=F)[[1]]
TDT_plot+geom_point(data=TDT_df, aes(x=T_c, y=logD), inherit.aes = F)+
  labs(x=expression(paste(italic("T")[c])), y=expression(paste("log"[10], italic("D"))))

# Code for Figure 26, calculation of t for 6D reduction at 60 C:

cava_all_pooled_multi_post <- cava_all_pooled_multi_post %>% mutate(t6D60=6*(beta_ref)*10^((58.75-60)/Z)*(log(10))^(1/b_alpha_Intercept))
TDT_model_post <- TDT_model_post %>% mutate(t6D60=6*10^b_Intercept*10^(b_T_c*(60-58.75)))

(t60_density <- ggplot(data=cava_all_pooled_multi_post)+
    geom_density(fill="red", aes(x=t6D60))+
    geom_vline(xintercept=mean(cava_all_pooled_multi_post$t6D60), lty=2)+
    geom_density(data=TDT_model_post, fill="turquoise", alpha=0.3, aes(x=t6D60))+
    geom_vline(xintercept=mean(TDT_model_post$t6D60), lty=2)+
    labs(x="time for 6 decimal reductions at 60 °C in min")+
    xlim(0,50)
)

