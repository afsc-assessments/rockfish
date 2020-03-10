library(ggplot2)
library(tidyr)
library(dplyr)
theme_set(theme_bw())
library(adnuts)
packageVersion('adnuts')                # needs to be >1.0.1
library(snowfall)
library(rstan)
library(shinystan)
chains <- parallel::detectCores()-1 # chains to run in parallel
## Reproducible seeds are passed to ADMB
set.seed(352)
seeds <- sample(1:1e4, size=chains)

## model name and directory (path)
m <- 'pop'; d <- 'GOAPOP'

## Try RWM to see potential mismatches with covar matrix
setwd(d)
system(paste(m, "-version"))
system(paste(m, "-binp pop.bar -nox"))
setwd('..')

tt <- 2000 # thin rate
fit.rwm <- sample_admb(m, d, iter=750*tt, warmup=250*tt, thin=tt, algorithm='RWM',
                   chains=chains, cores=chains, parallel=TRUE)
ess <- monitor(fit.rwm$samples, print=FALSE)[,'n_eff']
barplot(log10(sort(ess)))
png('plots/goapop_pairs_slow_rwm.png', width=6, height=4, units='in', res=600)
slow <- names(sort(ess, FALSE))[1:8]
pairs_admb(fit.rwm, pars=slow)
dev.off()
fast <- names(sort(ess, TRUE))[1:8]
png('plots/goapop_pairs_fast_rwm.png', width=6, height=4, units='in', res=600)
pairs_admb(fit.rwm, pars=fast)
dev.off()
saveRDS(object=fit.rwm, file='fits/goapop_fit_rwm.rds')
## launch_shinyadmb(fit.rwm)


### Try NUTS sampling
setwd(d)
system(paste(m, '-hbf -ainp pop.par -nox'))
setwd('..')
inits <- NULL
fit.nuts0 <- sample_admb(m, d, iter=500, warmup=100, init=inits, algorithm='NUTS',
                   chains=chains, cores=chains, parallel=TRUE, seeds=1:chains,
                   control=list(adapt_delta=.9, metric='mle'))

## Restart using that estimated covar
inits <- sample_inits(fit.nuts0, chains)
fit.nuts <- sample_admb(m, d, iter=3000, warmup=600, init=inits, algorithm='NUTS',
                   chains=chains, cores=chains, parallel=TRUE, seeds=1:chains,
                   thin=3,
                   control=list(metric=fit.nuts0$covar.est, adapt_delta=.8))
saveRDS(object=fit.nuts, file='fits/goapop_fit_nuts.rds')
mon <- monitor(fit.nuts$samples, warmup=fit.nuts$warmup, print=FALSE)
ess <- mon[,'n_eff']

barplot(log10(sort(ess)))
fast <- names(sort(ess, TRUE))[1:8]
png('plots/goapop_pairs_fast_nuts.png', width=6, height=4, units='in', res=600)
pairs_admb(fit.nuts, pars=fast)
dev.off()
slow <- names(sort(ess, FALSE))[1:8]
png('plots/goapop_pairs_slow_nuts.png', width=6, height=4, units='in', res=600)
pairs_admb(fit.nuts, pars=slow)
dev.off()
## Weirdness with the recdevs
pars <- grep('log_rec_', x=names(ess))
png('plots/goapop_pairs_recdevs1.png', width=6, height=4, units='in', res=600)
pairs_admb(fit.nuts, pars=c(2,pars[22:28]), label.cex=.5)
dev.off()
png('plots/goapop_pairs_recdevs2.png', width=6, height=4, units='in', res=600)
pairs_admb(fit.nuts, pars=c(2,pars[70:78]), label.cex=.5)
dev.off()
## Some key fixed effects
pars <- c('log_mean_rec', "sigr", "a503", "delta3", "delta_srv1",
          "logm", "mF50", "mF40", "mF35")
png('plots/goapop_pairs_FE.png', width=6, height=4, units='in', res=600)
pairs_admb(fit.nuts, pars=pars, label.cex=.5)
dev.off()
## launch_shinyadmb(fit.nuts)


pdf('plots/goapop_marginals.pdf', onefile=TRUE, width=7,height=5)
plot_marginals(fit.nuts, mon)
dev.off()
