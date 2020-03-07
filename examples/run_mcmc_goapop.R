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
ess <- monitor(fit.nuts$samples, print=FALSE)[,'n_eff']
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

## Make marginal plots of posterior vs MLE estimates
posterior <- extract_samples(fit.nuts)
ff <- function(par.num){
  par <- names(posterior)[par.num]
  mle <- fit.nuts$mle$est[par.num]
  se <-  fit.nuts$mle$se[par.num]
  x1 <- seq(qnorm(.001, mle, se), qnorm(.999, mle, se), len=100)
  y1 <- dnorm(x1, mle, se)
  data.frame(par=par, mle=mle, x=x1, y=y1)
}
fff <- function(par.num){
  par <- names(posterior)[par.num]
  post <- posterior[,par.num]
 ## tmp <- hist(post, plot=FALSE, breaks=30)
  tmp <- density(post)
  data.frame(par=par, x=tmp$x, y=tmp$y, ESS=as.numeric(ess[par.num]))
}
posts <- mles <- list()
for(i in 1:ncol(posterior)){
  mles[[i]] <- ff(i)
  posts[[i]] <- fff(i)
}
mles <- do.call(rbind,mles)
posts <- do.call(rbind, posts)

## Since can be too many parameters, break them up into
## pages. Stolen from
## http://stackoverflow.com/questions/22996911/segment-facet-wrap-into-multi-page-pdf
par.names <- names(posterior)
noVars <- length(par.names); noPlots <- 20
plotSequence <- c(seq(0, noVars-1, by = noPlots), noVars)
pdf('plots/goapop_marginals.pdf', onefile=TRUE, width=7,height=7)
for(ii in 2:length(plotSequence)){
  start <- plotSequence[ii-1] + 1; end <- plotSequence[ii]
  z1 <- subset(mles, par %in% par.names[start:end])
  z2 <- subset(posts, par %in% par.names[start:end])
  ## tmp <- z2 %>% group_by(par) %>% summarize( ESS=paste0("ESS=",ESS[1])) %>%
  ##   ungroup() %>% as.data.frame()
  g <- ggplot() +
    geom_line(data=z1, aes(x,y), lwd=1) +
    geom_line(data=z2, aes(x,y), lwd=1, color='red') +
    facet_wrap('par', scales='free')
  g <- g+ theme(text=element_text(size=12),
                axis.text.y=element_blank(),
                panel.spacing = unit(0, "lines"))+ labs(x=NULL, y=NULL)
  print(g)
}
dev.off()
