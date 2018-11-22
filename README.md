# Bayesian-mixed-latent-Markov-models-for-binary-diagnostic-data

The 3 files are the JAGS model files for the basic time homogeneous, basic time heterogeneous and the mixed time homogeneous latent Markov models.

You will need to prep your data via an R script where you then call the jags model files via something like (assuming you use rjags; slightly different for R2jags etc):

`jagsTimeHomoModel <- jags.model('scripts/jagsModelFile_TimeHomo_general.jags', data=dat, n.chains = 4, n.adapt = 1000)`

`parsTimeHomoModel <- coda.samples(model=jagsTimeHomoModel,variable.names=c('pInfInit','pUninf2Inf','pInf2Inf','crp'),n.iter=10000)`

(and similarly for the bsic time heterogeneous and the mixed time homogeneous models)
