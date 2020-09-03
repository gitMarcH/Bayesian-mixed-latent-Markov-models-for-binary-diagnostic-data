# Bayesian mixed latent Markov models for binary diagnostic data

## JAGS model files

* jagsModelFile_TimeHomo_general.jags
* jagsModelFile_TimeHetero_general.jags
* jagsModelFile_TimeHomoMixed_general.jags

These 3 files are the JAGS model files for the basic time homogeneous, basic time heterogeneous and the mixed time homogeneous latent Markov models.

You will need to prep your data via an R script where you then call the jags model files via something like (assuming you use rjags; slightly different for R2jags etc):

`jagsTimeHomoModel <- jags.model('scripts/jagsModelFile_TimeHomo_general.jags', data=dat, n.chains = 4, n.adapt = 1000)`

`parsTimeHomoModel <- coda.samples(model=jagsTimeHomoModel,variable.names=c('pInfInit','pUninf2Inf','pInf2Inf','crp'),n.iter=10000)`

(and similarly for the basic time heterogeneous and the mixed time homogeneous models)

## R data processing, model fitting, visualisation and correlation computation scripts

* jagsModelling_general.R (model fitting)
* summaryRealDatasetLMM.R (visualising the model estimated sensitivities, specificities, PPVs and NPVs)
* lmm_Correlations.R (calculating correlations between Ct measurements from the various tests)

These are R scripts that were used to produce the results from an upcoming publication.
