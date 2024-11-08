



load(file = "data/temp_error.rdata")
nsamples <- 50
burnin <- 10
colnames(X) <- paste0("pat", 1:ncol(X))
rownames(X) <- rownames(PriorFixed1)
PriorFixed1
out <- CompressiveNMF(X, K=0, nsamples = nsamples, burnin = burnin, S = pmax(PriorFixed1, 1e-10),
               betah = rep(1, ncol(PriorFixed1)), swap_prior = FALSE, alpha = 0.5)






