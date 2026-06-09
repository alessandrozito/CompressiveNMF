#' @export
logFullCond_a <-  function(x, Mu, Theta, epsilon, c0, d0){
  K <- length(Mu)
  J <- ncol(Theta)
  K * J * x * log(x) + (K * J * x + K)*log(epsilon*x*J) - K * J * lgamma(x) - K * lgamma(x * J + 1) -
    x * (2 * J * sum(log(Mu)) - sum(log(Theta)) + sum(rowSums(Theta)/Mu) +
           epsilon * J * sum(1/Mu)) + (c0 - 1) * log(x) - x * d0
}

#' @export
sample_FullCond_a <- function(nsamples, Mu, Theta, epsilon, c0, d0,  a_start = 1) {
  UppU <- -nlminb(start = a_start, objective = function(x) -logFullCond_a(x, Mu, Theta, epsilon, c0, d0),
                  lower = 1e-6)$objective

  UppV <- -nlminb(start = a_start, objective = function(x) -logFullCond_a(x, Mu, Theta, epsilon, c0, d0) - 2 * log(x),
                  lower = 1e-6)$objective

  s <- exp(0.5 * (UppV - UppU))
  samples <- rep(NA, nsamples)
  for(n in 1:nsamples){
    accepted <- 0
    while(accepted == 0) {
      u <- runif(1); v <- runif(1)
      x <- s * v / u
      if(log(u) < 0.5 * (logFullCond_a(x, Mu, Theta, epsilon, c0, d0) - UppU)){
        samples[n] = x
        accepted = 1
      }
    }
  }
  samples
}

# This function runs the Gibbs sampler for the compressive non-negative matrix
# factorization algorithm. We consider the case when also known signatures are
# included in the analysis, with the possibility to update them.

#' Function for our compressive NMF methods
#' @param X The input data matrix.
#' @param K Number of signatures.
#' @param nsamples Number of samples in each MCMC chain.
#' @param burnin Number of burn-in iterations in each MCMC chain.
#' @param alpha Parameter for the Dirichlet prior in the signature
#' @param a Parameter for the gamma prior over the loadings
#' @param epsilon Parameter controlling the mean of the Compressive hyperprior
#' @param nchains Number of chains to run in parallel.
#' @param ncores Number of CPU cores to use.
#' @param a0 Shape parameter for the inverse gamma compressive hyperprior. Default is \code{a * ncol(X) + 1}
#' @param b0 Rate parameter for the inverse gamma compressive hyperprior. Default is \code{a * epsilon * ncol(X)}
#' @param random_a Whether a should be treated as a random quantity
#' @param c0 Prior for a
#' @param d0 Prior for a
#' @param S Matrix of informative prior for the signatures
#' @param cutoff_excluded Cutoff for excluded signatures
#' @param use_cosmic Logical, indicating whether to use Cosmic data.
#' @param swap_prior Logical, indicating whether to swap the informative prior at 2/3 of the burnin phase.
#' @param betah Parameter controlling the strenght of the prior for the informative priors in S. Length has to be equal to the number of columns of \code{S}
#' @param betah_optimal Logical, indicating whether to use an optimal betah. See \code{src/tune_betah.cpp}
#' @param progressbar Logical, indicating whether to show a progress bar.
#' @param verbose Logical, indicating whether to display verbose output.
#' @param seed Seed for reproducibility.
#' @importFrom foreach %dopar% foreach
#' @return Returns a list of results from CompressiveNMF.
#' @useDynLib CompressiveNMF
#' @export
CompressiveNMF_randa <- function(X,
                           K = 20,
                           nsamples = 1000,
                           burnin = 2000,
                           alpha = 0.5,
                           a = 1,
                           epsilon = 0.001,
                           nchains = 1,
                           ncores = 1,
                           a0 = a * ncol(X) + 1,
                           b0 = a * epsilon * ncol(X),
                           random_a = FALSE,
                           c0 = 1,
                           d0 = 1,
                           S = NULL,
                           cutoff_excluded = 5 * epsilon,
                           use_cosmic = FALSE,
                           swap_prior = TRUE,
                           betah = 100,
                           betah_optimal = TRUE,
                           verbose = TRUE,
                           seed = 10) {

  # Check if the matrix X is made of integers
  if(any(round(X) != X)){
    stop("The matrix X has to be made of integers! Consider using round(X).")
  }
  # Matrix dimensions
  I <- nrow(X); J <- ncol(X)

  # Record time of execution
  t_start <- Sys.time()

  # Check if the matrix S is correctly specified
  if(!is.null(S)){
    if(!is.matrix(S)){
      stop("The hyperparameter S must be a matrix!")
    }
  }

  # If we use the Cosmic database, then overwrite the matrix S
  if(use_cosmic){
    if (I == 96) {
      cat("Using informative prior for SBS mutations from COSMICv3.4, genome build GRCh37. Artifacts excluded \n")
      cosmic_data <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4
      S <- as.matrix(cosmic_data[, -c(1,2,3)])
      rownames(S) <- cosmic_data$Channel
      # Re-align the names of S according to the names of X
      S <- S[rownames(X), ]
    } else if (I == 83) {
      cat("Using informative prior for ID mutations from COSMICv3.4, genome build GRCh37. \n")
      S <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37
      S <- S[rownames(X), ]
    } else if (I == 78) {
      cat("Using informative prior for DBS mutations from COSMICv3.4, genome build GRCh37. \n")
      S <- CompressiveNMF::COSMIC_v3.4_DBS78_GRCh37
      S <- S[rownames(X), ]
    } else {
      stop("Number of rows of X is not equal to the length of any of the COSMIC signatures.")
    }
  }

  # Extract information from the data
  if(!is.null(S)){
    H <- ncol(S)
    Ktot <- K + H
    prior_names <- colnames(S)
    if(K > 0){
      prior_names <- c(prior_names, paste0("Sig_new", c(1:K)))
    }
  } else {
    Ktot <- K
    prior_names <- paste0("Sig_new", c(1:K))
  }


  # Prepare the hyperparameters in the signatures
  SignaturePrior <- matrix(alpha, nrow = I, ncol = Ktot)

  if(!is.null(S)){
    if(betah_optimal & use_cosmic){
      if(I == 96){
        betah <- CompressiveNMF::Betah_SBS96_GRCh37_v3.4
      } else if (I == 83){
        betah <- CompressiveNMF::Betah_ID83_GRCh37_v3.4
      } else if (I == 78) {
        betah <- CompressiveNMF::Betah_DBS78_GRCh37_v3.4
      }
      betah <- betah[colnames(S)]
    }
    # Check that betah and S have the same dimension
    if(length(betah) == 1){
      betah <- rep(betah, ncol(S))
    } else if (length(betah) != ncol(S)) {
      stop("length of betah must be the same as number of columns of S")
    }
    SignaturePrior[, 1:H] <- t(betah)[rep(1, I), ] * S
  }
  colnames(SignaturePrior) <- prior_names
  SignaturePrior <- pmax(SignaturePrior, 1e-10) # correction for Dirichlet
  # Find the number of non-zero entries to speed-up the sampler
  nonzero_ids <- which(X != 0, arr.ind = TRUE)

  # Pack the model parameters
  model_pars <- list()
  model_pars$a <- a
  model_pars$epsilon <- epsilon
  model_pars$random_a <- random_a
  model_pars$alpha <- alpha
  model_pars$a0 <- a0
  model_pars$b0 <- b0
  model_pars$c0 <- c0
  model_pars$d0 <- d0
  model_pars$cutoff_excluded <- cutoff_excluded
  model_pars$SignaturePrior <- SignaturePrior
  model_pars$nonzero_ids <- nonzero_ids
  model_pars$verbose <- verbose
  model_pars$swap_prior <- swap_prior
  model_pars$betah <- betah
  model_pars$S <- S

  # Run the sampler across different chains
  if(ncores > 1){
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    mcmc_out <- foreach(i = 1:nchains, .packages = c("CompressiveNMF")) %dopar% {
      library(CompressiveNMF)
      if(verbose){
        verbose.out <- paste0("chain_status_", i,".txt")
      } else {
        verbose.out <- NULL
      }
      # Step 1 - initialize the sampler
      init_pars <- Initialize_CompressiveNMF(X, model_pars)
      # Step 2 - run the sampler
      if(swap_prior & !is.null(S)){
        # Run the sampler for 3/4, and the perform the swap
        sample_preswap <- Run_mcmc_CompressiveNMF(X,
                                                  model_pars,
                                                  init_pars,
                                                  nsamples = 1,
                                                  burnin = round(burnin * 0.75),
                                                  verbose.out = verbose.out,
                                                  nchains = nchains)
        # Swap the prior
        model_pars$SignaturePrior <- swap_SignaturePrior(SignaturePrior = model_pars$SignaturePrior,
                                                         R = sample_preswap$Signatures[1, , ],
                                                         mu = sample_preswap$Mu[1, ],
                                                         S = model_pars$S,
                                                         cutoff = model_pars$cutoff_excluded)
        # Re-initialize the sampler
        init_pars$R <- sample_preswap$Signatures[1, ,]
        init_pars$Theta <- sample_preswap$Weights[1, ,]
        init_pars$mu <- sample_preswap$Mu[1, ]
        Run_mcmc_CompressiveNMF(X, model_pars, init_pars, nsamples,
                                burnin = round(burnin * 0.25),
                                verbose.out = verbose.out, nchains = nchains)
        #CompressiveNMF_cpp(X, nonzero_ids = model)
      } else {
        Run_mcmc_CompressiveNMF(X, model_pars, init_pars, nsamples, burnin,
                                verbose.out = verbose.out, nchains = nchains)
      }
    }
    parallel::stopCluster(cl)
  } else {
    mcmc_out <- lapply(1:nchains, function(i) {
      if(verbose){
        verbose.out <- paste0("chain_status_", i,".txt")
      } else {
        verbose.out <- NULL
      }
      # Step 1 - initialize the sampler
      init_pars <- Initialize_CompressiveNMF(X, model_pars)
      # Step 2 - run the sampler
      if(swap_prior & !is.null(S)){
        # Run the sampler for 3/4, and the perform the swap
        sample_preswap <- Run_mcmc_CompressiveNMF(X,
                                                  model_pars,
                                                  init_pars,
                                                  nsamples = 1,
                                                  burnin = round(burnin * 0.75),
                                                  verbose.out = verbose.out,
                                                  nchains = nchains)
        # Swap the prior
        model_pars$SignaturePrior <- swap_SignaturePrior(SignaturePrior = model_pars$SignaturePrior,
                                                         R = sample_preswap$Signatures[1, , ],
                                                         mu = sample_preswap$Mu[1, ],
                                                         S = model_pars$S,
                                                         cutoff = model_pars$cutoff_excluded)
        # Re-initialize the sampler
        init_pars$R <- sample_preswap$Signatures[1, ,]
        init_pars$Theta <- sample_preswap$Weights[1, ,]
        init_pars$mu <- sample_preswap$Mu[1, ]
        Run_mcmc_CompressiveNMF(X, model_pars, init_pars, nsamples,
                                burnin = round(burnin * 0.25),
                                verbose.out = verbose.out, nchains = nchains)
        #CompressiveNMF_cpp(X, nonzero_ids = model)
      } else {
        Run_mcmc_CompressiveNMF_randa(X, model_pars, init_pars, nsamples, burnin,
                                verbose.out = verbose.out, nchains = nchains)
      }
    })
  }



  # Select the chain with the highest log posterior as solution
  post <-lapply(mcmc_out, function(x) postprocess_mcmc_out(x, cutoff_excluded))

  # Calculate Posterior means
  logposterior <- sapply(1:nchains, function(i) post[[i]]$logpost)
  id_best <- which.max(logposterior)

  # Record wall time of execution
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]

  # Final output
  out <- list(
    Signatures = post[[id_best]]$Signatures,
    Weights = post[[id_best]]$Weights,
    RelWeights = post[[id_best]]$RelWeights,
    selected_chain = id_best,
    time = time,
    mcmc_out = mcmc_out)
  class(out) <- "CompressiveNMF"
  return(out)
}


# RunMCMC CompressiveNMF
#'@export
Run_mcmc_CompressiveNMF_randa <- function(X, model_pars,
                                          init_pars, nsamples, burnin,
                                          verbose.out, nchains, ...){

  # Matrix dimension
  I <- nrow(X); J <- ncol(X); Ktot <- ncol(model_pars$SignaturePrior)

  # Unpack model parameters
  a <- model_pars$a
  random_a <- model_pars$random_a
  alpha <- model_pars$alpha
  epsilon <- model_pars$epsilon
  a0 <- model_pars$a0
  b0 <- model_pars$b0
  c0 <- model_pars$c0
  d0 <- model_pars$d0
  cutoff_excluded <- model_pars$cutoff_excluded
  SignaturePrior <- model_pars$SignaturePrior
  nonzero_ids <- model_pars$nonzero_ids
  verbose <- model_pars$verbose
  betah <- model_pars$betah
  S <- model_pars$S

  # Unpack Prior initialization
  R <- init_pars$R
  Theta <- init_pars$Theta
  mu <- init_pars$mu

  # Store the output
  SIGN <- array(dim = c(nsamples, I, Ktot))
  Ysums <- array(dim = c(nsamples, Ktot, J))
  THETA <- array(dim = c(nsamples, Ktot, J))
  MU <- matrix(nrow = nsamples, ncol = Ktot)
  ASamples <- rep(NA, nsamples)
  LOGPOST <- rep(NA, nsamples)

  # Print the message, if needed
  verbose_step <- round((nsamples + burnin)/20)
  if(verbose & !is.null(verbose.out)){
    cat("Sampling \n", file = verbose.out)
  }

  if(verbose){
    pb <- txtProgressBar(style=3)
  }

  for(iter in 1:(nsamples + burnin)) {
    if(verbose){
      if(nchains == 1){
        setTxtProgressBar(pb, iter/(nsamples + burnin))
      } else {
        if(iter %% verbose_step == 0){
          message <- paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%] \n")
          cat(message)
          if(!is.null(verbose.out)){
            cat(message, file = verbose.out, append = TRUE)
          }
        }
      }
    }
    #------------------------------ 1. Sample the latent variables from multinomial
    Y <- sample_Y(X = X, R = R, Theta = Theta, nonzero_ids = nonzero_ids)
    #------------------------------ 2. Sample the weights
    shape_mat <- a + apply(Y, c(3, 2), sum)
    rate_mat <- 1 + matrix(a / mu) [, rep(1, J)]
    Theta <- sample_weights(shape_mat, rate_mat) + 1e-12
    #------------------------------ 3. Sample the signatures
    # if(swap_prior & !is.null(S)){
    #   if(iter == round(2/3 * burnin))
    #     # Swap the prior to correct for mis-alignment in the sampler. This is done
    #     SignaturePrior <- swap_SignaturePrior(SignaturePrior, K, R, mu, S, cutoff_excluded)
    # }
    Alpha <- SignaturePrior + apply(Y, c(1, 3), sum)
    R <- sample_signatures(Alpha)

    #------------------------------ 4. Sample the global column mean
    if(random_a){
      a0 <- J * a + 1
      b0 <- epsilon * a * J
    }
    mu <- 1/rgamma(Ktot, rep(a0 + J * a, Ktot), b0 + a * rowSums(Theta))
    names(mu) <- colnames(SignaturePrior)

    #------------------------------ OPTIONAL - sample a
    if(random_a){
      a <- sample_FullCond_a(nsamples = 1, Mu = mu, Theta = Theta,
                             epsilon = epsilon, c0 = c0, d0 = d0, a_start = a)
    }

    #------------------------------ 5. Store the output
    if(iter > burnin) {
      Ysums[iter - burnin, ,] <-  apply(Y, c(3, 2), sum)
      SIGN[iter - burnin, ,] <- R
      THETA[iter - burnin, ,] <- Theta
      MU[iter - burnin, ] <- mu
      ASamples[iter - burnin] <- a
      LOGPOST[iter - burnin] <- calculate_logPosterior(X, R, Theta, mu, a, a0, b0, SignaturePrior)
    }
  }
  if(verbose){
    close(pb)
  }

  # Add dimension names
  colnames(MU) <- colnames(SignaturePrior)
  dimnames(SIGN) <- list(NULL, rownames(X), colnames(SignaturePrior))
  dimnames(THETA) <- list(NULL, colnames(SignaturePrior), colnames(X))
  dimnames(Ysums) <- list(NULL, colnames(SignaturePrior), colnames(X))

  # Return the values for each chain
  return(list(Signatures = SIGN, Weights = THETA, Mu = MU, Ysums = Ysums,
              logposterior = LOGPOST, a = ASamples, init_pars = init_pars))
}

