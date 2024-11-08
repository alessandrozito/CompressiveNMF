""" 
This file runs the python 2.7 method from Brouwer et al (2017) on the 21 breast 
cancer data
"""
import numpy as np
import scipy
import random
import os

# Load the Brouwer data
from BNMTF_ARD.code.models.bnmf_gibbs import bnmf_gibbs

# Set seed
random.seed(0); scipy.random.seed(0); np.random.seed(0)

''' Model settings. '''
nsamples = 1000
burn_in = 3000
iterations = nsamples + burn_in
thinning = 1

init_UV = 'random'
K = 20 # number of latent factors
ARD = True

lambdaU, lambdaV = 0.1, 0.1
alphatau, betatau = 1., 1.
alpha0, beta0 = 1., 1.
hyperparams = { 'alphatau':alphatau, 
                'betatau':betatau, 
                'alpha0':alpha0, 
                'beta0':beta0, 
                'lambdaU':lambdaU, 
                'lambdaV':lambdaV }
                
X = np.genfromtxt(os.path.expanduser("~/CompressiveNMF/R/run_BayesNMF_python/X_21breast.txt"), skip_header=1)[:, 1:]
# Run the method
M = np.ones(X.shape) 
BNMF = bnmf_gibbs(X, M, K, ARD, hyperparams) 
BNMF.initialise(init_UV)
BNMF.run(iterations)
# Save the output
U = BNMF.all_U[burn_in:,:,: ].transpose(2,0,1).reshape(K, -1)
V = BNMF.all_V[burn_in:,:,: ].transpose(2,0,1).reshape(K, -1)
lambdak = BNMF.all_lambdak[burn_in:,:]
tau = BNMF.all_tau[burn_in:]

output_filder = os.path.expanduser("~/CompressiveNMF/output/Application_21brca/BayesNMF_brouwer")

np.savetxt(output_filder + "/Signatures.txt", U.transpose(1, 0))
np.savetxt(output_filder + "/Loadings.txt", V.transpose(1, 0))
np.savetxt(output_filder + "/tau.txt", tau)
np.savetxt(output_filder +  "/lambda.txt", lambdak)
np.savetxt(output_filder + "/times.txt", BNMF.all_times)






