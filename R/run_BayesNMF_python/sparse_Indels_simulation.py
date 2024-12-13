""" 
This file runs the python 2.7 method from Brouwer et al (2017) on simulated data
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
burn_in = 2000
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

'''Load the data from all repositories'''
overd_list = [0, 0.15]
theta = 50
simulation_dir = "~/CompressiveNMF/output/sparse_indel_simulation/"

# Iterate across all scenarios un background
for overd in overd_list:
    file_path = "{}Indel_theta_{}_overd_{}/BayesNMF_brouwer/".format(simulation_dir, theta, overd)
    print(file_path)
    # Set imput and output files
    input_folder = os.path.expanduser(file_path + "data/")
    output_filder = os.path.expanduser(file_path + "output/")
    # Iterate across datasets
    for i in range(1, 21):
      # Load data
      X = np.genfromtxt(input_folder + "X_" + str(i) + ".txt",  skip_header=1)[:, 1:]
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
      
      np.savetxt(output_filder + "sim_" + str(i) + "/Signatures.txt", U.transpose(1, 0))
      np.savetxt(output_filder + "sim_" + str(i) + "/Loadings.txt", V.transpose(1, 0))
      np.savetxt(output_filder + "sim_" + str(i) + "/tau.txt", tau)
      np.savetxt(output_filder + "sim_" + str(i) + "/lambda.txt", lambdak)
      np.savetxt(output_filder + "sim_" + str(i) + "/times.txt", BNMF.all_times)
        
        

