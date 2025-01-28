from SigProfilerExtractor import sigpro as sig
import os

overd_list = [0, 0.15]

for overd in overd_list:
  for i in range(1,21):
    
    # Run SigProfilerExtractor on the simulated data
    project_name = "~/CompressiveNMF/output/sparse_indel_simulation/Indel_theta_50_overd_{}/SigProfiler/output/sim_{}".format(overd, i)  # Output directory name
    input_type = "matrix"             # Specify the input type
    input_data = "~/CompressiveNMF/output/sparse_indel_simulation/Indel_theta_50_overd_{}/SigProfiler/data/X_{}.txt".format(overd, i) # Load data
    
    # Call the sigprofiler function
    sig.sigProfilerExtractor(input_type, os.path.expanduser(project_name), os.path.expanduser(input_data), 
                             minimum_signatures=2, maximum_signatures=20, nmf_replicates=10)



