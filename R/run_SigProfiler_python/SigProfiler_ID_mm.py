
from SigProfilerExtractor import sigpro as sig

# Run SigProfilerExtractor on the MM data
project_name = "output/Application_ID_mm/sigprofiler_out/"  # Output directory name
input_type = "matrix"             # Specify the input type
input_data = "data/Indels_mutliple_myeloma.txt" 

# Call the sigprofiler function
sig.sigProfilerExtractor(input_type, project_name, input_data, 
                         minimum_signatures=2, maximum_signatures=15, 
                         nmf_replicates=10)



