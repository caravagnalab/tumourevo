#!/usr/bin/env python3

import pandas as pd
import os
import shutil
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
#from SigProfilerMatrixGenerator import install as genInstall


#import input data
joint_table = "/u/cdslab/kdavydzenka/mutationsTable.tsv"
input_data = pd.read_csv(joint_table, sep = '\t')
input_path_sigprof = "/u/cdslab/kdavydzenka/CLL/input_multisample/"

output_path = "output/SBS/CLL.SBS96.all"
output_folder_sigprof = "results/SBS96/"

#input data preprocessing
def input_processing(data):
    new_columns = {'Project': "CLL", 'Genome': 'GRCh37', 'Type': "SOMATIC", 'mut_type': "SNP"}
    df = data.assign(**new_columns)
    df['chr'] = df['chr'].astype(str).str[3:]
    df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
    df["ID"] = df["Sample"]
    df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
    #df = df.style.hide(axis='index') #omit row indexes
    return df

input_data = input_processing(input_data)

#saving input matrix to txt
input_data.to_csv('/u/cdslab/kdavydzenka/CLL/input_multisample/input_data.txt', sep='\t', index=False, header=True)

#Install your desired reference genome
#genInstall.install('GRCh37', rsync=False, bash=True)


#mutation's counts matrix generation
input_matrix = matGen.SigProfilerMatrixGeneratorFunc(project = "CLL", 
                                                     reference_genome = "GRCh37", 
                                                     path_to_input_files = input_path_sigprof)


#SBS_file = "/u/cdslab/kdavydzenka/CLL/input_multisample/output/SBS/CLL.SBS96.all"
output_path = "output/SBS/CLL.SBS96.all"

# Perform model fitting
sig.sigProfilerExtractor(input_type = "matrix", 
                         output = "results", 
                         input_data = input_path_sigprof+output_path,  
                         exome = False,
                         minimum_signatures = 1,
                         maximum_signatures = 10,
                         nmf_replicates = 100, 
                         resample = True,
                         matrix_normalization = "gmm", 
                         seeds= "random",
                         nmf_init = "random", 
                         min_nmf_iterations = 10000, 
                         max_nmf_iterations = 1000000, 
                         nmf_test_conv = 10000, 
                         nmf_tolerance = 1e-15,
                         cpu = 6,
                         gpu = False,
                         cosmic_version = 3.4,
                         make_decomposition_plots = True, 
                         collapse_to_SBS96 = True, 
                         get_all_signature_matrices = True,
                         export_probabilities = True)


#Generate a stacked bar plot showing activities in individuals
plotActivity(activity_file, 
             output_file = "Activity_in_samples.pdf", 
             bin_size = 50, 
             log = False)


#save the output results
dest_dir = "signature_deconvolution/Sigprofiler/CLL/"
source_dir = "results"
shutil.copytree(source_dir, dest_dir, dirs_exist_ok=True)
