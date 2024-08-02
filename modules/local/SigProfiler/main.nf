process SIG_PROFILER {
    tag "$meta.id"
    container = 'docker://katiad/sigprofiler:dev1'

    input:
       //tuple val(datasetID), val(patientID), val(sampleID), path(joint_table) //from formatter output
       tuple val(meta), path(joint_table)

    output:
       tuple val(datasetID), path("signature_deconvolution/Sigprofiler/$datasetID/results/SBS96/SBS96_selection_plot.pdf"),
       path("signature_deconvolution/Sigprofiler/$datasetID/results/SBS96/Suggested_Solution/"), 
       path("signature_deconvolution/Sigprofiler/$datasetID/results/SBS96/Samples.txt")
      
    script:
    
      def args                              = task.ext.args                                 ?: ''
      def reference_genome                  = args!='' && args.reference_genome             ? "$args.reference_genome" : ""
      def exome                             = args!='' && args.exome                        ? "$args.exome" : ""
      def bed_file                          = args!='' && args.bed_file                     ? "$args.bed_file" : ""
      def chrom_based                       = args!='' && args.chrom_based                  ? "$args.chrom_based" : ""
      def plot                              = args!='' && args.plot                         ? "$args.plot" : ""
      def tsb_stat                          = args!='' && args.tsb_stat                     ? "$args.tsb_stat" : ""
      def seqInfo                           = args!='' && args.seqInfo                      ? "$args.seqInfo" : ""
      def cushion                           = args!='' && args.cushion                      ? "$args.cushion" : ""
      def volume                            = args!='' && args.volume                       ? "$args.volume" : ""
      def input_type                        = args!='' && args.input_type                   ? "$args.input_type" : ""
      def context_type                      = args!='' && args.context_type                 ? "$args.context_type" : ""
      def minimum_signatures                = args!='' && args.minimum_signatures           ? "$args.minimum_signatures" : ""
      def maximum_signatures                = args!='' && args.maximum_signatures           ? "$args.maximum_signatures" : ""
      def nmf_replicates                    = args!='' && args.nmf_replicates               ? "$args.nmf_replicates" : ""
      def resample                          = args!='' && args.resample                     ? "$args.resample" : ""
      def seeds                             = args!='' && args.seeds                        ? "$args.seeds" : ""
      def matrix_normalization              = args!='' && args.matrix_normalization         ? "$args.matrix_normalization" : ""
      def nmf_init                          = args!='' && args.nmf_init                     ? "$args.nmf_init" : "random"
      def min_nmf_iterations                = args!='' && args.min_nmf_iterations           ? "$args.min_nmf_iterations" : ""
      def max_nmf_iterations                = args!='' && args.max_nmf_iterations           ? "$args.max_nmf_iterations" : ""
      def nmf_test_conv                     = args!='' && args.nmf_test_conv                ? "$args.nmf_test_conv" : ""
      def nmf_tolerance                     = args!='' && args.nmf_tolerance                ? "$args.nmf_tolerance" : ""
      def cpu                               = args!='' && args.cpu                          ? "$args.cpu" : ""
      def gpu                               = args!='' && args.gpu                          ? "$args.gpu" : ""
      def batch_size                        = args!='' && args.batch_size                   ? "$args.batch_size" : ""
      def stability                         = args!='' && args.stability                    ? "$args.stability" : ""
      def min_stability                     = args!='' && args.min_stability                ? "$args.min_stability" : ""
      def combined_stability                = args!='' && args.combined_stability           ? "$args.combined_stability" : ""   
      def cosmic_version                    = args!='' && args.cosmic_version               ? "$args.cosmic_version" : ""
      def de_novo_fit_penalty               = args!='' && args.de_novo_fit_penalty          ? "$args.de_novo_fit_penalty" : ""
      def nnls_add_penalty                  = args!='' && args.nnls_add_penalty             ? "$args.nnls_add_penalty" : ""
      def nnls_remove_penalty               = args!='' && args.nnls_remove_penalty          ? "$args.nnls_remove_penalty" : ""
      def initial_remove_penalty            = args!='' && args.initial_remove_penalty       ? "$args.initial_remove_penalty" : ""
      def refit_denovo_signatures           = args!='' && args.refit_denovo_signatures      ? "$args.refit_denovo_signatures" : ""
      def make_decomposition_plots          = args!='' && args.make_decomposition_plots     ? "$args.make_decomposition_plots" : ""
      def collapse_to_SBS96                 = args!='' && args.collapse_to_SBS96            ? "$args.collapse_to_SBS96" : ""
      def get_all_signature_matrices        = args!='' && args.get_all_signature_matrices   ? "$args.get_all_signature_matrices" : ""
      def export_probabilities              = args!='' && args.export_probabilities         ? "$args.export_probabilities": ""   

    
      """
      #!/usr/bin/env python3
     
      import os
      import shutil
      import pandas as pd
      from SigProfilerExtractor import sigpro as sig
      from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
      #from SigProfilerMatrixGenerator import install as genInstall
  
    
      #if os.path.exists(output_path):
       #shutil.rmtree(output_path)
    
      os.mkdir("$datasetID")
      os.mkdir("signature_deconvolution/Sigprofiler/$datasetID/")
   
      input_path = "$datasetID/"
      output_path = "output/SBS/CLL.SBS96.all"
      output_folder_sigprof = "results/SBS96/"

      input_data = pd.read_csv("$joint_table", sep = "\t")
      #'/orfeo/LTS/LADE/LT_storage/lvaleriani/nextflow_modules/work/28/a52d0fb3d52c3a96126331f9b9226c/joint_table.tsv'
    
      #input data preprocessing
      def input_processing(data):
         new_columns = {'Project': "$datasetID", 'Genome': '$reference_genome', 'Type': "SOMATIC", 'mut_type': "SNP"}
         df = data.assign(**new_columns)
         df['chr'] = df['chr'].astype(str).str[3:]
         df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
         df["ID"] = df["Sample"]
         df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
         return df
    
      input_data = input_processing(input_data)

      #saving input matrix to txt
      input_data.to_csv("input_path/input_data.txt", sep="\t", index=False, header=True)

      #download desired reference genome
      #genInstall.install('$reference_genome', rsync=False, bash=True)

      #mutation's counts matrix generation
      input_matrix = matGen.SigProfilerMatrixGeneratorFunc(
              project = "$datasetID", 
              reference_genome = "$reference_genome", 
              path_to_input_files = input_path,
              exome = bool("$exome"),
              bed_file = eval("$bed_file"),
              chrom_based = bool("$chrom_based"),
              plot = bool("$plot"),
              tsb_stat = bool("$tsb_stat"),
              seqInfo = bool("$seqInfo),
              cushion = int("$cushion))

      # Perform model fitting
      sig.sigProfilerExtractor(input_type = "$input_type", 
                               out_put = "results", 
                               input_data = input_path+output_path,  
                               context_type = "$context_type",  
                               exome = "$exome",
                               minimum_signatures = int("$minimum_signatures"),  
                               maximum_signatures = int("$maximum_signatures"), 
                               nmf_replicates = int("$nmf_replicates"),
                               resample = bool("$resample"),
                               matrix_normalization = "$matrix_normalization", 
                               seeds= "$seeds",
                               nmf_init = "$nmf_init", 
                               min_nmf_iterations = int("$min_nmf_iterations"), 
                               max_nmf_iterations = int("$max_nmf_iterations"),
                               nmf_test_conv = int("$nmf_test_conv"), 
                               nmf_tolerance = float("$nmf_tolerance"), 
                               cpu = int("$cpu"),
                               gpu = bool("$gpu"),
                               batch_size = int("$batch_size"),
                               stability = float("$stability"),
                               min_stability = float("$min_stability"),
                               combined_stability = float("$combined_stability"),
                               cosmic_version = float("$cosmic_version"),
                               de_novo_fit_penalty = float("$de_novo_fit_penalty"),
                               nnls_add_penalty = float("$nnls_add_penalty"),
                               nnls_remove_penalty = float("$nnls_remove_penalty"),
                               initial_remove_penalty = float("$initial_remove_penalty"),
                               refit_denovo_signatures = bool("$refit_denovo_signatures"),
                               make_decomposition_plots = bool("$make_decomposition_plots"), 
                               collapse_to_SBS96 = bool("$collapse_to_SBS96"), 
                               get_all_signature_matrices = bool("$get_all_signature_matrices"),
                               export_probabilities = bool("$export_probabilities"))
    
    

      #save the output results

      source_dir = "output_folder_sigprof"
      dest_dir = "signature_deconvolution/Sigprofiler/$datasetID/"
      shutil.copytree(source_dir, dest_dir)
   """
}
