//
// MUTATIONAL SIGNATURES DECONVOLUTION WORKFLOW
//

include { FORMATTER as FORMATTER_RDS} from "../../../subworkflows/local/formatter/main"
include { SPARSE_SIGNATURES } from "../../../modules/local/SparseSignatures/main"
include { SIGPROFILER } from "../../modules/local/SigProfiler/SigProfiler/main"
include { DOWNLOAD_GENOME_SIGPROFILER } from "../../modules/local/SigProfiler/download/main"


workflow SIGNATURE_DECONVOLUTION {
    take:
    meta 
    joint_table

    main:
    plot_pdf = null
    plot_rds = null
    signatures_nmfOut = null 
    bestConf = null
    sign_cv = null
    Sigprofiler_out = null


    if (params.tools && params.tools.split(',').contains('sparsesignatures')) {
        out_sparse = FORMATTER_RDS(joint_table, "rds")
        SPARSE_SIGNATURES(out_sparse.groupTuple(by: 0)) // run SparseSignatures
        
        plot_pdf = SPARSE_SIGNATURES.out.signatures_plot_pdf
        plot_rds = SPARSE_SIGNATURES.out.signatures_plot_rds
        signatures_nmfOut = SPARSE_SIGNATURES.out.signatures_nmfOut_rds
        bestConf = SPARSE_SIGNATURES.out.signatures_bestConf_rds
        sign_cv = SPARSE_SIGNATURES.out.signatures_cv_rds
    } 
  

    if (params.tools && params.tools.split(',').contains('sigprofiler')) {
        // Check if we should download SigProfiler genome
        if (params.download_genome_sigprofiler) {
            
            genome = DOWNLOAD_GENOME_SIGPROFILER(meta).genome
        
        } else {
           
            // Use pre-installed genome path 
            def installed_genome_path = "signature_deconvolution/Sigprofiler/genome/"

            // Check if the pre-installed genome path exists
            if (!file(installed_genome_path).exists()) {
                error "The pre-installed genome path $installed_genome_path does not exist! Please install the genome or set download_genome_sigprofiler to true."
            } else {
                genome = installed_genome_path
            }
        }
            
       
        out_sigprof = FORMATTER_RDS(joint_table, "rds")
        Sigprofiler_out = SIGPROFILER(out_sigprof, genome) // run SigProfiler
      
        
    }

    emit:
    plot_pdf
    plot_rds
    signatures_nmfOut
    bestConf
    sign_cv
    Sigprofiler_out
    
}
