// Re-run SUBCLONAL_INTERPRETATION from a finished tumourevo results folder.
// Run from the repository root so that nextflow.config (containers, interpreter_* params, publishDir) is loaded:
//   nextflow run rerun_subclonal_interpretation.nf -profile docker --results /path/to/outdir --outdir /path/to/new_outdir
// Optional: --patients P1,P2   --viber_heuristic true   --tools mobster,pyclone-vi,viber,sigprofiler   --interpreter_<param> <value>

include { SUBCLONAL_INTERPRETATION } from './modules/local/subclonal_interpretation/main'

workflow {
    if (!params.containsKey('results') || !params.results) {
        error "Missing --results: path to the output folder of a finished tumourevo run"
    }
    def res = params.results
    def wanted = params.containsKey('patients') && params.patients ? params.patients.toString().tokenize(',') : null
    def viber_fit = params.viber_heuristic ? '*_viber_best_heuristic_fit.rds' : '*_viber_best_fit.rds'

    channel.fromPath("${res}/assign_signature/prepare_cluster/*/*", type: 'dir')
        .filter { dir -> !wanted || dir.name in wanted }
        .map { dir ->
            def dataset = dir.parent.name
            def patient = dir.name
            def meta = [dataset: dataset, patient: patient, id: "${dataset}_${patient}"]

            def tables  = files("${res}/assign_signature/prepare_cluster/${dataset}/${patient}/*.txt").sort { f -> f.name }
            def results = files("${res}/assign_signature/assign_cluster/${dataset}/${patient}/*", type: 'dir').sort { f -> f.name }
            def fits    = files("${res}/subclonal_deconvolution/pyclonevi/${dataset}/${patient}/*_best_fit.txt") +
                          files("${res}/subclonal_deconvolution/viber/${dataset}/${patient}/${viber_fit}") +
                          files("${res}/subclonal_deconvolution/mobster/${dataset}/${patient}/*/*_mobster_best_fit.rds")
            def joint   = files("${res}/formatter/cnaqc2tsv/${dataset}/${patient}/*_joint_table.tsv")

            log.info "${meta.id}: ${tables.size()} tables, ${results.size()} SigProfiler dirs, ${fits.size()} fits, ${joint.size()} joint table(s)"
            [meta, tables, results, fits.sort { f -> f.name }, joint]
        }
        .filter { row -> row[1] && row[2] }
        .set { ch_input }

    SUBCLONAL_INTERPRETATION(ch_input)
}
