//
// QC SUB-WORKFLOW
//

include { TINC } from '../../../modules/local/tinc/main'
include { CNAQC } from '../../../modules/local/CNAqc/main'
include { JOIN_CNAQC } from '../../../modules/local/join_CNAqc/main'


workflow QC {
    take: 
        input

    main:

        // TINC(cna, vcf)
        CNAQC(input)
        in_join = CNAQC.out.qc_rds.map{ meta, rds -> 
            meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            sample = meta.tumour_sample
            [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
            | groupTuple
        JOIN_CNAQC(in_join)
    
    emit:
        rds_cnaqc = CNAQC.out.qc_rds
        plot_cnaqc_rds = CNAQC.out.plot_rds
        plot_cnaqc_data = CNAQC.out.plot_pdf_data
        plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

        // plot_rds_tinc = TINC.out.plot_rds
        // rds_tinc = TINC.out.rds
        // pdf_tinc = TINC.out.plot_pdf
        
        rds_join = JOIN_CNAQC.out.rds
}
