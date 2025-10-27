//
// QC SUB-WORKFLOW
//

include { TINC } from '../../../modules/nf-core/tinc/main'
include { CNAQC } from '../../../modules/nf-core/cnaqc/main'
include { JOIN_CNAQC } from '../../../modules/local/join_cnaqc/main'


workflow QC {
    take:
        input

    main:
        TINC(input)
        // contamination = TINC.out.tinc_csv
        //     .splitCsv( header: true )
        //     .map{ meta, csv ->
        //     meta = meta + [id: "${meta.dataset}_${meta.patient}"]
        //     normal_contamination = csv.normal_contamination_flag
        //     [meta.subMap('dataset', 'patient', 'id'), normal_contamination ]}
        //     .unique()
        //     | groupTuple
        //     | map{ meta, normal_contamination ->
        //         normal_contamination = normal_contamination.max()
        //         [ meta, normal_contamination]
        //     }
        input_cnaqc = input.map{meta, cna, snv ->
                def sample =  meta.tumour_sample
                [meta, snv, cna, sample]
        }

        CNAQC(input_cnaqc)
        in_join_cnaqc = CNAQC.out.qc_rds.map{ meta, rds ->
            def sample = meta.tumour_sample
            meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
            | groupTuple

        // out_tinc = in_join_cnaqc.map{  meta, rds, sample, normal_contamination ->
        //    meta = meta + [nc: normal_contamination]
        //        [meta, rds, sample] }
        //        .branch { meta, rds, sample ->
        //                pass: meta.nc == '0'
        //                not_pass: meta.nc == '1'
        //        }
        JOIN_CNAQC(in_join_cnaqc)

    emit:
        rds_cnaqc = CNAQC.out.qc_rds
        plot_cnaqc_data_rds = CNAQC.out.data_plot_rds
        plot_cnaqc_qc_rds = CNAQC.out.qc_plot_rds
        plot_cnaqc_data = CNAQC.out.plot_pdf_data
        plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

        plot_rds_tinc = TINC.out.plot_rds
        rds_tinc = TINC.out.rds
        pdf_tinc = TINC.out.plot_pdf
        csv_tinc = TINC.out.tinc_csv

        join_cnaqc_ALL = JOIN_CNAQC.out.rds_all
        join_cnaqc_PASS = JOIN_CNAQC.out.rds_pass
}
