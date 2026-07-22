//
// QC SUB-WORKFLOW
//

include { TINC } from '../../../modules/nf-core/tinc/main'
include { CNAQC } from '../../../modules/nf-core/cnaqc/main'
include { JOIN_CNAQC } from '../../../modules/local/join_cnaqc/main'
include { QC_BY_CHR } from '../../../modules/local/cnaqc_by_chr/main'
include { JOIN_CNAQC_BY_CHR } from '../../../modules/local/join_cnaqc_by_chr/main'


workflow QC {
    take:
        input

    main:
        ch_versions = Channel.empty()

        plot_rds_tinc = null
        rds_tinc = null
        pdf_tinc = null
        csv_tinc = null

        plot_cnaqc_by_chr_qc = null
        plot_cnaqc_by_chr_qc_pdf = null
        rds_cnaqc_by_chr = null

        if (params.tools && params.tools.split(',').contains('tinc')) {
            TINC(input)
            ch_versions = ch_versions.mix(TINC.out.versions_tinc)
            plot_rds_tinc = TINC.out.plot_rds
            rds_tinc = TINC.out.rds
            pdf_tinc = TINC.out.plot_pdf
            csv_tinc = TINC.out.tinc_csv
        }

        input_cnaqc = input.map{meta, cna, snv ->
                def sample =  meta.tumour_sample
                [meta, snv, cna, sample]
        }

        CNAQC(input_cnaqc)
        ch_versions = ch_versions.mix(CNAQC.out.versions)
        
        if(params.qc_chr == true) {

            in_qc_by_chr = CNAQC.out.qc_rds.map{ meta, rds ->
                def sample = meta.tumour_sample
                [meta, rds, sample]
            }
            QC_BY_CHR(in_qc_by_chr)

            plot_cnaqc_by_chr_qc_pdf = QC_BY_CHR.out.plot_pdf_qc_by_chr
            plot_cnaqc_by_chr_qc = QC_BY_CHR.out.plot_qc_by_chr
            rds_cnaqc_by_chr = QC_BY_CHR.out.qc_rds
            ch_versions = ch_versions.mix(QC_BY_CHR.out.versions)
            
            in_join_cnaqc = QC_BY_CHR.out.qc_rds.map{ meta, rds ->
                def sample = meta.tumour_sample
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
                .groupTuple()

            JOIN_CNAQC_BY_CHR(in_join_cnaqc)
            join_cnaqc_ALL = JOIN_CNAQC_BY_CHR.out.rds_all
            join_cnaqc_PASS = JOIN_CNAQC_BY_CHR.out.rds_pass

            ch_versions = ch_versions.mix(JOIN_CNAQC_BY_CHR.out.versions)

        } else {
            in_join_cnaqc = CNAQC.out.qc_rds.map{ meta, rds ->
                def sample = meta.tumour_sample
                meta = meta + [id: "${meta.dataset}_${meta.patient}"]
                [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
                .groupTuple()

            JOIN_CNAQC(in_join_cnaqc)
            join_cnaqc_ALL = JOIN_CNAQC.out.rds_all
            join_cnaqc_PASS = JOIN_CNAQC.out.rds_pass

            ch_versions = ch_versions.mix(JOIN_CNAQC.out.versions)
        }

    emit:
        
        // emits always the same elements, but slightly different based on if the qc is run by chr or not
        rds_cnaqc = CNAQC.out.qc_rds
        plot_cnaqc_data_rds = CNAQC.out.data_plot_rds
        plot_cnaqc_qc_rds = CNAQC.out.qc_plot_rds
        plot_cnaqc_data = CNAQC.out.plot_pdf_data
        plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

        plot_rds_tinc
        rds_tinc
        pdf_tinc
        csv_tinc

        plot_cnaqc_by_chr_qc
        plot_cnaqc_by_chr_qc_pdf
        rds_cnaqc_by_chr

        join_cnaqc_ALL
        join_cnaqc_PASS
        
        versions = ch_versions
}
