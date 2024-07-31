//
// LIFTER SUB-WORKFLOW
//

include { GET_POSITIONS } from "../../../modules/local/get_positions/main"
include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { JOIN_POSITIONS } from "../../../modules/local/join_positions/main"


workflow LIFTER {
    take:
        data
        fasta

    main:
        out = Channel.empty()

        rds = data.map{ meta, rds, bam, bai -> 
            meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            sample = meta.tumour_sample
            [meta.subMap('dataset', 'patient', 'id', 'normal_sample'), rds, bam, bai, sample]}
            | groupTuple

        GET_POSITIONS(rds)
        all_pos = GET_POSITIONS.out.all_pos.transpose().map{meta, rds, sample -> 
            [meta, rds]
        }

        bed = GET_POSITIONS.out.bed.transpose().map{ meta, bed, bam, bai, sample ->
            meta = meta + [tumour_sample: "${sample}".toInteger()]
            meta = meta + [id: "${meta.dataset}_${meta.patient}_${meta.tumour_sample}"]
            [meta, bam, bed]
        }

        BCFTOOLS_MPILEUP(bed, fasta, false)
        tmp_data = data.map{ meta, rds, bam, bai -> 
            [meta.subMap('dataset', 'patient', 'id', 'normal_sample', 'tumour_sample'), rds]
        }

        join = tmp_data.join(BCFTOOLS_MPILEUP.out.vcf)
        out = JOIN_POSITIONS(join, all_pos)

    emit:
        out 

}
