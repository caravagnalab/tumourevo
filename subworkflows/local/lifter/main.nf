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
            [meta, rds]
        }

        tumour_bam = data.map{ meta, rds, bam, bai ->
            [meta, bam, bai]
        }

        GET_POSITIONS(rds.groupTuple(by: [0,1]))
        //BCFTOOLS_MPILEUP([meta, tumour_bam, GET_POSITIONS.out.bed.transpose(by: [2,3])], fasta, false)
        //pileup = rds.join(BCFTOOLS_MPILEUP.out.vcf, by: [0,1,2])
        //out = JOIN_POSITIONS(pileup, GET_POSITIONS.out.pos.transpose(by: 2))

    emit:
        out 

}
