//
// LIFTER SUB-WORKFLOW
//

include { GET_POSITIONS } from "../../modules/local/get_positions/main"
// include { BCFTOOLS_MPILEUP } from "../../modules/local/mpileup/main"
include { BCFTOOLS_MPILEUP } from '../modules/nf-core/bcftools/mpileup/main'
include { JOIN_POSITIONS } from "../../modules/local/join_positions/main"


workflow LIFTER {
    take:
        rds
        tumor_bam


    main:

        GET_POSITIONS(rds.groupTuple(by: [0,1]))
        BCFTOOLS_MPILEUP(GET_POSITIONS.out.bed.transpose(by: [2,3]), tumor_bam)
        pileup = rds.join(BCFTOOLS_MPILEUP.out.vcf, by: [0,1,2])
        JOIN_POSITIONS(pileup, GET_POSITIONS.out.pos.transpose(by: 2))

    emit:
        JOIN_POSITIONS.out.rds

}
