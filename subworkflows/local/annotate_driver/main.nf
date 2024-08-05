//
// DRIVER_ANNOTATION SUB-WORKFLOW
// 
// Should work with also with prebuilt RefCDS
//
include { BUILD_REFERENCE } from '../../../modules/local/build_reference/main.nf'
include { DNDSCV } from '../../../modules/local/dndscv/main.nf'


workflow DRIVER_ANNOTATION {
    take:        
        rds
        driver_list
        cds
        genome
    
    main:

        if (params.dndscv_refcds_rda) {
            println(params.dndscv_refcds_rda)
            rda = Channel.from(file(params.dndscv_refcds_rda, checkIfExists: true))
            dndscv_ch = rds.combine(driver_list).combine(rda)
        } else {
            // create reference for variant annotation, meta here is not needed
            BUILD_REFERENCE(cds, genome)
            rda = BUILD_REFERENCE.out.dnds_reference
            dndscv_ch = rds.combine(driver_list).combine(rda)
        }
        
        DNDSCV(dndscv_ch) // add dndsCV statistics and columns "known_driver" (based on driver list, default IntoGen) and "potential_driver"

    emit:
        DNDSCV.out.dnds_rds
}
