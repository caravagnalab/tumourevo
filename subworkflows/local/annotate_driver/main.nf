//
// DRIVER_ANNOTATION SUB-WORKFLOW
// 

include { BUILD_REFERENCE } from '../../../modules/local/annotate_driver/main'
include { DNDSCV } from '../../../modules/local/annotate_driver/main'


workflow DRIVER_ANNOTATION {
    take:
        rds
        driver_list
    
    main:
        ref = BUILD_REFERENCE() // create reference for variant annotation
        DNDSCV(rds, driver_list, ref) // add dndsCV statistics and columns "known_driver" (based on driver list, default IntoGen) and "potential_driver"

    emit:
        ANNOTATE_DRIVER.out.rds

}
