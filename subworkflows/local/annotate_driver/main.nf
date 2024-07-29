//
// DRIVER_ANNOTATION SUB-WORKFLOW
//

include { ANNOTATE_DRIVER } from '../../../modules/local/annotate_driver/main'


workflow DRIVER_ANNOTATION {
    take:
        rds
    
    main:
        ANNOTATE_DRIVER(rds)


    emit:
        ANNOTATE_DRIVER.out.rds

}
