add_library(fabm_models_memg OBJECT
            memg_model_library.F90
            bio_fennel/memg_bio_fennel_model_library.F90
            bio_fennel/1p1z.F90
            bio_fennel/2p2z.F90
            bio_fennel/2p2z_sink.F90
            bio_fennel/light.F90
            bio_fennel/oxygen.F90
           )

source_group(bio_fennel REGULAR_EXPRESSION bio_fennel/)

add_dependencies(fabm_models_memg fabm_base)

