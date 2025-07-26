target_sources(qchem PRIVATE Imp/SCFIterator.C)
target_sources(qchem PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    SCFIterator.C
    SCFParams.C
)

