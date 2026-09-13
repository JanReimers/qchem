target_sources(qcGaussian_BS PRIVATE
Imp/BasisSet.C
Imp/IrrepBasisSet.C
)
target_sources(qcGaussian_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    BasisSet.C
)
target_include_directories(qcGaussian_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
