# Command-line applications (not tests).

# scfrun: command-line SCF driver for accelerator tuning (no recompile).  Drives qchem::AtomCalculation
# directly (no QchemTester); the RelativeHF/DHFError oracle helpers come from qchem.PeriodicTable (R2.20)
# -- no test module imported.
add_executable(scfrun scfrun.C)
target_compile_definitions(scfrun PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET scfrun PROPERTY CXX_STANDARD 20)
target_include_directories(scfrun PRIVATE
    ${QCHEM_BSPLINE_INCLUDE}
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_TABULATE_INCLUDE}
)
target_link_libraries(scfrun PUBLIC qchem qcMath qcCommon blas)

# valgen: command-line valence-basis generator (GPWPlan1 step 1).  A thin arg-parser over the
# qchem.ValenceBasisGen library; dogfoods qchem.Reporting for its console output.
add_executable(valgen valgen.C)
target_compile_definitions(valgen PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET valgen PROPERTY CXX_STANDARD 20)
target_include_directories(valgen PRIVATE
    ${QCHEM_BSPLINE_INCLUDE}
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_TABULATE_INCLUDE}
)
target_link_libraries(valgen PUBLIC qchem qcMath qcCommon blas)

# gpwprobe: the GPW campaign INSTRUMENTS (supercell ladder, k sweep, NaF GDM-smear probe, the V2.6 Becke
# ladders, the MnO AFM-II campaign run) -- promoted out of IntegrationTests/GPW_SCF_UT.C's DISABLED_ tests
# (doc/TestSuitePlan.md §8).  Shares the gates' harness (qcGPW_Harness) so a probe and its gate build the
# same cell the same way.
add_executable(gpwprobe gpwprobe.C)
target_compile_definitions(gpwprobe PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET gpwprobe PROPERTY CXX_STANDARD 20)
target_include_directories(gpwprobe PRIVATE
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_BLAZE_DIR}
)
target_link_libraries(gpwprobe PUBLIC qcGPW_Harness qchem qcMath qcCommon blas)

# rundeck: run a periodic SCF from an input DECK (D-ENV step 6): structure NAME + options, --set overrides, a numbered resolved
# revision written before the run starts.  Thin: all of it is qchem.Deck.
# qcFrontEnd: the RUN-CONTROL FRONT END rundeck and ITMain share (qchem.FrontEnd.RunControl: signals, the pid-keyed control directory,
# the heartbeat).  Front-end code only: libqchem never links it (the library sees qchem.RunControl's abstract face alone).
add_library(qcFrontEnd)
target_sources(qcFrontEnd PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES FrontEnd/RunControlFrontEnd.C)
target_compile_definitions(qcFrontEnd PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET qcFrontEnd PROPERTY CXX_STANDARD 20)
target_link_libraries(qcFrontEnd PUBLIC qcCommon)

add_executable(rundeck rundeck.C)
target_compile_definitions(rundeck PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET rundeck PROPERTY CXX_STANDARD 20)
target_include_directories(rundeck PRIVATE
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_BLAZE_DIR}
)
target_link_libraries(rundeck PUBLIC qcFrontEnd qchem qcMath qcCommon blas)

# gth2upf: write a QE UPF for one of our GTH pseudopotentials with PP_CHI = our pseudo-atom's orbitals, so
# pw.x/hp.x (the U oracle, DFT+U increment 3) run the same PP and the same `atomic` Hubbard projector.
add_executable(gth2upf gth2upf.C)
target_compile_definitions(gth2upf PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET gth2upf PROPERTY CXX_STANDARD 20)
target_include_directories(gth2upf PRIVATE ${QCHEM_BLAZE_DIR} ${QCHEM_JSON_INCLUDE})
target_link_libraries(gth2upf PUBLIC qchem qcMath qcCommon blas)

# q2probe: the Q2 molecular Hartree-route table's driver (lane S 2026-11-25): a molecule deck converged, its wall split setup / SCF.
add_executable(q2probe q2probe.C)
target_compile_definitions(q2probe PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET q2probe PROPERTY CXX_STANDARD 20)
target_include_directories(q2probe PRIVATE ${QCHEM_JSON_INCLUDE} ${QCHEM_BLAZE_DIR})
target_link_libraries(q2probe PUBLIC qchem qcMath qcCommon blas)
