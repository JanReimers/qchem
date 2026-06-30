# qcUmbrella: the `qchem` umbrella module -- a thin convenience aggregator above the whole stack.
# Deliberately its OWN library so the umbrella never re-exports a same-library sibling (which CMake's
# C++-module scanner would read as a bogus cycle -- the qcMath-split lesson).  It re-exports modules
# that all live in already-linked libraries (the qchem aggregate + qcMath + qcStructure ...), so there
# is no module back-edge.  No #includes here -> no extra include dirs needed.
add_library(qcUmbrella)
target_sources(qcUmbrella PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Umbrella.C
)
set_property(TARGET qcUmbrella PROPERTY CXX_STANDARD 20)
target_compile_options(qcUmbrella PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcUmbrella PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcUmbrella PUBLIC qchem)
