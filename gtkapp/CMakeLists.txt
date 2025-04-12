cmake_minimum_required(VERSION 3.5.0)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GTKMM REQUIRED gtkmm-4.0)

#file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>: ../../oml/Debug>
    $<$<CONFIG:Release>: ../../oml/Release>
)
add_executable(qchem-gtk qchem-gtk.cpp)

target_compile_options(qchem-gtk PRIVATE -Wall -pg -fopenmp ${GTKMM_CFLAGS_OTHER}
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_link_options(qchem-gtk PRIVATE -fopenmp ${GTKMM_LDFLAGS_OTHER})
target_include_directories(qchem-gtk PRIVATE ../include ../../oml/include ${GTKMM_INCLUDE_DIRS})
target_link_libraries(qchem-gtk qchem ${GTKMM_LIBRARIES})
