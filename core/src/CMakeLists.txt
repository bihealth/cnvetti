# Create version.cpp file with appropriate version
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/version.cpp.in
                ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

find_library(HTS_LIBRARY NAMES hts)

add_executable(
    cnvetti
    contig_selection.cpp
    main.cpp
    main_background.cpp
    main_coverage.cpp
    main_normalize.cpp
    main_segment.cpp
    histo_stats.cpp
    program_options.cpp
    segmentation_haar.cpp
    ${version_file})

target_link_libraries(cnvetti ${HTS_LIBRARY})

target_compile_features(
    cnvetti
    PRIVATE
    cxx_strong_enums)

install(
    TARGETS cnvetti
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
