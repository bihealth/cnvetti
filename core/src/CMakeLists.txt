# Create version.cpp file with appropriate version
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/version.cpp.in
                ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

add_executable(
    cnvetti
    main.cpp
    program_options.cpp
    main_coverage.cpp
    main_summaries.cpp
    ${version_file})

target_compile_features(
    cnvetti
    PRIVATE
    cxx_strong_enums)

install(
    TARGETS cnvetti
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
