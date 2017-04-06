# Create version.cpp file with appropriate version
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/version.cpp.in
                ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

add_executable(cnvetti main.cpp ${version_file})

install(
    TARGETS cnvetti
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
