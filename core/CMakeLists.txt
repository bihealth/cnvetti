# Require CMake >=3.2
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# Get version from Git tag
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake ${CMAKE_MODULE_PATH})
include(GetGitRevisionDescription)
git_describe(VERSION --tags --dirty=-dirty)

# Parse the version information into its components
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Define project and version
message(STATUS "VERSION ${VERSION}")
message(STATUS "VERSION_SHORT ${VERSION_SHORT}")
project(cnvetti_core VERSION "${VERSION_SHORT}" LANGUAGES CXX)

# Setup output directories for bin/lib
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include sub directories to get the build system bootstrapped
include_directories(include)
add_subdirectory(src)
