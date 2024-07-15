cmake_minimum_required(VERSION 3.2)
project(orfanage CXX)

include(ExternalProject)
include(FetchContent)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DOPENMP_AVAILABLE)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DNOCURL=1")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")
    add_compile_options(-fsanitize=address)
    add_compile_options(-fsanitize=leak)
    add_compile_options(-fsanitize=undefined)
    add_link_options(-fsanitize=address)
    add_link_options(-fsanitize=leak)
    add_link_options(-fsanitize=undefined)
endif()

include_directories(src)

option(ORFANAGE_NATIVE_BUILD "Architecture-specific optimizations, i.e., g++ -march=native." OFF)
if (ORFANAGE_NATIVE_BUILD)
    add_definitions(-DORFANAGE_NATIVE_BUILD=1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif()
endif()

option(ORFANAGE_STATIC_BUILD "Static build." OFF)
if (ORFANAGE_STATIC_BUILD)
    add_definitions(-DORFANAGE_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        message(WARNING "WARNING: Builds on Mac are never fully static.")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
    endif()
endif()

# GCLIB
FetchContent_Declare(
  gclib
  GIT_REPOSITORY https://github.com/gpertea/gclib.git
  GIT_TAG dcb0141c1faa589e87755dbbc7f658d1ac58b3b4
)

FetchContent_MakeAvailable(gclib)

include_directories("${gclib_SOURCE_DIR}/")

add_executable(orfanage
    ${gclib_SOURCE_DIR}/gff.h
    ${gclib_SOURCE_DIR}/gff.cpp
    ${gclib_SOURCE_DIR}/GStr.h
    ${gclib_SOURCE_DIR}/GStr.cpp
    ${gclib_SOURCE_DIR}/GBase.h
    ${gclib_SOURCE_DIR}/GBase.cpp
    ${gclib_SOURCE_DIR}/codons.h
    ${gclib_SOURCE_DIR}/codons.cpp
    ${gclib_SOURCE_DIR}/gdna.h
    ${gclib_SOURCE_DIR}/gdna.cpp
    ${gclib_SOURCE_DIR}/GStr.h
    ${gclib_SOURCE_DIR}/GStr.cpp
    ${gclib_SOURCE_DIR}/GVec.hh
    ${gclib_SOURCE_DIR}/GList.hh
    ${gclib_SOURCE_DIR}/GHash.hh
    ${gclib_SOURCE_DIR}/GFaSeqGet.h
    ${gclib_SOURCE_DIR}/GFaSeqGet.cpp
    ${gclib_SOURCE_DIR}/GFastaIndex.h
    ${gclib_SOURCE_DIR}/GFastaIndex.cpp
    ${gclib_SOURCE_DIR}/GFastaFile.h
    ${CMAKE_SOURCE_DIR}/src/arg_parse.hpp
    ${CMAKE_SOURCE_DIR}/src/Transcriptome.cpp
    ${CMAKE_SOURCE_DIR}/src/Transcriptome.h
    ${CMAKE_SOURCE_DIR}/src/common.hpp
    ${CMAKE_SOURCE_DIR}/src/gff_utils.h
    ${CMAKE_SOURCE_DIR}/src/gff_utils.cpp
    ${CMAKE_SOURCE_DIR}/main.cpp)

target_link_libraries(orfanage z)

install(TARGETS orfanage DESTINATION bin)
