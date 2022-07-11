cmake_minimum_required(VERSION 3.2)
project(orfanage)

include(ExternalProject)
include(FetchContent)

find_package (OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 11)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DNOCURL=1")

add_compile_options(-fsanitize=address)
add_compile_options(-fsanitize=leak)
add_compile_options(-fsanitize=undefined)
add_link_options(-fsanitize=address)
add_link_options(-fsanitize=leak)
add_link_options(-fsanitize=undefined)

include_directories(src)

option(ORFANAGE_BUILD_LIBBIGWIG "Turn on/off building step of libBigWig (used for bioconda)." ON)
if (NOT ORFANAGE_BUILD_LIBBIGWIG)
    add_definitions(-DORFANAGE_BUILD_LIBBIGWIG=0)
    set(LIBBIGWIG_MAKE_CMD "")
else ()
    add_definitions(-DORFANAGE_BUILD_LIBBIGWIG=1)
    set(LIBBIGWIG_MAKE_CMD "make")
endif ()

option(ORFANAGE_NATIVE_BUILD "Architecture-specific optimizations, i.e., g++ -march=native." OFF)
if (ORFANAGE_NATIVE_BUILD)
    add_definitions(-DORFANAGE_NATIVE_BUILD=1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif ()
endif ()

option(ORFANAGE_STATIC_BUILD "Static build." OFF)
if (ORFANAGE_STATIC_BUILD)
    add_definitions(-DORFANAGE_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        message(WARNING "WARNING: Builds on Mac are never fully static.")
    else ()
        set (CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
    endif ()
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    # if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    #     set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    # endif ()
endif ()

ExternalProject_Add(libBigWig
        GIT_REPOSITORY "https://github.com/alevar/libBigWig"
        UPDATE_COMMAND ""
        DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/include/"
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/libBigWig"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND "${LIBBIGWIG_MAKE_CMD}"
        INSTALL_COMMAND ""
        )

include_directories(${CMAKE_SOURCE_DIR}/include/libBigWig)
link_directories(${CMAKE_SOURCE_DIR}/include/libBigWig)

# KLIB
FetchContent_Declare(klib
        GIT_REPOSITORY https://github.com/lh3/ksw2.git
        )
FetchContent_MakeAvailable(klib)
message(STATUS "klib source dir: ${klib_SOURCE_DIR}")
message(STATUS "klib binary dir: ${klib_BINARY_DIR}")
include_directories(${klib_SOURCE_DIR})

# GCLIB
include_directories("${CMAKE_SOURCE_DIR}/include/")

# HTSLIB
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoreconf --install COMMAND ./configure --host=x86_64 --without-libdeflate --disable-libcurl
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(orfanage
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.h
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GVec.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GList.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFastaIndex.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GFastaIndex.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFastaFile.h
        ${CMAKE_SOURCE_DIR}/src/arg_parse.hpp
        ${CMAKE_SOURCE_DIR}/src/Transcriptome.cpp
        ${CMAKE_SOURCE_DIR}/src/Transcriptome.h
        ${CMAKE_SOURCE_DIR}/src/common.hpp
        ${CMAKE_SOURCE_DIR}/src/gff_utils.h
        ${CMAKE_SOURCE_DIR}/src/gff_utils.cpp
        ${klib_SOURCE_DIR}/ksw2.h
        ${klib_SOURCE_DIR}/ksw2_extz2_sse.c
        ${klib_SOURCE_DIR}/ksw2_extz.c
        ${klib_SOURCE_DIR}/ksw2_gg2_sse.c
        ${klib_SOURCE_DIR}/ksw2_gg.c
        ${CMAKE_SOURCE_DIR}/src/aln.hpp
        ${CMAKE_SOURCE_DIR}/main.cpp)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(orfanage PUBLIC DEBUG=1)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_definitions(orfanage PUBLIC DEBUG=0)
endif()

add_dependencies(orfanage htslib libBigWig)

if (APPLE AND ORFANAGE_STATIC_BUILD)
    target_link_libraries(orfanage ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
else ()
    target_link_libraries(orfanage ${CMAKE_SOURCE_DIR}/include/libBigWig/libBigWig.a ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
endif ()

install(TARGETS orfanage DESTINATION ${CMAKE_INSTALL_PREFIX})