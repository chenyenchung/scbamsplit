cmake_minimum_required(VERSION 3.18)
project(scbamsplit C)
include(CheckIPOSupported)

set(CMAKE_C_STANDARD 99)
set(CMAKE_BUILD_TYPE Debug)
#########
#https://github.com/luntergroup/octopus/blob/develop/build/cmake/modules/FindHTSlib.cmake
#########
# - Try to find htslib
# Once done, this will define
#
#  htslib_FOUND - system has htslib
#  htslib_INCLUDE_DIRS - the htslib include directories
#  htslib_LIBRARIES - link these to use htslib
#
# This code was modified from https://github.com/genome/build-common/blob/master/cmake/FindHTSlib.cmake
#

# A simple wrapper to make pkg-config searches a bit easier.
# Works the same as CMake's internal pkg_check_modules but is always quiet.
macro (libfind_pkg_check_modules)
    find_package(PkgConfig QUIET)
    if (PKG_CONFIG_FOUND)
        pkg_check_modules(${ARGN} QUIET)
    endif()
endmacro()

macro (libfind_package PREFIX)
    set (LIBFIND_PACKAGE_ARGS ${ARGN})
    if (${PREFIX}_FIND_QUIETLY)
        set (LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} QUIET)
    endif (${PREFIX}_FIND_QUIETLY)
    if (${PREFIX}_FIND_REQUIRED)
        set (LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} REQUIRED)
    endif (${PREFIX}_FIND_REQUIRED)
    find_package(${LIBFIND_PACKAGE_ARGS})
endmacro (libfind_package)

macro (libfind_process PREFIX)
    # Skip processing if already processed during this run
    if (NOT ${PREFIX}_FOUND)
        # Start with the assumption that the library was found
        set (${PREFIX}_FOUND TRUE)

        # Process all includes and set _FOUND to false if any are missing
        foreach (i ${${PREFIX}_PROCESS_INCLUDES})
            if (${i})
                set (${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIRS} ${${i}})
                mark_as_advanced(${i})
            else (${i})
                set (${PREFIX}_FOUND FALSE)
            endif (${i})
        endforeach (i)

        # Process all libraries and set _FOUND to false if any are missing
        foreach (i ${${PREFIX}_PROCESS_LIBS})
            if (${i})
                set (${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARIES} ${${i}})
                mark_as_advanced(${i})
            else (${i})
                set (${PREFIX}_FOUND FALSE)
            endif (${i})
        endforeach (i)

        # Print message and/or exit on fatal error
        if (${PREFIX}_FOUND)
            if (NOT ${PREFIX}_FIND_QUIETLY)
                message (STATUS "Found ${PREFIX} ${${PREFIX}_VERSION}")
            endif (NOT ${PREFIX}_FIND_QUIETLY)
        else (${PREFIX}_FOUND)
            if (${PREFIX}_FIND_REQUIRED)
                foreach (i ${${PREFIX}_PROCESS_INCLUDES} ${${PREFIX}_PROCESS_LIBS})
                    message("${i}=${${i}}")
                endforeach (i)
                message (FATAL_ERROR "Required library ${PREFIX} NOT FOUND.\nInstall the library (dev version) and try again. If the library is already installed, use ccmake to set the missing variables manually.")
            endif (${PREFIX}_FIND_REQUIRED)
        endif (${PREFIX}_FOUND)
    endif (NOT ${PREFIX}_FOUND)
endmacro (libfind_process)

set(HTSLIB_SEARCH_DIRS
        ${HTSLIB_SEARCH_DIRS}
        $ENV{HTSLIB_ROOT}
        ${HTSLIB_ROOT}
        )

if(NOT HTSlib_NO_SYSTEM_PATHS)
    set(HTSLIB_SEARCH_DIRS
            ${HTSLIB_SEARCH_DIRS}
            /usr
            /usr/local
            /opt/homebrew
            )
endif()

set(_htslib_ver_path "htslib-${htslib_FIND_VERSION}")

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(HTSLIB_PKGCONF htslib)

# Include dir
find_path(HTSlib_INCLUDE_DIR
        NAMES ${HTSLIB_ADDITIONAL_HEADERS} htslib/sam.h
        PATHS ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_INCLUDE_DIRS}
        PATH_SUFFIXES include htslib/${_htslib_ver_path}
        NO_DEFAULT_PATH
        )



if (HTSlib_USE_STATIC_LIBS)
    # Dependencies
    set(ZLIB_ROOT ${HTSLIB_ROOT})
    libfind_package(HTSlib ZLIB)
    set(BZip2_ROOT ${HTSLIB_ROOT})
    libfind_package(HTSlib BZip2)
    set(LibLZMA_ROOT ${HTSLIB_ROOT})
    libfind_package(HTSlib LibLZMA)
    set(CURL_ROOT ${HTSLIB_ROOT})
    libfind_package(HTSlib CURL)
    if (NOT APPLE)
        set(OpenSSL_ROOT ${HTSLIB_ROOT})
        libfind_package(HTSlib OpenSSL)
    endif()
    set(HTSlib_LIBRARY_names libhts.a)
else()
    set(HTSlib_LIBRARY_names libhts.so.3 libhts.so.2 libhts.so libhts.3.dylib libhts.2.dylib libhts.dylib)
endif()

# rpath dir
find_path(HTSlib_LIB_DIR
        NAMES ${HTSlib_LIBRARY_names}
        PATHS ${HTSlib_INCLUDE_DIR} ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_LIBRARY_DIRS}
        NO_DEFAULT_PATH
        PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu lib/aarch64-linux-gnu ${_htslib_ver_path}
)

# Finally the library itself
find_library(HTSlib_LIBRARY
        NAMES ${HTSlib_LIBRARY_names}
        PATHS ${HTSlib_INCLUDE_DIR} ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_LIBRARY_DIRS}
        NO_DEFAULT_PATH
        PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu lib/aarch64-linux-gnu ${_htslib_ver_path}
        )

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(HTSlib_PROCESS_INCLUDES HTSlib_INCLUDE_DIR)
set(HTSlib_PROCESS_LIBS HTSlib_LIBRARY)

if (HTSlib_USE_STATIC_LIBS)
    set(HTSlib_PROCESS_INCLUDES ${HTSlib_PROCESS_INCLUDES}
            ZLIB_INCLUDE_DIR
            BZIP2_INCLUDE_DIR
            LIBLZMA_INCLUDE_DIRS
            CURL_INCLUDE_DIRS)
    set(HTSlib_PROCESS_LIBS ${HTSlib_PROCESS_LIBS}
            ZLIB_LIBRARIES
            BZIP2_LIBRARIES
            LIBLZMA_LIBRARIES
            CURL_LIBRARIES)
    if (NOT APPLE)
        set(HTSlib_PROCESS_INCLUDES ${HTSlib_PROCESS_INCLUDES} OPENSSL_INCLUDE_DIR)
        set(HTSlib_PROCESS_LIBS ${HTSlib_PROCESS_LIBS} OPENSSL_LIBRARIES)
    endif()
endif()

libfind_process(HTSlib)

list(APPEND CMAKE_BUILD_RPATH ${HTSlib_LIB_DIR})
set(CMAKE_INSTALL_RPATH "${CMAKE_BUILD_RPATH}")
find_package(Threads REQUIRED)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)

## Check if compiler supports inter-procedural optimization
## https://stackoverflow.com/questions/31355692/how-do-i-enable-link-time-optimization-lto-with-cmake
check_ipo_supported(RESULT IPO_SUPPORT OUTPUT IPO_ERROR)

message(STATUS "   Building with ${CMAKE_BUILD_TYPE} configuration...")
message(STATUS "   HTSlib include dirs: ${HTSlib_INCLUDE_DIRS}")
message(STATUS "   HTSlib libraries: ${HTSlib_LIBRARIES}")
message(STATUS "   HTSlib link path: ${HTSlib_LIB_DIR}")
message(STATUS "   CMake Build Rpath: ${CMAKE_BUILD_RPATH}")
message(STATUS "   CMake Install Rpath: ${CMAKE_INSTALL_RPATH}")

#########
# add extra include directories

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    message(STATUS "   With sanitizer")
#    add_compile_options(-fsanitize=thread)
#    add_link_options(-fsanitize=thread)
#    add_compile_options(-fsanitize=address)
#    add_link_options(-fsanitize=address)
endif()
add_executable(${PROJECT_NAME}
        src/main.c
        src/hash.c
        src/utils.c
        src/sort.c
        src/thread_pool.c
        src/thread_pool.h)
if ( IPO_SUPPORT )
    if (NOT CMAKE_BUILD_TYPE MATCHES "Debug")
        message(STATUS "Enabling link-time optimization")
        # Using link-time optimization prevents setting breakpoints for debugging
        set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
else()
    message(STATUS "IPO not supported")
endif()
target_include_directories(${PROJECT_NAME} PUBLIC ${HTSlib_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSlib_LIBRARIES} Threads::Threads)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
