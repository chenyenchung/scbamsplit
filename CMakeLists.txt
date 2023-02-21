cmake_minimum_required(VERSION 3.22)
project(scbamsplit C)

set(CMAKE_C_STANDARD 99)

#include_directories(/usr/local/Cellar/htslib/1.16/include)
#link_directories(/usr/local/Cellar/htslib/1.16/lib/)
find_library(htslib SHARED IMPORTED) # or STATIC instead of SHARED
#set_target_properties(htslib PROPERTIES
#        IMPORTED_LOCATION "/usr/local/Cellar/htslib/1.16/lib/libhts.1.16.dylib"
#        INTERFACE_INCLUDE_DIRECTORIES "/usr/local/Cellar/htslib/1.16/include"
#        )

add_executable(scbamsplit main.c)
# target_link_libraries(scbamsplit htslib)