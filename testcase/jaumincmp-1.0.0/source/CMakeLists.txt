cmake_minimum_required(VERSION 3.3)
project(JAUMIN_example CXX C Fortran)

if(NOT JAUMIN_ROOT)
  set(JAUMIN_ROOT "/home/user/hipro-devel/jaumin")
endif()
if (NOT JAUMIN_NDIM)
  set(JAUMIN_NDIM 2)
endif()
if (NOT JAUMIN_CREL)
  set(JAUMIN_CREL dbg)
endif()

list(APPEND CMAKE_MODULE_PATH "${JAUMIN_ROOT}/share/cmake/Modules/scns")
find_package(JAUMIN REQUIRED)
include_directories(${JAUMIN_INCLUDE_DIRS})
add_definitions(${JAUMIN_DEFINITIONS})

file(GLOB_RECURSE src_files ./*.C)
add_executable(jaumincmp${JAUMIN_NDIM}d ${src_files})
target_link_libraries(jaumincmp${JAUMIN_NDIM}d ${JAUMIN_LIBRARIES})
#set_target_properties(jaumincmp${JAUMIN_NDIM}d PROPERTIES
#  LINK_FLAGS "-Wl,--as-needed")
