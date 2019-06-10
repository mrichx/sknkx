cmake_minimum_required(VERSION 2.8)

# Project name
project(sknkx_opt)

# environment variables
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:PKG_CONFIG_PATH=/home/maxr/Repos/Ipopt/build/lib64/pkgconfig:/home/maxr/Repos/Ipopt/build/lib/pkgconfig:/home/maxr/Repos/Ipopt/build/share/pkgconfig: pkg-config --cflags ipopt")

# enable C++11 Standard
set(CMAKE_CXX_STANDARD 14)

# find IPOPT package
find_package(PkgConfig REQUIRED)
pkg_search_module(IPOPT REQUIRED ipopt)
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()


# include directories
include_directories(
  include
  ${IPOPT_INCLUDE_DIRS}
)

# additional link directories
link_directories(
  ${IPOPT_LIBRARY_DIRS}
)

# message(STATUS $ENV{IPOPT_INCLUDE_DIRS})

# define library
add_library(sknkx_opt
  # src/sknkxCore.cpp
  src/sknkxNlp.cpp
  src/sknkxBlp.cpp
  src/sknkxUtil.cpp
  src/modelfunc.cpp
  src/constraintfunc.cpp
  src/defect.cpp
)

# Required on Unix OS family to be able to be linked into shared libraries.
# set_target_properties(${PROJECT_NAME}
#                       PROPERTIES POSITION_INDEPENDENT_CODE ON)

# link libraries
target_link_libraries(sknkx_opt
  ${IPOPT_LIBRARIES}
)

# Expose project's public includes
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include
                                 ${IPOPT_INCLUDE_DIRS}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

# install executable/libraries
install(TARGETS
  sknkx_opt DESTINATION ${CMAKE_BINARY_DIR}/lib
)

# install option files for
install(DIRECTORY config/ DESTINATION ${CMAKE_BINARY_DIR}/config
        FILES_MATCHING PATTERN "*.opt")
