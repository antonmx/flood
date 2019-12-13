cmake_minimum_required (VERSION 3.1)

project(flood)

set(VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${flood_SOURCE_DIR})

add_custom_target(dist
    COMMAND git archive --prefix=${PROJECT_NAME}.${VERSION}/ HEAD
        | bzip2 > ${CMAKE_BINARY_DIR}/${PROJECT_NAME}.${VERSION}.tar.bz2
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "$ENV{CPPFLAGS} ${CMAKE_CXX_FLAGS} -ffast-math")

include(CheckIncludeFile)
include(CheckLibraryExists)

add_subdirectory(blitz-long)
include_directories("blitz-long")

# check ImageMagick
FIND_PACKAGE(ImageMagick REQUIRED)
FIND_PACKAGE(ImageMagick COMPONENTS Magick++ REQUIRED)
FIND_PACKAGE(ImageMagick COMPONENTS MagickCore REQUIRED)
include_directories("${ImageMagick_INCLUDE_DIRS}")



FIND_PACKAGE(TIFF REQUIRED)

FIND_PACKAGE(Readline REQUIRED)

find_library(HAVE_PTHREAD pthread)
if( "${HAVE_PTHREAD}" MATCHES "HAVE_PTHREAD-NOTFOUND" )
  message(FATAL_ERROR "Pthread library not found.")
endif()

check_include_file(pthread.h HAVE_PTHREAD_H ${CMAKE_CXX_FLAGS})
if( "${HAVE_PTHREAD_H}" MATCHES "HAVE_PTHREAD_H-NOTFOUND" )
  message(FATAL_ERROR "pthread.h header not found.")
endif()



find_library(HAVE_POPTMX poptmx)
if( "${HAVE_POPTMX}" MATCHES "HAVE_POPTMX-NOTFOUND" )
  message(FATAL_ERROR "Poptmx library not found: https://github.com/antonmx/poptmx")
endif()

check_include_file(poptmx.h HAVE_POPTMX_H ${CMAKE_CXX_FLAGS})
if( "${HAVE_POPTMX_H}" MATCHES "HAVE_POPTMX_H-NOTFOUND" )
  message(FATAL_ERROR "Poptmx header not found. https://github.com/antonmx/poptmx")
endif()




add_library(flood SHARED tools.h tools.cpp algorithm3d.h algorithm3d.cpp)
target_link_libraries(flood
  ${ImageMagick_LIBRARIES}
  tiff
  pthread
)
add_dependencies(flood blitz-long)
install(TARGETS flood LIBRARY DESTINATION lib)


add_executable(flood3d flood3d.cpp)
#add_executable(flood3d flood3d.cpp tools.h tools.cpp algorithm3d.h algorithm3d.cpp)
target_link_libraries(flood3d
  flood
  #${ImageMagick_LIBRARIES}
  #tiff
  #pthread
  poptmx
  readline
)
add_dependencies(flood blitz-long)
install(TARGETS flood3d RUNTIME DESTINATION bin)


#add_executable(flood2d flood2d.cpp tools.cpp)
#target_link_libraries(flood2d
#  ${ImageMagick_LIBRARIES}
#  tiff
#  poptmx
#  pthread
#)
#add_dependencies(flood2d blitz-long)
#install(TARGETS flood2d RUNTIME DESTINATION bin)


