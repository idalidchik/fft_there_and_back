cmake_minimum_required(VERSION 2.6)
project(fft_there_and_back C CXX)

if(UNIX)
set(CMAKE_CXX_FLAGS "-Wall -std=c++14 ")
endif(UNIX)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
message(FATAL_ERROR "In-source builds not allowed.")
endif()

# Default build is in Release mode
if(NOT CMAKE_BUILD_TYPE)
    if(MSVC)
        set(CMAKE_BUILD_TYPE "Debug")
    else()
        set(CMAKE_BUILD_TYPE "Release")
    endif(MSVC)
endif(NOT CMAKE_BUILD_TYPE)

if(WIN32)
  message(STATUS "WIN32: ${WIN32}")
  add_definitions(-DSDL_OS_WIN32 -DWIN32)
endif(WIN32)

if(UNIX)
  message(STATUS "UNIX: ${UNIX}")
  add_definitions(-DSDL_OS_UNIX)
endif(UNIX) 

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DSDL_DEBUG=1)
else()
    add_definitions(-DSDL_DEBUG=0)
endif()

#set the default path for built executables to the "bin" directory
#set the default path for built libraries to the "lib" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories(${ADDITIONAL_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

set( SDL_SOURCE_FILES
  src/main.cpp 
  src/fft_real.h
  src/fft_real.cpp
  src/fft_t.h
  )

source_group(src FILES ${SDL_SOURCE_FILES} )

add_executable(fft_there_and_back ${SDL_SOURCE_FILES} )
