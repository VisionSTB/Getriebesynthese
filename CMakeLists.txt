cmake_minimum_required(VERSION 3.9)
project(synthese
  LANGUAGES CXX
  VERSION 0.1.17
  )

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Werror -pedantic-errors")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

option(BUILD_DOCS "Build the documentation" OFF)

#add_executable(${PROJECT_NAME} source/main.cpp)

#find_package(Threads)
find_package(SDL2)

include_directories(/usr/include/SDL2)
include_directories(include)

add_subdirectory(source)
#add_subdirectory(include)

set(SDL2_USE_STATIC_LIBS        ON) # only find static libs
set(SDL2_USE_MULTITHREADED      OFF)
find_package(SDL2      REQUIRED sdl2)

include(FindPkgConfig)

pkg_search_module(SDL2 REQUIRED sdl2)

include_directories(${SDL2_INCLUDE_DIRS})

if(SDL2_FOUND)
  include_directories(${SDL2_INCLUDE_DIRS})
endif()

if(${BUILD_DOCS})
  add_subdirectory(docs)
endif()

