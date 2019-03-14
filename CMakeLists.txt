#--------------------------------------
# Project settings
#--------------------------------------

cmake_minimum_required(VERSION 3.1)

project(HoaLibrary CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_SUPPRESS_REGENERATION true)

#--------------------------------------
# Options
#--------------------------------------

option(HOA_BUILD_TESTS "Build HoaLibrary tests" ON)
option(GCOV_SUPPORT "Build for gcov" OFF)

#--------------------------------------
# Options
#--------------------------------------

set(HOALIBRARY_TARGET_NAME ${PROJECT_NAME})
set(HOALIBRARY_INCLUDE_BUILD_DIR "${PROJECT_SOURCE_DIR}/Sources/")

#--------------------------------------
# Compiler settings
#--------------------------------------

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.9)
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")

if(UNIX)
  if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
  endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

#--------------------------------------
# Create Library
#--------------------------------------

add_library(${HOALIBRARY_TARGET_NAME} INTERFACE)
add_library(
  ${PROJECT_NAME}::${HOALIBRARY_TARGET_NAME}
  ALIAS ${HOALIBRARY_TARGET_NAME}
)

target_include_directories(
  ${HOALIBRARY_TARGET_NAME}
  INTERFACE
    $<BUILD_INTERFACE:${HOALIBRARY_INCLUDE_BUILD_DIR}>
    $<INSTALL_INTERFACE:include>
)

file(GLOB HOALIBRARY_SOURCES ${HOALIBRARY_INCLUDE_BUILD_DIR}*.hpp)

target_sources(
  ${HOALIBRARY_TARGET_NAME}
  INTERFACE
    $<BUILD_INTERFACE:${HOALIBRARY_SOURCES}>
)

#--------------------------------------
# Tests
#--------------------------------------

if(${HOA_BUILD_TESTS})
  enable_testing()
  file(GLOB TESTSOURCES ${PROJECT_SOURCE_DIR}/Tests/*.cpp)
  source_group(Tests FILES ${TESTSOURCES})
  include_directories(${PROJECT_SOURCE_DIR}/Tests)

  add_executable(hoatest ${TESTSOURCES} ${HOALIBRARY_SOURCES})
  target_link_libraries(hoatest ${HOALIBRARY_TARGET_NAME})
  add_test(hoatestall hoatest)

  if(${GCOV_SUPPORT})
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
      set(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
      target_link_libraries(hoatest gcov)
  endif()
endif()
