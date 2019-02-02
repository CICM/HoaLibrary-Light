cmake_minimum_required(VERSION 3.1)

project(Hoa CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MACOSX_RPATH On)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.9)
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")

set(CMAKE_SUPPRESS_REGENERATION true)
enable_testing()
option(GCOV_SUPPORT "Build for gcov" Off)

if(UNIX)
  if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat -Wno-weak-vtables -Wno-padded -Wno-float-equal -Wno-double-promotion -Wno-disabled-macro-expansion")
  endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

file(GLOB HOASOURCES ${PROJECT_SOURCE_DIR}/Sources/*.hpp)
source_group(Hoa FILES ${HOASOURCES})
include_directories(${PROJECT_SOURCE_DIR}/Sources)

file(GLOB TESTSOURCES ${PROJECT_SOURCE_DIR}/Tests/*.cpp)
source_group(Tests FILES ${TESTSOURCES})
include_directories(${PROJECT_SOURCE_DIR}/Tests)

add_executable(hoatest ${HOASOURCES} ${TESTSOURCES})
add_test(hoatestall hoatest)

if(${GCOV_SUPPORT} STREQUAL "On")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
    target_link_libraries(hoatest gcov)
endif()
