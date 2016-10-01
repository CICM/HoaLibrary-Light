cmake_minimum_required(VERSION 2.8)

project(Hoa)

set(CMAKE_MACOSX_RPATH On)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.7)
enable_testing()
option(GCOV_SUPPORT "Build for gcov" Off)

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
