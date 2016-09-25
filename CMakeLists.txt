cmake_minimum_required(VERSION 2.8)

project(Hoa)
option(GCOV_SUPPORT "Build for gcov" Off)

file(GLOB HOASOURCES ${PROJECT_SOURCE_DIR}/Sources/*.cpp ${PROJECT_SOURCE_DIR}/Sources/*.hpp)
source_group(Hoa FILES ${HOASOURCES})
include_directories(${PROJECT_SOURCE_DIR}/Tests)
include_directories(${PROJECT_SOURCE_DIR}/Sources)

set(CMAKE_MACOSX_RPATH On)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.7)

add_executable(hoatest ${HOASOURCES} ${PROJECT_SOURCE_DIR}/Tests/Test.cpp)
set_target_properties(hoatest PROPERTIES OUTPUT_NAME test)
set_target_properties(hoatest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

if(${GCOV_SUPPORT} STREQUAL "On")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS="-fprofile-arcs -ftest-coverage")
    target_link_libraries(hoatest gcov)
endif()
