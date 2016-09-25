cmake_minimum_required(VERSION 2.8)

project(Hoa)
option(COVERALLS "Build with coveralls" OFF)

file(GLOB HOASOURCES ${PROJECT_SOURCE_DIR}/Sources/*.cpp ${PROJECT_SOURCE_DIR}/Sources/*.hpp)
source_group(Hoa FILES ${HOASOURCES})
include_directories(${PROJECT_SOURCE_DIR}/Test)
include_directories(${PROJECT_SOURCE_DIR}/Sources)

set(CMAKE_MACOSX_RPATH On)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.7)

if(${COVERALLS})
    message(STATUS "Build with coveralls")
    ENABLE_TESTING()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

add_executable(hoatest ${HOASOURCES} ${PROJECT_SOURCE_DIR}/Test/Test.cpp)
set_target_properties(hoatest PROPERTIES OUTPUT_NAME test)
set_target_properties(hoatest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)
