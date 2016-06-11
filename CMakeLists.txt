cmake_minimum_required(VERSION 2.8)

project(Hoa)
option(COVERALLS "Build with coveralls" OFF)

set(HOASOURCES
  ${PROJECT_SOURCE_DIR}/Sources/Decoder.hpp
  ${PROJECT_SOURCE_DIR}/Sources/HrirIrc1002C3D.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Rotate.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Defs.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Math.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Scope.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Encoder.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Meter.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Signal.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Exchanger.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Optim.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Source.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Harmonics.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Planewaves.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Tools.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Hoa.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Processor.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Vector.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Hrir.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Projector.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Voronoi.hpp
  ${PROJECT_SOURCE_DIR}/Sources/HrirIrc1002C2D.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Recomposer.hpp
  ${PROJECT_SOURCE_DIR}/Sources/Wider.hpp)

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
