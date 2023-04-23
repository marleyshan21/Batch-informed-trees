cmake_minimum_required(VERSION 3.0)

project(BitStar)

set(CMAKE_CXX_STANDARD 11)
include_directories(include)


find_package(OpenCV 4 REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})



add_executable(BitStar unit_test2.cpp )

target_link_libraries(BitStar ${OpenCV_LIBS})