cmake_minimum_required(VERSION 3.31)
find_package(samurai)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(lsburgers_amr lsburgers_amr.cpp)
target_link_libraries(lsburgers_amr samurai::samurai)

add_executable(lsburgers_mra lsburgers_mra.cpp)
target_link_libraries(lsburgers_mra samurai::samurai)

add_executable(MRA_compression MRA_compression.cpp)
target_link_libraries(MRA_compression samurai::samurai)
