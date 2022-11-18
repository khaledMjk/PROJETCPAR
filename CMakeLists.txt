cmake_minimum_required(VERSION 3.22)

project(parallele)

find_package(MPI REQUIRED)

message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

add_executable(sandbox.exe sandbox.cpp)
target_link_libraries(sandbox.exe PUBLIC MPI::MPI_CXX)

add_library(runge_kutta runge_kutta.cpp)
target_link_libraries(sandbox.exe PUBLIC runge_kutta)

add_library(cartesian_grid_of_speed cartesian_grid_of_speed.cpp)
target_link_libraries(sandbox.exe PUBLIC cartesian_grid_of_speed)

add_library(cloud_of_points cloud_of_points.cpp)
target_link_libraries(sandbox.exe PUBLIC cloud_of_points)

add_library(vortex vortex.cpp)
target_link_libraries(sandbox.exe PUBLIC vortex)

target_link_libraries(sandbox.exe PUBLIC MPI::MPI_CXX)
