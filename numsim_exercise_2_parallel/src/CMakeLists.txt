cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim_parallel)

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).
# Also specify the source files.
add_executable(${PROJECT_NAME}
   main.cpp
   output_writer/output_writer_paraview.cpp
   output_writer/output_writer_text.cpp
   output_writer/output_writer.cpp
   settings_parser/settings.cpp
   array2d/array2d.cpp
   discretization/central_differences.cpp
   discretization/discretization.cpp
   discretization/donor_cell.cpp
   discretization/staggered_grid.cpp
   pressure_solver/pressure_solver.cpp
#   pressure_solver/gauss_seidel.cpp
   pressure_solver/SOR_parallel.cpp
   message_passing/message_passer.cpp
)

# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

# Search for the external package "VTK"
find_package(VTK)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
  include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
  target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

# install numsim executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)

# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# for MPI
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#copy settings files to buid/src/
file(COPY lid_driven_cavity.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY driven_cavity.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#copy settings files to build/
file(COPY lid_driven_cavity.txt DESTINATION ${PROJECT_SOURCE_DIR}/../build)
file(COPY driven_cavity.txt DESTINATION ${PROJECT_SOURCE_DIR}/../build)

#gprof option
option(PROFILE "Add flags to profile the program with gprof." OFF)
if(PROFILE)
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pg -g")
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg -g")
endif()
