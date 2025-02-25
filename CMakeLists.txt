# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( project2 )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5 )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED COMPONENTS system filesystem )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()



# Creating entries for target: project
# Link the executable to CGAL and third-party libraries
# ############################


add_executable(opt_triangulation src/main.cpp src/implementation.cpp src/simAnnealing.cpp src/localSearch.cpp src/antColony.cpp)
target_link_libraries(opt_triangulation PRIVATE CGAL::CGAL )


if (CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER)
  target_link_libraries(opt_triangulation PRIVATE CGAL::CGAL_Qt5)
endif()

#add_to_cached_list( CGAL_EXECUTABLE_TARGETS tringulation )
