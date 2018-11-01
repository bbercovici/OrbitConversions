# MIT License
# Copyright (c) 2018 Benjamin Bercovici

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   May, 2018
# @brief  CMake listing enabling compilation and installation of the OrbitConversions library

if (EXISTS /home/bebe0705/.am_fortuna)
	set(IS_FORTUNA ON)
	message("-- This is Fortuna")
	set(RBK_LOC "/home/bebe0705/libs/local/lib/cmake/RigidBodyKinematics")
endif()
cmake_minimum_required(VERSION 2.8)

set(LIB_NAME "OrbitConversions" CACHE STRING "Name of library to be created.")

project(${LIB_NAME})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

add_definitions(-Wall -O2 )
set(CMAKE_CXX_FLAGS "-std=c++14")

# Find Armadillo 
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Find RBK 
find_package(RigidBodyKinematics REQUIRED PATHS ${RBK_LOC})
include_directories(${RBK_INCLUDE_DIR} )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(${LIB_NAME} SHARED 
	source/State.cpp
	source/CartState.cpp
	source/KepState.cpp)

target_link_libraries(${LIB_NAME} ${ARMADILLO_LIBRARIES} ${RBK_LIBRARY})

file(GLOB includes include/OrbitConversions/*.hpp)

# Installing
if(${IS_FORTUNA})

	install (TARGETS ${LIB_NAME} DESTINATION /home/bebe0705/libs/local/lib)
	install(FILES ${includes} DESTINATION /home/bebe0705/libs/local/include/OrbitConversions)
	install (FILES OrbitConversions/CMake/OrbitConversionsConfig.cmake DESTINATION /home/bebe0705/libs/local/lib/cmake/OrbitConversions/)

else()
	
	install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib)
	install(FILES ${includes} DESTINATION /usr/local/include/OrbitConversions)
	install (FILES OrbitConversions/CMake/OrbitConversionsConfig.cmake DESTINATION /usr/local/share/OrbitConversions/)

endif()



