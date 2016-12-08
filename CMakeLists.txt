# Permission is  hereby  granted, free  of charge, to any person
# obtaining a copy of this software and associated documentation
# files  ( the "Software" ),  to  deal in  the  Software without
# restriction, including  without limitation the rights to  use,
# copy, modify,  merge, publish, distribute, sublicense,  and/or
# sell copies of the Software, and to permit persons to whom the
# Software is  furnished  to do  so,  subject  to  the following
# conditions:
#
# The above  copyright notice  and  this permission notice shall
# be included  in  all copies  or  substantial  portions  of the
# Software.
#
# THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
# KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
# WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.

cmake_minimum_required(VERSION 2.8)

project(ImplicitFunction)

include(CTest)
option(BUILD_TESTING "Build with tests." ON)
if(BUILD_TESTING)
  find_package(GTest 1.7 REQUIRED)
endif()

set(EIGEN3_INSTALL_DIR "" CACHE PATH "Eigen 3 include directory.")
if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
  find_path(EIGEN3_INCLUDE_DIR signature_of_eigen3_matrix_library
    PATHS ${EIGEN3_INSTALL_DIR} ${EIGEN3_INSTALL_DIR}/include
    PATH_SUFFIXES eigen3 eigen
    NO_DEFAULT_PATH
  )
  if (NOT EIGEN3_INCLUDE_DIR)
    MESSAGE(ERROR "${EIGEN3_INCLUDE_DIR} not valid Eigen include directory.")
  endif()
endif()

set(Tetgen_DIR "" CACHE PATH "Tetgen directory.")
option(TravisCI_BUILD "Enable if Travis-CI build." OFF)

# have to hardcode some setting for the Travis-CI build
# Tetgen cmake files downloaded by Travis-CI build
# TODO: is there a better way to handle this?
if(TravisCI_BUILD)
  set(TETGEN_INCLUDE ${Tetgen_DIR})
  set(TETGEN_LIBRARY_DIR ${Tetgen_DIR})
  set(TETGEN_USE_FILE ${Tetgen_DIR}/UseTetgen.cmake)
  set(TETGEN_LIBRARY "tet")
  configure_file(${Tetgen_DIR}/TetgenConfig.cmake.in ${Tetgen_DIR}/TetgenConfig.cmake @ONLY)
  configure_file(${Tetgen_DIR}/UseTetgen.cmake ${TETGEN_USE_FILE} COPYONLY)
endif()

find_package(Tetgen CONFIGS TetgenConfig.cmake HINTS ${Tetgen_DIR} NO_SYSTEM_ENVIRONMENT_PATH)
if(NOT Tetgen_FOUND)
  message(FATAL_ERROR "Tetgen library not found in ${Tetgen_DIR}")
endif()
include(${TETGEN_USE_FILE})

###########################################
# *Nix C++ compiler flags
###########################################

if(UNIX)
  option(ENABLE_COVERAGE "Enable coverage compiler options." OFF)
  if(ENABLE_COVERAGE)
    #if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    #  set(COVERAGE_FLAGS "")
    #else()
    #  # probably GCC
    #  set(COVERAGE_FLAGS " -fprofile-arcs -ftest-coverage")
    #endif()
    set(COVERAGE_FLAGS "--coverage")
    set(CTEST_CUSTOM_COVERAGE_EXCLUDE "/Tests/")
  else()
    set(COVERAGE_FLAGS "")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall ${COVERAGE_FLAGS}")
  if(APPLE)
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -ftemplate-depth=256")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wshorten-64-to-32 ${CMAKE_CXX_FLAGS_DEBUG}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -ldl -lrt")
  endif()
endif()

set(ImplicitFunction_LIB_NAME "rbf")

set(ImplicitFunction_SRCS
  BBox.cpp
  BBox.h
  BHNode.cpp
  BHNode.h
  ETSP.cpp
  ETSP.h
  FMM.cpp
  FMM.h
  LinearSolver.cpp
  LinearSolver.h
  RBF.cpp
  RBF.h
  ScatteredData.cpp
  ScatteredData.h
  SparseMatrix.cpp
  SparseMatrix.h
  vec3.cpp
  vec3.h
  RBFInterface.cpp
  RBFInterface.h
)

include_directories(${EIGEN3_INCLUDE_DIR})

add_library(${ImplicitFunction_LIB_NAME}
  ${ImplicitFunction_SRCS}
)

target_link_libraries(${ImplicitFunction_LIB_NAME}
  ${TETGEN_LIBRARY}
)

if(BUILD_TESTING)
  add_subdirectory(Tests)
endif()

set(ImplicitFunction_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR})
set(ImplicitFunction_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(ImplicitFunction_USE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/UseImplicitFunction.cmake")

configure_file(ImplicitFunctionConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ImplicitFunctionConfig.cmake @ONLY)
