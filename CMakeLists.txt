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

set(QHULL_LIBRARY_DIR "" CACHE PATH "QHull library directory.")
if (IS_DIRECTORY ${QHULL_LIBRARY_DIR})
  # TODO: this probably won't work in VS or Xcode
  find_library(QHULL_STATIC_LIBRARY
    NAMES qhullstatic qhullstatic_d
    PATHS ${QHULL_LIBRARY_DIR}
    DOC "QHull static library."
  )
  find_library(QHULL_CPP_LIBRARY
    NAMES qhullcpp qhullcpp_d
    PATHS ${QHULL_LIBRARY_DIR}
    DOC "QHull C++ library."
  )
  message(STATUS "QHull libs (debug): ${QHULL_STATIC_LIBRARY} ${QHULL_CPP_LIBRARY}")
endif()

###########################################
# *Nix C++ compiler flags
###########################################

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
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

add_library(${ImplicitFunction_LIB_NAME}
  ${ImplicitFunction_SRCS}
)

set(ImplicitFunction_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR})
set(ImplicitFunction_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(ImplicitFunction_USE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/UseImplicitFunction.cmake")

configure_file(ImplicitFunctionConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ImplicitFunctionConfig.cmake @ONLY)
