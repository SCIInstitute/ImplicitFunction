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

set(ImplicitFunction_Tests_SRCS
  ConvexHull2DTests.cpp
  ConvexHull3DTests.cpp
  RBFTests.cpp
  RBFInterfaceTests.cpp
  ScatteredDataTests.cpp
  Vec3Tests.cpp
  Seg3DIntegrationTest.cpp
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ImplicitFunction_BINARY_DIR})

include_directories(
  ${GTEST_INCLUDE_DIR}
  ${ImplicitFunction_SOURCE_DIR}
)

add_executable(ImplicitFunction_Tests ${ImplicitFunction_Tests_SRCS})

set_property(TARGET ImplicitFunction_Tests PROPERTY CXX_STANDARD 17)

target_compile_definitions(ImplicitFunction_Tests
  PUBLIC -DREGRESSION_DIR="${CMAKE_CURRENT_SOURCE_DIR}/regression")

target_link_libraries(ImplicitFunction_Tests
  ${ImplicitFunction_LIB_NAME}
  ${TETGEN_LIBRARY}
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
)

IF(UNIX)
  TARGET_LINK_LIBRARIES(ImplicitFunction_Tests -lpthread -lm)
ENDIF()

ADD_TEST("ImplicitFunction_Tests" ${ImplicitFunction_BINARY_DIR}/ImplicitFunction_Tests)
