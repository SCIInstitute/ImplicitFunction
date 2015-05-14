#  Permission is  hereby  granted, free  of charge, to any person
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

project(ImplicitFunction)

set(ImplicitFunction_LIB_NAME "rbf")

set(ImplicitFunction_SRCS
  BBox.cpp
  BBox.h
  BHNode.cpp
  BHNode.h
  BoundingBox.cpp
  BoundingBox.h
  ETSP.cpp
  ETSP.h
  FMM.cpp
  FMM.h
  LinearSolver.cpp
  LinearSolver.h
  RBF.cpp
  RBF.h
  SampleData.cpp
  SampleData.h
  ScatteredData.cpp
  ScatteredData.h
  SparseMatrix.cpp
  SparseMatrix.h
  fault.cpp
  fault.h
  fileIO.h
  horizon.cpp
  horizon.h
  structure.cpp
  structure.h
  tear.cpp
  tear.h
  tree.cpp
  tree.h
  vec3.cpp
  vec3.h
  RBFInterface.cpp
  RBFInterface.h
)

add_library(${ImplicitFunction_LIB_NAME}
  ${ImplicitFunction_SRCS}
)
