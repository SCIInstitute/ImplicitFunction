//-------------------------------------------------------------------
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include "RBFInterface.h"
#include "vec3.h"

class ConvexHull3DTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    // TODO: new coordinates needed to match Seg3D
//    tet.push_back( vec3( 1.0,  1.0,  1.0) );
//    tet.push_back( vec3(-1.0, -1.0,  1.0) );
//    tet.push_back( vec3(-1.0,  1.0, -1.0) );
//    tet.push_back( vec3( 1.0, -1.0, -1.0) );
//
//    cube.push_back( vec3(-1.0, -1.0,  1.0) );
//    cube.push_back( vec3( 1.0, -1.0,  1.0) );
//    cube.push_back( vec3( 1.0,  1.0,  1.0) );
//    cube.push_back( vec3(-1.0,  1.0,  1.0) );
//    cube.push_back( vec3(-1.0, -1.0, -1.0) );
//    cube.push_back( vec3( 1.0, -1.0, -1.0) );
//    cube.push_back( vec3( 1.0,  1.0, -1.0) );
//    cube.push_back( vec3(-1.0,  1.0, -1.0) );
//
//    prism.push_back( vec3(-1.0, -0.58, -1.0) );
//    prism.push_back( vec3(-1.0, -0.58,  1.0) );
//    prism.push_back( vec3( 1.0, -0.58, -1.0) );
//    prism.push_back( vec3( 1.0, -0.58,  1.0) );
//    prism.push_back( vec3(   0,  1.15, -1.0) );
//    prism.push_back( vec3(   0,  1.15,  1.0) );
//
//    prismWithInternalPoint.push_back( vec3(-1.0, -0.58, -1.0) );
//    prismWithInternalPoint.push_back( vec3(-1.0, -0.58,  1.0) );
//    prismWithInternalPoint.push_back( vec3( 1.0, -0.58, -1.0) );
//    prismWithInternalPoint.push_back( vec3( 1.0, -0.58,  1.0) );
//    prismWithInternalPoint.push_back( vec3(   0,  1.15, -1.0) );
//    prismWithInternalPoint.push_back( vec3(   0,  1.15,  1.0) );
//    prismWithInternalPoint.push_back( vec3(   0,     0,    0) );
//
//    prism2.push_back( vec3(    0, 1.0,   0) );
//    prism2.push_back( vec3(  1.0,   0,   0) );
//    prism2.push_back( vec3( -1.0,   0,   0) );
//    prism2.push_back( vec3(    0, 1.0, 1.0) );
//    prism2.push_back( vec3(  1.0,   0, 1.0) );
//    prism2.push_back( vec3( -1.0,   0, 1.0) );

    prismWithInternalPoint2.push_back( vec3(    0, 1.0,   0) );
    prismWithInternalPoint2.push_back( vec3(  1.0,   0,   0) );
    prismWithInternalPoint2.push_back( vec3( -1.0,   0,   0) );
    prismWithInternalPoint2.push_back( vec3(    0, 1.0, 1.0) );
    prismWithInternalPoint2.push_back( vec3(  1.0,   0, 1.0) );
    prismWithInternalPoint2.push_back( vec3( -1.0,   0, 1.0) );
    prismWithInternalPoint2.push_back( vec3(    0, 0.5, 0.5) );

    gridSize50[0] = 50.0; gridSize50[1] = 50.0; gridSize50[2] = 50.0;
    gridSpacing1[0] = 1.0; gridSpacing1[1] = 1.0; gridSpacing1[2] = 1.0;
    normalOffset = 2.0;
    axisDataZ.insert(axisDataZ.begin(), 3, axis_t::Z);
    thinPlateKernel = ThinPlate;
    gaussianPlateKernel = Gaussian;
    multiQuadKernel = MultiQuadratic;
  }

// virtual void TearDown() {}

//  std::vector<vec3> tet;
//  std::vector<vec3> cube;
//  std::vector<vec3> prism;
//  std::vector<vec3> prismWithInternalPoint;
//  std::vector<vec3> prism2;
  std::vector<vec3> prismWithInternalPoint2;

  vec3 origin0; // (0, 0, 0)
  vec3 gridSize50; // (50, 50, 50)
  vec3 gridSpacing1; // (1, 1, 1)
  double normalOffset;
  std::vector<axis_t> axisDataZ; // #axis entries min 3
  Kernel thinPlateKernel;
  Kernel gaussianPlateKernel;
  Kernel multiQuadKernel;
};

TEST_F(ConvexHull3DTest, PrismWithInternalPoint2BasicInterfaceTest)
{
  RBFInterface rbfInterface( prismWithInternalPoint2, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset, axisDataZ,
                             true, false, ThinPlate );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
}

//TEST_F(ConvexHull3DTest, Prism2BasicInterfaceTest)
//{
//  RBFInterface rbfInterface( prism2, origin0,
//                             gridSize50, gridSpacing1,
//                             normalOffset10, axisDataZ,
//                             true, false, ThinPlate );
//  double threshold = rbfInterface.getThresholdValue();
//  ASSERT_EQ( threshold, 0 ); // default
//}
//
//TEST_F(ConvexHull3DTest, PrismWithInternalPointBasicInterfaceTest)
//{
//  RBFInterface rbfInterface( prismWithInternalPoint, origin0,
//                             gridSize50, gridSpacing1,
//                             normalOffset10, axisDataZ,
//                             true, false, ThinPlate );
//  double threshold = rbfInterface.getThresholdValue();
//  ASSERT_EQ( threshold, 0 ); // default
//}
//
//TEST_F(ConvexHull3DTest, PrismBasicInterfaceTest)
//{
//  RBFInterface rbfInterface( prism, origin0,
//                             gridSize50, gridSpacing1,
//                             normalOffset10, axisDataZ,
//                             true, false, ThinPlate );
//  double threshold = rbfInterface.getThresholdValue();
//  ASSERT_EQ( threshold, 0 ); // default
//}
//
//TEST_F(ConvexHull3DTest, CubeBasicInterfaceTest)
//{
//  RBFInterface rbfInterface( cube, origin0,
//                             gridSize50, gridSpacing1,
//                             normalOffset10, axisDataZ,
//                             true, false, ThinPlate );
//  double threshold = rbfInterface.getThresholdValue();
//  ASSERT_EQ( threshold, 0 ); // default
//}
//
//TEST_F(ConvexHull3DTest, TetBasicInterfaceTest)
//{
//  RBFInterface rbfInterface( tet, origin0,
//                             gridSize50, gridSpacing1,
//                             normalOffset10, axisDataZ,
//                             true, false, ThinPlate );
//  double threshold = rbfInterface.getThresholdValue();
//  ASSERT_EQ( threshold, 0 ); // default
//}