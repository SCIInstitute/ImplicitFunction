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
#include "ScatteredData.h"
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
    prismWithInternalPoint2.push_back( vec3(    0, 1.0, 10.0) );
    prismWithInternalPoint2.push_back( vec3(  1.0,   0, 10.0) );
    prismWithInternalPoint2.push_back( vec3( -1.0,   0, 10.0) );
    prismWithInternalPoint2.push_back( vec3(    0, 0.5, 5.0) );

    gridSize50[0] = 50.0; gridSize50[1] = 50.0; gridSize50[2] = 50.0;
    gridSpacing1[0] = 1.0; gridSpacing1[1] = 1.0; gridSpacing1[2] = 1.0;
    normalOffset1 = 1.0;
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
  double normalOffset1;
  std::vector<axis_t> axisDataZ; // #axis entries min 3
  Kernel thinPlateKernel;
  Kernel gaussianPlateKernel;
  Kernel multiQuadKernel;
};

TEST_F(ConvexHull3DTest, PrismWithInternalPoint2BasicInterfaceTest)
{
  RBFInterface rbfInterface( prismWithInternalPoint2, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset1, axisDataZ,
                             true, false, ThinPlate );

  const ScatteredData *data = rbfInterface.getSurfaceData();
  ASSERT_TRUE( data != nullptr );

  // check scattered data setup, normal calculation for stability
  std::vector<double> x = { -1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, -0.2265412489744628, -1.7734587510255371, 0.2058958485078267, 1.7941041514921734, 0.0, 0.0, -0.4384836277759048, -1.5615163722240952, 0.6029479320322761, 1.3970520679677239, 0.0, 0.0 };
  std::vector<double> y = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.3203771272650137, -0.3203771272650137, -0.2325877093742827, 0.2325877093742827, 0.1835037736389938, 1.8164962263610063, 0.2325877133161922, -0.2325877133161922, 0.7259807942675245, -0.7259807942675245, 0.0571912277167201, 1.9428087722832799};
  std::vector<double> z = { 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 5.0, 9.4530823069676728, 10.5469176930323272, 9.4384838105643816, 10.5615161894356184, 9.4226502293769805, 10.5773497706230195, 0.7941037281936285, -0.7941037281936285, 0.5615161784353572, -0.5615161784353572, 0.3333330950255075, -0.3333330950255075 };
  std::vector<double> f = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0 };

  ASSERT_EQ( x.size(), data->surfacePoints_[0].size() );
  ASSERT_EQ( y.size(), data->surfacePoints_[1].size() );
  ASSERT_EQ( z.size(), data->surfacePoints_[2].size() );
  ASSERT_EQ( f.size(), data->fnc_.size() );

  for (size_t i = 0; i < x.size(); ++i )
  {
    ASSERT_NEAR( x[i], data->surfacePoints_[0][i], 1.0e-10 );
  }

  for (size_t i = 0; i < y.size(); ++i )
  {
    ASSERT_NEAR( y[i], data->surfacePoints_[1][i], 1.0e-10 );
  }

  for (size_t i = 0; i < z.size(); ++i )
  {
    ASSERT_NEAR( z[i], data->surfacePoints_[2][i], 1.0e-10 );
  }

  for (size_t i = 0; i < f.size(); ++i )
  {
    ASSERT_NEAR( f[i], data->fnc_[i], 1.0e-10 );
  }

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