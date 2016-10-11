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

// Note: empty point vector results in segfault
// TODO: should be an error...

class RBFInterfaceTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    pointsTriangleClockwise.push_back( vec3(10.7, 38.5, 10.0) );
    pointsTriangleClockwise.push_back( vec3(34.7, 37.8, 10.0) );
    pointsTriangleClockwise.push_back( vec3(22.3, 15.0, 10.0) );
    gridSize50[0] = 50.0; gridSize50[1] = 50.0; gridSize50[2] = 50.0;
    gridSpacing1[0] = 1.0; gridSpacing1[1] = 1.0; gridSpacing1[2] = 1.0;
    normalOffset10 = 10;
    axisDataZ.insert(axisDataZ.begin(), 3, axis_t::Z);
    thinPlateKernel = ThinPlate;
    gaussianKernel = Gaussian;
    multiQuadKernel = MultiQuadratic;
  }

  // virtual void TearDown() {}
  std::vector<vec3> pointsTriangleClockwise; // triangle
  vec3 origin0; // (0, 0, 0)
  vec3 gridSize50; // (50, 50, 50)
  vec3 gridSpacing1; // (1, 1, 1)
  double normalOffset10;
  std::vector<axis_t> axisDataZ; // #axis entries min 3
  Kernel thinPlateKernel;
  Kernel gaussianKernel;
  Kernel multiQuadKernel;
};

TEST_F(RBFInterfaceTest, BasicInterfaceTestThinPlate)
{
  RBFInterface rbfInterface( pointsTriangleClockwise, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset10, axisDataZ,
                             false, false, thinPlateKernel );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
  const DataStorage rasterData = rbfInterface.getRasterData();
  // at least check that values in threshold range were generated
  // TODO: check linear system numerics
  int counter = 0;
  for (size_t i = 0; i < rasterData.size(); ++i)
  {
    for (size_t j = 0; j < rasterData[i].size(); ++j)
    {
      for (size_t k = 0; k < rasterData[i][j].size(); ++k)
      {
        if ( rasterData[i][j][k] > threshold )
          ++counter;
      }
    }
  }
  EXPECT_GT(counter, 0);
}

TEST_F(RBFInterfaceTest, BasicInterfaceTestGaussian)
{
  RBFInterface rbfInterface( pointsTriangleClockwise, origin0,
                            gridSize50, gridSpacing1,
                            normalOffset10, axisDataZ,
                            false, false, gaussianKernel );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
  const DataStorage rasterData = rbfInterface.getRasterData();
  // at least check that values in threshold range were generated
  // TODO: check linear system numerics
  int counter = 0;
  for (size_t i = 0; i < rasterData.size(); ++i)
  {
    for (size_t j = 0; j < rasterData[i].size(); ++j)
    {
      for (size_t k = 0; k < rasterData[i][j].size(); ++k)
      {
        if ( rasterData[i][j][k] > threshold )
          ++counter;
      }
    }
  }
  EXPECT_GT(counter, 0);
}

TEST_F(RBFInterfaceTest, BasicInterfaceTestMultiQuad)
{
  RBFInterface rbfInterface( pointsTriangleClockwise, origin0,
                            gridSize50, gridSpacing1,
                            normalOffset10, axisDataZ,
                            false, false, multiQuadKernel );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
  const DataStorage rasterData = rbfInterface.getRasterData();
  // at least check that values in threshold range were generated
  // TODO: check linear system numerics
  int counter = 0;
  for (size_t i = 0; i < rasterData.size(); ++i)
  {
    for (size_t j = 0; j < rasterData[i].size(); ++j)
    {
      for (size_t k = 0; k < rasterData[i][j].size(); ++k)
      {
        if ( rasterData[i][j][k] > threshold )
          ++counter;
      }
    }
  }
  EXPECT_GT(counter, 0);
}