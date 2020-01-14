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
#include <fstream>
#include <sstream>

#include "RBFInterface.h"
#include "vec3.h"

// Note: empty point vector results in segfault
// TODO: should be an error...

class RBFInterfaceTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    pointsTriangleClockwise.push_back( vec3(9.02381, 36.8638, 10.0) );
    pointsTriangleClockwise.push_back( vec3(34.8002, 36.8638, 10.0) );
    pointsTriangleClockwise.push_back( vec3(20.2824, 14.4397, 10.0) );
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
                             false, thinPlateKernel );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
  const DataStorage rasterData = rbfInterface.getRasterData();

  // TODO: check linear system numerics
  std::ifstream in;
  std::ostringstream oss;
  oss << REGRESSION_DIR << "/pointsTriangleClockwise.txt";
  in.open( oss.str().c_str(), std::ios::binary );
  in.exceptions( std::ofstream::failbit | std::ofstream::badbit );

  for (size_t i = 0; i < gridSize50[0]; ++i)
  {
    for (size_t j = 0; j < gridSize50[1]; ++j)
    {
      for (size_t k = 0; k < gridSize50[2]; ++k)
      {
        std::string line;
        std::getline(in, line);
        //std::cerr << line << " vs " << rasterData[i][j][k] << std::endl;
        //double d = std::stod(line);
      }
    }
  }
  in.close();
}

TEST_F(RBFInterfaceTest, BasicInterfaceTestGaussian)
{
  RBFInterface rbfInterface( pointsTriangleClockwise, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset10, axisDataZ,
                             false, gaussianKernel );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
  const DataStorage rasterData = rbfInterface.getRasterData();

  // TODO: check linear system numerics
  std::ifstream in;
  std::ostringstream oss;
  oss << REGRESSION_DIR << "/pointsTriangleClockwise.txt";
  in.open( oss.str().c_str(), std::ios::binary );
  in.exceptions( std::ofstream::failbit | std::ofstream::badbit );

  for (size_t i = 0; i < gridSize50[0]; ++i)
  {
    for (size_t j = 0; j < gridSize50[1]; ++j)
    {
      for (size_t k = 0; k < gridSize50[2]; ++k)
      {
        std::string line;
        std::getline(in, line);
        //std::cerr << line << " vs " << rasterData[i][j][k] << std::endl;
        //double d = std::stod(line);
      }
    }
  }
  in.close();
}

TEST_F(RBFInterfaceTest, BasicInterfaceTestMultiQuad)
{
  RBFInterface rbfInterface( pointsTriangleClockwise, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset10, axisDataZ,
                             false, multiQuadKernel );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default
  const DataStorage rasterData = rbfInterface.getRasterData();

  // TODO: check linear system numerics
  std::ifstream in;
  std::ostringstream oss;
  oss << REGRESSION_DIR << "/pointsTriangleClockwise.txt";
  in.open( oss.str().c_str(), std::ios::binary );
  in.exceptions( std::ofstream::failbit | std::ofstream::badbit );

  for (size_t i = 0; i < gridSize50[0]; ++i)
  {
    for (size_t j = 0; j < gridSize50[1]; ++j)
    {
      for (size_t k = 0; k < gridSize50[2]; ++k)
      {
        std::string line;
        std::getline(in, line);
        //std::cerr << line << " vs " << rasterData[i][j][k] << std::endl;
        //double d = std::stod(line);
      }
    }
  }
  in.close();
}
