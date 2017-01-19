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
//#include <iomanip>
//#include <limits>

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

    prism.push_back( vec3(  20.0, 30.0, 20.0) );
    prism.push_back( vec3(  30.0, 30.0, 20.0) );
    prism.push_back( vec3(  25.0, 20.0, 20.0) );
    prism.push_back( vec3(  20.0, 30.0, 30.0) );
    prism.push_back( vec3(  30.0, 30.0, 30.0) );
    prism.push_back( vec3(  25.0, 20.0, 30.0) );

    prismWithInternalPoint.push_back( vec3(  20.0, 30.0, 20.0) );
    prismWithInternalPoint.push_back( vec3(  30.0, 30.0, 20.0) );
    prismWithInternalPoint.push_back( vec3(  25.0, 20.0, 20.0) );
    prismWithInternalPoint.push_back( vec3(  20.0, 30.0, 30.0) );
    prismWithInternalPoint.push_back( vec3(  30.0, 30.0, 30.0) );
    prismWithInternalPoint.push_back( vec3(  25.0, 20.0, 30.0) );
    prismWithInternalPoint.push_back( vec3(  25.0, 25.0, 25.0) );


    // prism with point just outside convex hull (v. close to triangle)

    gridSize50[0] = 50.0; gridSize50[1] = 50.0; gridSize50[2] = 50.0;
    gridSpacing1[0] = 1.0; gridSpacing1[1] = 1.0; gridSpacing1[2] = 1.0;
    normalOffset2 = 2.0;
    axisDataZ.insert(axisDataZ.begin(), 3, axis_t::Z);
    thinPlateKernel = ThinPlate;
    gaussianPlateKernel = Gaussian;
    multiQuadKernel = MultiQuadratic;

    nrrdHeader50_0_1.push_back("NRRD0001");
    nrrdHeader50_0_1.push_back("# Complete NRRD file format specification at:");
    nrrdHeader50_0_1.push_back("# http://teem.sourceforge.net/nrrd/format.html");
    nrrdHeader50_0_1.push_back("type: double");
    nrrdHeader50_0_1.push_back("dimension: 3");
    std::ostringstream sizes;
    sizes << "sizes: " << gridSize50[0] << " " << gridSize50[1] << " " << gridSize50[2];
    nrrdHeader50_0_1.push_back( sizes.str() );
    std::ostringstream mins;
    mins << "axis mins: " << origin0[0] << ", " << origin0[1] << ", " << origin0[2];
    nrrdHeader50_0_1.push_back( mins.str() );
    std::ostringstream spacings;
    spacings << "spacings: " << gridSpacing1[0] << " " << gridSpacing1[1] << " " << gridSpacing1[2];
    nrrdHeader50_0_1.push_back( spacings.str() );
    //nrrdHeader50_0_1.push_back("centerings: cell cell cell");
    nrrdHeader50_0_1.push_back("centerings: node node node");
    nrrdHeader50_0_1.push_back("endian: little");
    //nrrdHeader50_0_1.push_back("encoding: ascii");
    nrrdHeader50_0_1.push_back("encoding: raw");
  }

// virtual void TearDown() {}

//  std::vector<vec3> tet;
//  std::vector<vec3> cube;
  std::vector<vec3> prism;
  std::vector<vec3> prismWithInternalPoint;

  vec3 origin0; // (0, 0, 0)
  vec3 gridSize50; // (50, 50, 50)
  vec3 gridSpacing1; // (1, 1, 1)
  double normalOffset2;
  std::vector<axis_t> axisDataZ; // #axis entries min 3
  Kernel thinPlateKernel;
  Kernel gaussianPlateKernel;
  Kernel multiQuadKernel;

  std::vector< std::string > nrrdHeader50_0_1;
};

TEST_F(ConvexHull3DTest, PrismBasicInterfaceTest)
{
  RBFInterface rbfInterface( prism, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset2, axisDataZ,
                             false, ThinPlate );
  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default

  const DataStorage rasterData = rbfInterface.getRasterData();

//  std::string filename("testPrism.nrrd");
//  std::cout << "Writing file '" << filename << "'" << std::endl;
//  std::ofstream nrrdFile(filename.c_str(), std::ofstream::binary);
//  nrrdFile.exceptions( std::ofstream::failbit | std::ofstream::badbit );
//
//  if ( nrrdFile.is_open() )
//  {
////    nrrdFile << "NRRD0001" << std::endl;
////    nrrdFile << "# Complete NRRD file format specification at:" << std::endl;
////    nrrdFile << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
////    nrrdFile << "type: double" << std::endl;
////    nrrdFile << "dimension: 3" << std::endl;
////    nrrdFile << "sizes: " << gridSize50[0] << " " << gridSize50[1] << " " << gridSize50[2] << std::endl;
////    nrrdFile << "axis mins: " << origin0[0] << ", " << origin0[1] << ", " << origin0[2] << std::endl;
////    nrrdFile << "spacings: " << gridSpacing1[0] << " " << gridSpacing1[1] << " " << gridSpacing1[2] << std::endl;
//////    nrrdFile << "centerings: cell cell cell" << std::endl;
////    nrrdFile << "centerings: node node node" << std::endl;
////    nrrdFile << "endian: little" << std::endl;
////    nrrdFile << "encoding: raw" << std::endl;
//////    nrrdFile << "encoding: ascii" << std::endl;
////    nrrdFile << std::endl;
//
//    // write data portion
//    for (size_t i = 0; i < gridSize50[0]; ++i)
//    {
//      for (size_t j = 0; j < gridSize50[1]; ++j)
//      {
//        for (size_t k = 0; k < gridSize50[2]; ++k)
//        {
//          double val = rasterData[i][j][k];
//          nrrdFile.write((char*)&val, sizeof(double));
////          nrrdFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << val << " ";
//        }
//      }
//    }
//
//    nrrdFile.close();
//  }
//
//  std::ifstream baselineNrrdFile;
//  std::ostringstream oss;
//  oss << REGRESSION_DIR << "/testPrism.nrrd";
//  baselineNrrdFile.open( oss.str().c_str(), std::ios::binary );
//  baselineNrrdFile.exceptions( std::ofstream::failbit | std::ofstream::badbit );
//
//  std::string line;
//  if ( baselineNrrdFile.is_open() )
//  {
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//    std::cout << baselineNrrdFile.tellg() << std::endl;
//  }
}

TEST_F(ConvexHull3DTest, PrismWithInternalPointBasicInterfaceTest)
{
  RBFInterface rbfInterface( prismWithInternalPoint, origin0,
                             gridSize50, gridSpacing1,
                             normalOffset2, axisDataZ,
                             true, false, ThinPlate );

  const ScatteredData *data = rbfInterface.getSurfaceData();
  ASSERT_TRUE( data != nullptr );

//  // check scattered data setup, normal calculation for stability
//  std::vector<double> x = { -1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, -0.2265412489744628, -1.7734587510255371, 0.2058958485078267, 1.7941041514921734, 0.0, 0.0, -0.4384836277759048, -1.5615163722240952, 0.6029479320322761, 1.3970520679677239, 0.0, 0.0 };
//  std::vector<double> y = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.3203771272650137, -0.3203771272650137, -0.2325877093742827, 0.2325877093742827, 0.1835037736389938, 1.8164962263610063, 0.2325877133161922, -0.2325877133161922, 0.7259807942675245, -0.7259807942675245, 0.0571912277167201, 1.9428087722832799};
//  std::vector<double> z = { 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 5.0, 9.4530823069676728, 10.5469176930323272, 9.4384838105643816, 10.5615161894356184, 9.4226502293769805, 10.5773497706230195, 0.7941037281936285, -0.7941037281936285, 0.5615161784353572, -0.5615161784353572, 0.3333330950255075, -0.3333330950255075 };
//  std::vector<double> f = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, -10.0 };
//
//  ASSERT_EQ( x.size(), data->surfacePoints_[0].size() );
//  ASSERT_EQ( y.size(), data->surfacePoints_[1].size() );
//  ASSERT_EQ( z.size(), data->surfacePoints_[2].size() );
//  ASSERT_EQ( f.size(), data->fnc_.size() );
//
//  for (size_t i = 0; i < x.size(); ++i )
//  {
//    ASSERT_NEAR( x[i], data->surfacePoints_[0][i], 1.0e-10 );
//  }
//
//  for (size_t i = 0; i < y.size(); ++i )
//  {
//    ASSERT_NEAR( y[i], data->surfacePoints_[1][i], 1.0e-10 );
//  }
//
//  for (size_t i = 0; i < z.size(); ++i )
//  {
//    ASSERT_NEAR( z[i], data->surfacePoints_[2][i], 1.0e-10 );
//  }
//
//  for (size_t i = 0; i < f.size(); ++i )
//  {
//    ASSERT_NEAR( f[i], data->fnc_[i], 1.0e-10 );
//  }


  double threshold = rbfInterface.getThresholdValue();
  ASSERT_EQ( threshold, 0 ); // default

  const DataStorage rasterData = rbfInterface.getRasterData();

//  std::string filename("testPrismWithInternalPoint.nrrd");
//  std::cout << "Writing file '" << filename << "'" << std::endl;
//  std::ofstream nrrdFile(filename.c_str(), std::ofstream::binary);
//  nrrdFile.exceptions( std::ofstream::failbit | std::ofstream::badbit );
//
//  if ( nrrdFile.is_open() )
//  {
//    nrrdFile << "NRRD0001" << std::endl;
//    nrrdFile << "# Complete NRRD file format specification at:" << std::endl;
//    nrrdFile << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
//    nrrdFile << "type: double" << std::endl;
//    nrrdFile << "dimension: 3" << std::endl;
//    nrrdFile << "sizes: " << gridSize50[0] << " " << gridSize50[1] << " " << gridSize50[2] << std::endl;
//    nrrdFile << "axis mins: " << origin0[0] << ", " << origin0[1] << ", " << origin0[2] << std::endl;
//    nrrdFile << "spacings: " << gridSpacing1[0] << " " << gridSpacing1[1] << " " << gridSpacing1[2] << std::endl;
////    nrrdFile << "centerings: cell cell cell" << std::endl;
//    nrrdFile << "centerings: node node node" << std::endl;
//    nrrdFile << "endian: little" << std::endl;
//    nrrdFile << "encoding: raw" << std::endl;
//    nrrdFile << std::endl;
//
//    // write data portion
//    for (size_t i = 0; i < gridSize50[0]; ++i)
//    {
//      for (size_t j = 0; j < gridSize50[1]; ++j)
//      {
//        for (size_t k = 0; k < gridSize50[2]; ++k)
//        {
//          double val = rasterData[i][j][k];
//          nrrdFile.write((char*)&val, sizeof(double));
//        }
//      }
//    }
//
//    nrrdFile.close();
//  }
//
//  std::ifstream baselineNrrdFile;
//  std::ostringstream oss;
//  oss << REGRESSION_DIR << "/testPrismWithInternalPoint.nrrd";
//  baselineNrrdFile.open( oss.str().c_str(), std::ios::binary );
//  baselineNrrdFile.exceptions( std::ofstream::failbit | std::ofstream::badbit );
//
//  std::string line;
//  if ( baselineNrrdFile.is_open() )
//  {
//    std::getline(baselineNrrdFile, line);
//    std::cout << line << std::endl;
//  }
}

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