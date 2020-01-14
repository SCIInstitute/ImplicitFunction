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

class Seg3DIntegrationTest : public ::testing::Test
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

TEST_F(Seg3DIntegrationTest, ImplicitModel)
{
  FAIL() << "todo";
  #if 0
  DataLayerHandle srcDataLayer = boost::dynamic_pointer_cast<DataLayer>(this->actionInternal_->srcLayer_);
  DataLayerHandle dstDataLayer = boost::dynamic_pointer_cast<DataLayer>(this->actionInternal_->dstLayer_);
  GridTransform srcGridTransform = srcDataLayer->get_grid_transform();

  std::vector<vec3> modelPointData;
  for ( auto &vertex : this->actionInternal_->vertices_ )
  {
    modelPointData.push_back( vec3(vertex.x(), vertex.y(), vertex.z()) );
  }

  std::vector<axis_t> axisData;
  for ( auto &mode : this->actionInternal_->view_modes_ )
  {
    if ( mode == Viewer::SAGITTAL_C ) // X
    {
      axisData.push_back(axis_t::X);
    }
    else if (mode == Viewer::CORONAL_C ) // Y
    {
      axisData.push_back(axis_t::Y);
    }
    else if ( mode == Viewer::AXIAL_C ) // Z
    {
      axisData.push_back(axis_t::Z);
    }
    else
    {
      std::ostringstream oss;
      oss << "Invalid viewer mode " << mode;
      this->report_error( oss.str() );
      return;
    }
  }

  // origin and size from source data layer
  Point origin = srcGridTransform.get_origin();
  vec3 modelOrigin(origin.x(), origin.y(), origin.z());
  vec3 modelGridSize(srcGridTransform.get_nx(), srcGridTransform.get_ny(), srcGridTransform.get_nz());
  vec3 modelGridSpacing(srcGridTransform.spacing_x(), srcGridTransform.spacing_y(), srcGridTransform.spacing_z());

  // From RBF class. ThinPlate is the default kernel.
  Kernel kernel = ThinPlate;
  if (this->actionInternal_->kernel_ == "gaussian")
  {
    kernel = Gaussian;
  }
  else if (this->actionInternal_->kernel_ == "multi_quadratic")
  {
    kernel = MultiQuadratic;
  }

    RBFInterface modelAlgo( modelPointData, modelOrigin, modelGridSize, modelGridSpacing,
                          this->actionInternal_->normalOffset_, axisData,
                          this->actionInternal_->compute2DConvexHull_,
                          this->actionInternal_->invertSeedOrder_, kernel );

    this->actionInternal_->thresholdValue_ = modelAlgo.getThresholdValue();

    Core::DataBlockHandle dstDataBlock = Core::StdDataBlock::New( srcGridTransform, Core::DataType::DOUBLE_E );
    if ( ! dstDataBlock )
    {
      this->report_error( "Could not allocate enough memory." );
      return;
    }

    const DataStorage rasterData = modelAlgo.getRasterData();
    for (size_t i = 0; i < dstDataBlock->get_nx(); ++i)
    {
      for (size_t j = 0; j < dstDataBlock->get_ny(); ++j)
      {
        for (size_t k = 0; k < dstDataBlock->get_nz(); ++k)
        {
          dstDataBlock->set_data_at( i, j, k, rasterData[i][j][k] );
        }
      }
    }
    dstDataBlock->update_histogram();

#endif

  #if 0
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
        double d = std::stod(line);
      }
    }
  }
  in.close();
  #endif
}
