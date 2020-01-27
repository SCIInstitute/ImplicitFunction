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
#include <iomanip>

#include "RBFInterface.h"
#include "vec3.h"

typedef std::vector< std::string > ViewModeList;

class Seg3DIntegrationTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    view_modes_.assign(modelPointData.size(), "axial");

    for ( const auto& mode : view_modes_ )
    {
      if ( mode == "sagittal" ) // X
      {
        axisData_.push_back(axis_t::X);
      }
      else if (mode == "coronal" ) // Y
      {
        axisData_.push_back(axis_t::Y);
      }
      else if ( mode == "axial" ) // Z
      {
        axisData_.push_back(axis_t::Z);
      }
      else
      {
        throw "Invalid viewer mode";
      }
    }
  }

  //TODO: loop this [0,20]
  double normalOffset_ {2.0};
  //TODO: loop over these bools
  bool compute2DConvexHull_ {false}, invertSeedOrder_ {false};
  //TODO: loop over kernels
  //enum Kernel { Gaussian, ThinPlate, MultiQuadratic };
  Kernel kernel_ {ThinPlate};

  std::vector<vec3> modelPointData = {
    {329.27408562074709, -3.1975149796730840, 477.29595999999998 }
    , {312.11090641613799, -17.909120331180826, 477.29595999999998 }
    , {311.62052986743487, -40.630377485176155, 477.29595999999998 }
    , {324.86069668241902, -53.707360019849709, 477.29595999999998 }
    , {326.65874402766383, -46.515019625779246, 477.29595999999998 }
    , {324.69723783285133, -25.918772133668394, 477.29595999999998 }
    , {311.78398871700256, -39.486141513392205, 478.19439699999998 }
    , {313.09165951354419, -19.380280866331621, 478.19439699999998 }
    , {330.58175641728877, -4.6686755148238497, 478.19439699999998 }
    , {325.35107323112214, -25.264923006934715, 478.19439699999998 }
    , {326.65874402766383, -44.880396808945051, 478.19439699999998 }
    , {323.71648473544508, -53.543897738166287, 478.19439699999998 }
    , {329.92792101901790, -5.8129114866077956, 479.09283399999998 }
    , {313.90895376138275, -21.341828246532643, 479.09283399999998 }
    , {313.58203606224731, -36.707282724774082, 479.09283399999998 }
    , {323.71648473544508, -49.784265259447636, 479.09283399999998 }
    , {326.82220287723152, -43.409236273794278, 479.09283399999998 }
    , {324.86069668241902, -25.918772133668394, 479.09283399999998 }
  };

  ViewModeList view_modes_;
  std::vector<axis_t> axisData_;

  vec3 modelOrigin_ {229.25999999999999, -126.41900000000000, 405.42099999999999 };
  vec3 modelGridSize_ {160.00000000000000, 232.00000000000000, 160.00000000000000 };
  vec3 modelGridSpacing_ {1.0000000000000000, 0.89843700000000004, 0.89843700000000004 };
};

TEST_F(Seg3DIntegrationTest, ImplicitModel)
{
  RBFInterface modelAlgo( modelPointData,
    modelOrigin_, modelGridSize_, modelGridSpacing_,
    normalOffset_, axisData_,
    compute2DConvexHull_, invertSeedOrder_, kernel_ );

  auto thresholdValue = modelAlgo.getThresholdValue();

  EXPECT_EQ(0.0, thresholdValue);

  const auto rasterData = modelAlgo.getRasterData();

  EXPECT_EQ(rasterData.size(), 160);
  EXPECT_EQ(rasterData[0].size(), 232);
  EXPECT_EQ(rasterData[0][0].size(), 160);

  EXPECT_NEAR(rasterData[0][0][0], -2701.018664769188, 1.0e-7);
  EXPECT_NEAR(rasterData[3][12][80], -1880.7781485604355, 1.0e-7);
  EXPECT_NEAR(rasterData[7][31][92], -1491.0771375730401, 1.0e-7);
  EXPECT_NEAR(rasterData[39][57][39], -604.78687230695505, 1.0e-7);
  EXPECT_NEAR(rasterData[79][115][79], -14.501352507455977, 1.0e-7);
  EXPECT_NEAR(rasterData[119][173][119], -631.7780323674524, 1.0e-7);
  EXPECT_NEAR(rasterData[159][231][159], -2885.3566775373765, 1.0e-7);

  //TODO: convert to move semantics for seg3d datablock usage
  #if 0
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
  #endif
}
