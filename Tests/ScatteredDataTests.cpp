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

#include "ScatteredData.h"

class ScatteredDataTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    for ( size_t i = 0; i < 10; ++i )
    {
      xCoords10.push_back( static_cast<double>(i) );
      yCoords10.push_back( static_cast<double>(i) );
      zCoords10.push_back( static_cast<double>(i) );
    }

    func10.insert(func10.begin(), 10, 0);
    axisInfo10.insert(axisInfo10.begin(), 10, axis_t::Z);
  }

// virtual void TearDown() {}

  std::vector<double> xCoords10, yCoords10, zCoords10, func10;
  std::vector<axis_t> axisInfo10;
};

// TODO: what happens if vector sizes don't match???

TEST_F(ScatteredDataTest, DefaultConstructor)
{
  ScatteredData data;
  ASSERT_EQ( data.surfacePoints_[0].size(), 0 );
  ASSERT_EQ( data.surfacePoints_[1].size(), 0 );
  ASSERT_EQ( data.surfacePoints_[2].size(), 0 );
  ASSERT_EQ( data.leftovers_[0].size(), 0 );
  ASSERT_EQ( data.leftovers_[1].size(), 0 );
  ASSERT_EQ( data.leftovers_[2].size(), 0 );
  ASSERT_EQ( data.fnc_.size(), 0 );
  ASSERT_EQ( data.axisInformation_.size(), 0 );
  ASSERT_EQ( data.updatedAxisInformation_.size(), 0 );
  ASSERT_EQ( data.origSize_, 0 );
}

TEST_F(ScatteredDataTest, SimpleSetup10)
{
  ScatteredData data( xCoords10, yCoords10, zCoords10, func10, axisInfo10 );
  ASSERT_EQ( data.surfacePoints_[0].size(), 10 );
  ASSERT_EQ( data.surfacePoints_[1].size(), 10 );
  ASSERT_EQ( data.surfacePoints_[2].size(), 10 );
  ASSERT_EQ( data.leftovers_[0].size(), 0 );
  ASSERT_EQ( data.leftovers_[1].size(), 0 );
  ASSERT_EQ( data.leftovers_[2].size(), 0 );
  ASSERT_EQ( data.fnc_.size(), 10 );
  ASSERT_EQ( data.axisInformation_.size(), 10 );
  ASSERT_EQ( data.updatedAxisInformation_.size(), 10 );
  ASSERT_EQ( data.origSize_, 10 );

  for (size_t i = 0; i < 10; ++i)
  {
    ASSERT_EQ( data.surfacePoints_[0][i], static_cast<double>(i) );
    ASSERT_EQ( data.surfacePoints_[1][i], static_cast<double>(i) );
    ASSERT_EQ( data.surfacePoints_[2][i], static_cast<double>(i) );
    ASSERT_EQ( data.fnc_[i], 0 );
    ASSERT_EQ( data.axisInformation_[i], axis_t::Z );
    ASSERT_EQ( data.updatedAxisInformation_[i], axis_t::Z );
  }
}

TEST_F(ScatteredDataTest, SimpleSetData10)
{
  ScatteredData data( xCoords10, yCoords10, zCoords10, func10, axisInfo10 );
  ASSERT_EQ( data.surfacePoints_[0].size(), 10 );
  ASSERT_EQ( data.surfacePoints_[1].size(), 10 );
  ASSERT_EQ( data.surfacePoints_[2].size(), 10 );
  ASSERT_EQ( data.leftovers_[0].size(), 0 );
  ASSERT_EQ( data.leftovers_[1].size(), 0 );
  ASSERT_EQ( data.leftovers_[2].size(), 0 );
  ASSERT_EQ( data.fnc_.size(), 10 );
  ASSERT_EQ( data.axisInformation_.size(), 10 );
  ASSERT_EQ( data.updatedAxisInformation_.size(), 10 );
  ASSERT_EQ( data.origSize_, 10 );

  for (size_t i = 0; i < 10; ++i)
  {
    ASSERT_EQ( data.surfacePoints_[0][i], static_cast<double>(i) );
    ASSERT_EQ( data.surfacePoints_[1][i], static_cast<double>(i) );
    ASSERT_EQ( data.surfacePoints_[2][i], static_cast<double>(i) );
    ASSERT_EQ( data.fnc_[i], 0 );
    ASSERT_EQ( data.axisInformation_[i], axis_t::Z );
    ASSERT_EQ( data.updatedAxisInformation_[i], axis_t::Z );
  }

  double newVal = 501.0;
  std::vector<double> newXCoords10(xCoords10);
  newXCoords10[5] = newVal;

  data.setData( newXCoords10, yCoords10, zCoords10, func10 );
  ASSERT_EQ( data.surfacePoints_[0].size(), 10 );
  ASSERT_EQ( data.surfacePoints_[1].size(), 10 );
  ASSERT_EQ( data.surfacePoints_[2].size(), 10 );
  ASSERT_EQ( data.leftovers_[0].size(), 0 );
  ASSERT_EQ( data.leftovers_[1].size(), 0 );
  ASSERT_EQ( data.leftovers_[2].size(), 0 );
  ASSERT_EQ( data.fnc_.size(), 10 );
  ASSERT_EQ( data.axisInformation_.size(), 10 );
  ASSERT_EQ( data.updatedAxisInformation_.size(), 10 );
  ASSERT_EQ( data.origSize_, 10 );

  for (size_t i = 0; i < 10; ++i)
  {
    if (i == 5)
    {
      ASSERT_EQ( data.surfacePoints_[0][i], newVal );
    }
    else
    {
      ASSERT_EQ( data.surfacePoints_[0][i], static_cast<double>(i) );
    }
    ASSERT_EQ( data.surfacePoints_[1][i], static_cast<double>(i) );
    ASSERT_EQ( data.surfacePoints_[2][i], static_cast<double>(i) );
    ASSERT_EQ( data.fnc_[i], 0 );
    ASSERT_EQ( data.axisInformation_[i], axis_t::Z );
    ASSERT_EQ( data.updatedAxisInformation_[i], axis_t::Z );
  }
}