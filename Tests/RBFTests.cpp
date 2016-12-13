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

#include "RBF.h"

class RBFTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    std::vector<double> xCoords10, yCoords10, zCoords10, func10;
    std::vector<axis_t> axisInfo10;
    for ( size_t i = 0; i < 10; ++i )
    {
      xCoords10.push_back( static_cast<double>(i) );
      yCoords10.push_back( static_cast<double>(i) );
      zCoords10.push_back( static_cast<double>(i) );
    }

    func10.insert(func10.begin(), 10, 0);
    axisInfo10.insert(axisInfo10.begin(), 10, axis_t::Z);
    data = new ScatteredData( xCoords10, yCoords10, zCoords10, func10, axisInfo10 );
  }

  virtual void TearDown()
  {
    delete data;
  }

  ScatteredData *data;
};

TEST_F(RBFTest, BasicSetup)
{
  // TODO: need to expose getters to test RBF state
  RBF rbf(data, ThinPlate);
}