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

#include "vec3.h"

class Vec3Test : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    sampleVec1[0] = 1.0, sampleVec1[1] = 2.0, sampleVec1[2] = 3.0;
    sampleVec2[0] = 1.5, sampleVec2[1] = 2.0, sampleVec2[2] = 3.0;
    sampleVec3[0] = 1.01, sampleVec3[1] = 2.01, sampleVec3[2] = 3.01;
  }

  // virtual void TearDown() {}
  vec3 emptyVec, sampleVec1, sampleVec2, sampleVec3;
};

TEST_F(Vec3Test, EmptyVec3)
{
  ASSERT_EQ( emptyVec[0], 0 );
  ASSERT_EQ( emptyVec[1], 0 );
  ASSERT_EQ( emptyVec[2], 0 );
}

TEST_F(Vec3Test, SampleVec3)
{
  ASSERT_EQ( sampleVec1[0], 1.0 );
  ASSERT_EQ( sampleVec1[1], 2.0 );
  ASSERT_EQ( sampleVec1[2], 3.0 );
}

TEST_F(Vec3Test, NotEqual)
{
  ASSERT_TRUE( sampleVec1 != sampleVec2 );
}

TEST_F(Vec3Test, GreaterThanEq)
{
  ASSERT_TRUE( sampleVec2 >= sampleVec1 );
}

TEST_F(Vec3Test, GreaterThan)
{
  ASSERT_TRUE( sampleVec3 > sampleVec1 );
  ASSERT_FALSE( sampleVec2 > sampleVec1 );
}

TEST_F(Vec3Test, LessThanEq)
{
  ASSERT_TRUE( sampleVec1 <= sampleVec2 );
}

TEST_F(Vec3Test, LessThan)
{
  ASSERT_TRUE( sampleVec1 < sampleVec3 );
  ASSERT_FALSE( sampleVec1 < sampleVec2 );
}

