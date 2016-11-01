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
    sample1Vec[0] = 1.0, sample1Vec[1] = 2.0, sample1Vec[2] = 3.0;
    sample2Vec[0] = 1.5, sample2Vec[1] = 2.0, sample2Vec[2] = 3.0;
    sample3Vec[0] = 1.01, sample3Vec[1] = 2.01, sample3Vec[2] = 3.01;
  }

  // virtual void TearDown() {}
  vec3 emptyVec, sampleSetterVec, sample1Vec, sample2Vec, sample3Vec;
};

TEST_F(Vec3Test, EmptyVec)
{
  ASSERT_EQ( emptyVec[0], 0 );
  ASSERT_EQ( emptyVec[1], 0 );
  ASSERT_EQ( emptyVec[2], 0 );
}

TEST_F(Vec3Test, ZeroVec)
{
  ASSERT_EQ( vec3::zero[0], 0 );
  ASSERT_EQ( vec3::zero[1], 0 );
  ASSERT_EQ( vec3::zero[2], 0 );
}

TEST_F(Vec3Test, UnitXVec)
{
  ASSERT_EQ( vec3::unitX[0], 1 );
  ASSERT_EQ( vec3::unitX[1], 0 );
  ASSERT_EQ( vec3::unitX[2], 0 );
}

TEST_F(Vec3Test, UnitYVec)
{
  ASSERT_EQ( vec3::unitY[0], 0 );
  ASSERT_EQ( vec3::unitY[1], 1 );
  ASSERT_EQ( vec3::unitY[2], 0 );
}

TEST_F(Vec3Test, UnitZVec)
{
  ASSERT_EQ( vec3::unitZ[0], 0 );
  ASSERT_EQ( vec3::unitZ[1], 0 );
  ASSERT_EQ( vec3::unitZ[2], 1 );
}

TEST_F(Vec3Test, Sample1Vec)
{
  ASSERT_EQ( sample1Vec[0], 1.0 );
  ASSERT_EQ( sample1Vec[1], 2.0 );
  ASSERT_EQ( sample1Vec[2], 3.0 );
}

TEST_F(Vec3Test, Sample2VecGetters)
{
  ASSERT_EQ( sample2Vec.x(), 1.5 );
  ASSERT_EQ( sample2Vec.y(), 2.0 );
  ASSERT_EQ( sample2Vec.z(), 3.0 );
}

TEST_F(Vec3Test, Sample2VecSetters)
{
  sampleSetterVec.x( 141.0 );
  sampleSetterVec.y( 400.0 );
  sampleSetterVec.z( 999.98 );

  ASSERT_EQ( sampleSetterVec.x(), 141.0 );
  ASSERT_EQ( sampleSetterVec.y(), 400.0 );
  ASSERT_EQ( sampleSetterVec.z(), 999.98 );
}

TEST_F(Vec3Test, NotEqual)
{
  ASSERT_TRUE( sample1Vec != sample2Vec );
}

TEST_F(Vec3Test, GreaterThanEq)
{
  ASSERT_TRUE( sample2Vec >= sample1Vec );
}

TEST_F(Vec3Test, GreaterThan)
{
  ASSERT_TRUE( sample3Vec > sample1Vec );
  ASSERT_FALSE( sample2Vec > sample1Vec );
}

TEST_F(Vec3Test, LessThanEq)
{
  ASSERT_TRUE( sample1Vec <= sample2Vec );
}

TEST_F(Vec3Test, LessThan)
{
  ASSERT_TRUE( sample1Vec < sample3Vec );
  ASSERT_FALSE( sample1Vec < sample2Vec );
}

