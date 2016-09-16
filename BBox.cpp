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

#include "BBox.h"

#include <limits>

BBox::BBox() :
  empty_(true)
{
  reset();
}

BBox::~BBox()
{
}

void BBox::setMin(const vec3& m)
{
  this->min_ = m;
  isEmpty();
}

void BBox::setMin(float m)
{
  this->min_ = vec3(m, m, m);
  isEmpty();
}

void BBox::setMax(const vec3& m)
{
  this->max_ = m;
  isEmpty();
}

void BBox::setMax(float m)
{
  this->max_ = vec3(m, m, m);
  isEmpty();
}

void BBox::reset()
{
  this->min_ = vec3(std::numeric_limits<float>::max(),
                    std::numeric_limits<float>::max(),
                    std::numeric_limits<float>::max());
  this->max_ = vec3(std::numeric_limits<float>::min(),
                    std::numeric_limits<float>::min(),
                    std::numeric_limits<float>::min());
  isEmpty();
}

inline
bool BBox::isEmpty()
{
  this->empty_ = !(this->min_.x_ < this->max_.x_ &&
                   this->min_.y_ < this->max_.y_ &&
                   this->min_.z_ < this->max_.z_);
  return this->empty_;
}

inline
bool BBox::inside(const vec3 &pos)
{
  return (pos.x_ >= this->min_.x_ && pos.x_ <= this->max_.x_ &&
          pos.y_ >= this->min_.y_ && pos.y_ <= this->max_.y_ &&
          pos.z_ >= this->min_.z_ && pos.z_ <= this->max_.z_);
}
