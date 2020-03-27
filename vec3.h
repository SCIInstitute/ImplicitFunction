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


#ifndef VEC3_H
#define VEC3_H


#include <iostream>


#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif


class vec3
{
public:
  vec3();
  vec3(double x, double y, double z);

private:
  double data_[3];

public:
  bool operator!=(const vec3 &a) const;
  bool operator==(const vec3 &a) const;
  bool operator<=(const vec3 &a) const;
  bool operator>=(const vec3 &a) const;
  bool operator<(const vec3 &a) const;
  bool operator>(const vec3 &a) const;
  vec3& operator=(const vec3 &a);
  vec3& operator+=(const vec3 &a);
  vec3& operator*=(double c);
  vec3& operator/=(double c);

  const double* getPtr() const { return data_; }

  double& operator[](const size_t);
  double  operator[](const size_t) const;

  double dot(const vec3 &b) const;
  vec3 cross(const vec3 &b);

  const double x() const { return data_[0]; }
  double& x() { return data_[0]; }
  void x( double b ) { data_[0] = b; }

  const double y() const { return data_[1]; }
  double& y() { return data_[1]; }
  void y( double b ) { data_[1] = b; }

  const double z() const { return data_[2]; }
  double& z() { return data_[2]; }
  void z( double b ) { data_[2] = b; }

  static vec3 zero;
  static vec3 unitX;
  static vec3 unitY;
  static vec3 unitZ;
  static vec3 min(const vec3 &a, const vec3 &b);
  static vec3 max(const vec3 &a, const vec3 &b);
  static const double PI;

  std::string toString() const;
  friend std::ostream& operator<<(std::ostream &stream, const vec3 &v);
};

vec3 cross(const vec3 &a, const vec3 &b);
vec3 bisect(const vec3 &a, const vec3 &b);
double dot(const vec3 &a, const vec3 &b);
double length(const vec3 &a);
double lengthSquared(const vec3 &a);
double distance(const vec3 &a, const vec3 &b);
double L1(const vec3 &a);
double L2(const vec3 &a);
double clamp(double value, double min, double max);
vec3 normalize(const vec3 &v1);
vec3 normalize(const vec3 &v1, float epsilon);

double vec2polar(const vec3 &a);

double angleBetween(const vec3 &a, const vec3 &b);
double angleBetween(const vec3 &a, const vec3 &b, float epsilon);

vec3 operator+(const vec3 &a, const vec3 &b);
vec3 operator-(const vec3 &a, const vec3 &b);
vec3 operator*(const vec3 &a, double b);
vec3 operator*(double a, const vec3 &b);
vec3 operator/(const vec3 &a, double b);

#endif // VEC3_H
