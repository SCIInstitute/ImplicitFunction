#include "BBox.h"

#include <limits>

BBox::BBox()
{
    empty = 1;
}

BBox::~BBox()
{

}

void BBox::setMin(vec3 _m)
{
    min = _m;
    isEmpty();
}

void BBox::setMin(float _m)
{
    min = vec3(_m,_m,_m);
    isEmpty();
}

void BBox::setMax(vec3 _m)
{
    max = _m;
    isEmpty();
}

void BBox::setMax(float _m)
{
    max = vec3(_m,_m,_m);
    isEmpty();
}

void BBox::reset()
{
    min = vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    max = vec3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
    isEmpty();
}

bool BBox::isEmpty()
{
    empty = !(min.x < max.x && min.y < max.y && min.z < max.z);
    return empty;
}

bool BBox::inside(const vec3 &pos)
{
    return (pos.x >= min.x && pos.x <= max.x &&
            pos.y >= min.y && pos.y <= max.y &&
            pos.z >= min.z && pos.z <= max.z);
}
