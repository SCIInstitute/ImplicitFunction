#ifndef BBOX_H
#define BBOX_H
#include "vec3.h"

class BBox
{
public:
    BBox();
    ~BBox();
    vec3 getMin(){return min;};
    vec3 getMax(){return max;};

    void setMin(vec3 _m);
    void setMin(float _m);
    void setMax(vec3 _m);
    void setMax(float _m);
    void reset();
    bool isEmpty();
    bool inside(const vec3 &pos);
//private:
    vec3 min, max;
    bool empty;
};
#endif // BBOX_H
