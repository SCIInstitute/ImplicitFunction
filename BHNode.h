#ifndef BHNODE_H
#define BHNODE_H
#include "vec3.h"
#include <vector>
#include <limits>
#include "BBox.h"


class BHNode
{
public:
    BHNode();
    ~BHNode();
    bool leaf;
    double mass;
    double coeff;
    vec3 center;
    BHNode *nodes[8];
    BBox box;
    std::vector<int> pts;
    int index;

    double sum_coeff, sum_coeff_xi, sum_xi2;
};

#endif // BHNODE_H

