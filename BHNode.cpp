#include "BHNode.h"

BHNode::BHNode()
{
    leaf = 1;
    for (int i=0; i<8; i++){
        nodes[i] = 0;
    }
    mass = 1.0;
    center = vec3(0.0,0.0,0.0);
}

BHNode::~BHNode()
{

}
