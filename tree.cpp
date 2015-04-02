#include "tree.h"
#include "RBF.h"
#include "vec3.h"

double StructureTree::computeValue(vec3 x)
{
	if(fault!=NULL)
		return fault->computeValue(x);
	else if (horizon!=NULL)
		return horizon->computeValue(x);
	else if (tear!=NULL)
		return tear->computeValue(x);
	return 0;
}
