/*
 * SCC_RandOpDirichletVLG1d.h
 *
 *  Created on: Nov 14, 2022
 *      Author: anderson
 */

//
// Randomizes VLayeredGridFun1d with homogeneous boundary values
// and, if specified, continuity of values across layer boundaries
//

#ifndef RAND_OP_DIRICHLET_VLG_1D_
#define RAND_OP_DIRICHLET_VLG_1D_

#include "SCC_VLayeredGridFun1d.h"
#include "RandOpNd/RandOp1d.h"
namespace SCC
{
class RandOpDirichletVLG1d
{
    public :

	RandOpDirichletVLG1d(bool continuityFlag = true)
	{
		this->continuityFlag = continuityFlag;
	}

	RandOpDirichletVLG1d(const RandOpDirichletVLG1d& R)
	{
		this->continuityFlag = R.continuityFlag;
	}

	void initialize(bool continuityFlag = true)
	{
		this->continuityFlag = continuityFlag;
	}

    virtual ~RandOpDirichletVLG1d(){}

    virtual void randomize(VLayeredGridFun1d& v)
    {
    	long layerCount = v.getLayerCount();
    	for(long k = 0; k < layerCount; k++)
    	{
    	randOp1d.randomize(v.layer[k]);
    	}

    	v.setBoundaryValues(0.0);

    	// ensure continuity


    	long zPanels;
    	if(continuityFlag)
    	{
    		for(long k = 1; k <= v.getLayerCount()-1; k++)
    		{
    			zPanels         = v.layer[k-1].xPanels;
    			v.layer[k](0)   = v.layer[k-1](zPanels);
    		}
    	}
    }

    RandOp1d    randOp1d;
    bool continuityFlag;
};
} // namespace SCC

#endif /* RAND_OP_DIRICHLET_1D_VLG_ */
