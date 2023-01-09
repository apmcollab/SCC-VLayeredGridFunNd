/*
 * SCC_RandOpDirichletVLG2d.h
 *
 *  Created on: Nov 14, 2022
 *      Author: anderson
 */

//
// Randomizes VLayeredGridFun2d with homogeneous boundary values
// and, if specified, continuity of values across layer boundaries
//

#ifndef RAND_OP_DIRICHLET_VLG_2D_
#define RAND_OP_DIRICHLET_VLG_2D_

#include "SCC_VLayeredGridFun2d.h"
#include "RandOpNd/RandOp2d.h"
namespace SCC
{
class RandOpDirichletVLG2d
{
    public :

	RandOpDirichletVLG2d(bool continuityFlag = true)
	{
		this->continuityFlag = continuityFlag;
	}

	RandOpDirichletVLG2d(const RandOpDirichletVLG2d& R)
	{
		this->continuityFlag = R.continuityFlag;
	}

	void initialize(bool continuityFlag = true)
	{
		this->continuityFlag = continuityFlag;
	}

    virtual ~RandOpDirichletVLG2d(){}

    virtual void randomize(VLayeredGridFun2d& v)
    {
    	long layerCount = v.getLayerCount();
    	for(long k = 0; k < layerCount; k++)
    	{
    	randOp2d.randomize(v.layer[k]);
    	}

    	v.setBoundaryValues(0.0);

    	// ensure continuity

    	long zPanels;
    	long xPanels = v.getXpanelCount();

    	if(continuityFlag)
    	{
    		for(long k = 1; k <= v.getLayerCount()-1; k++)
    		{
    			for(long i = 0; i <= xPanels; i++)
    			{
    			zPanels           = v.layer[k-1].yPanels;
    			v.layer[k](i,0)   = v.layer[k-1](i,zPanels);
    			}
    		}
    	}
    }

    RandOp2d    randOp2d;
    bool continuityFlag;
};
} // namespace SCC

#endif /* RAND_OP_DIRICHLET_1D_VLG_ */
