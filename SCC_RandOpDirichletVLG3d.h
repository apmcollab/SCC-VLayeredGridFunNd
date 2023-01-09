/*
 * SCC_RandOpDirichletVLG3d.h
 *
 *  Created on: Nov 14, 2022
 *      Author: anderson
 */

//
// Randomizes VLayeredGridFun3d with homogeneous boundary values
// and, if specified, continuity of values across layer boundaries
//

#ifndef RAND_OP_DIRICHLET_VLG_3D_
#define RAND_OP_DIRICHLET_VLG_3D_

#include "SCC_VLayeredGridFun3d.h"
#include "RandOpNd/RandOp3d.h"
namespace SCC
{
class RandOpDirichletVLG3d
{
    public :

	RandOpDirichletVLG3d(bool continuityFlag = true)
	{
		this->continuityFlag = continuityFlag;
	}

	RandOpDirichletVLG3d(const RandOpDirichletVLG3d& R)
	{
		this->continuityFlag = R.continuityFlag;
	}

	void initialize(bool continuityFlag = true)
	{
		this->continuityFlag = continuityFlag;
	}

    virtual ~RandOpDirichletVLG3d(){}

    virtual void randomize(VLayeredGridFun3d& v)
    {
    	long layerCount = v.getLayerCount();
    	for(long k = 0; k < layerCount; k++)
    	{
    	randOp3d.randomize(v.layer[k]);
    	}

    	v.setBoundaryValues(0.0);

    	// ensure continuity

    	long zPanels;
    	long xPanels = v.getXpanelCount();
    	long yPanels = v.getYpanelCount();

    	if(continuityFlag)
    	{
    		for(long k = 1; k <= v.getLayerCount()-1; k++)
    		{
    			for(long i = 0; i <= xPanels; i++)
    			{
    			for(long j = 0; j <= yPanels; j++)
    			{
    			zPanels             = v.layer[k-1].zPanels;
    			v.layer[k](i,j,0)   = v.layer[k-1](i,j,zPanels);
    			}}
    		}
    	}
    }

    RandOp3d    randOp3d;
    bool continuityFlag;
};
} // namespace SCC

#endif /* RAND_OP_DIRICHLET_1D_VLG_ */
