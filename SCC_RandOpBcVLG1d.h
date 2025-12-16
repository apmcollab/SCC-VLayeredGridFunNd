/*
 * SCC_RandOpBcVLG1d.h
 *
 *  Created on:Dec 13, 2025
 *      Author: anderson
 */

//
// Randomizes VLayeredGridFun1d with either homogeneous boundary values,
// periodic boundary values, or unrestricted boundary values,
// and if specified, continuity of values across layer boundaries is enforced.
//
//  RandOpBC1d::DIRICHLET ==> sets grid boundary values to 0.0
//	RandOpBC1d::PERIODIC  ==> sets grid boundary values to be periodic
//	RandOpBC1d::NONE      ==> sets grid boundary values to be random values (no boundary constraint)
//
//

#include "SCC_VLayeredGridFun1d.h"
#include "RandOpNd/RandOpBC1d.h"

#ifndef RAND_OP_BC_VLG_1D_
#define RAND_OP_BC_VLG_1D_

namespace SCC
{
class RandOpBcVLG1d
{
    public :

	RandOpBcVLG1d(bool continuityFlag = true, int  Z_BC = RandOpBC1d::DIRICHLET)
	{
	initialize(continuityFlag,Z_BC);
	}

	RandOpBcVLG1d(const RandOpBcVLG1d& R)
	{
		this->continuityFlag = R.continuityFlag;
		this->Z_BC           = R.Z_BC;

		this->randOpBC1d.initialize(RandOpBC1d::NONE);
	}

	void initialize(bool continuityFlag = true,
	int  Z_BC = RandOpBC1d::DIRICHLET)
	{;
		this-> Z_BC = Z_BC;

		this->randOpBC1d.initialize(RandOpBC1d::NONE);

		this->continuityFlag = continuityFlag;
	}

    virtual ~RandOpBcVLG1d(){}

    virtual void randomize(VLayeredGridFun1d& v)
    {
    	long layerCount = v.getLayerCount();

    	for(long k = 0; k < layerCount; k++)
    	{
    	randOpBC1d.randomize(v.layer[k]);
    	}

    	long zPanels;

    	if(Z_BC == RandOpBC1d::DIRICHLET)
        {
    		v.layer[0](0) = 0.0;
    		zPanels  = v.layer[layerCount-1].xPanels;
    		v.layer[layerCount-1](zPanels) = 0.0;
        }

        if(Z_BC == RandOpBC1d::PERIODIC)
        {
        	zPanels                        = v.layer[layerCount-1].xPanels;
        	v.layer[layerCount-1](zPanels) = v.layer[0](0);
        }

    	// ensure continuity

    	if(continuityFlag)
    	{
    		for(long k = 1; k <= v.getLayerCount()-1; k++)
    		{
    			zPanels         = v.layer[k-1].xPanels;
    			v.layer[k](0)   = v.layer[k-1](zPanels);
    		}
    	}
    }

    int  Z_BC;

    RandOpBC1d    randOpBC1d;
    bool      continuityFlag;
};
} // namespace SCC

#endif /*  */
