/*
 * SCC_RandOpBcVLG3d.h
 *
 *  Created on:Dec 13, 2025
 *      Author: anderson
 */

//
// Randomizes VLayeredGridFun3d with either homogeneous boundary values,
// periodic boundary values, or unrestricted boundary values in each
// coordinate direction, and if specified, continuity of values
// across layer boundaries is enforced.
//
//  RandOpBC3d::DIRICHLET ==> sets grid boundary values to 0.0
//	RandOpBC3d::PERIODIC  ==> sets grid boundary values to be periodic
//	RandOpBC3d::NONE      ==> sets grid boundary values to be random values (no boundary constraint)
//
//
#include "SCC_VLayeredGridFun3d.h"
#include "RandOpNd/RandOpBC3d.h"

#ifndef RAND_OP_BC_VLG_3D_
#define RAND_OP_BC_VLG_3D_

namespace SCC
{
class RandOpBcVLG3d
{
    public :

    enum {DIRICHLET,PERIODIC,NONE}; // Enumeration consistent between all RandOpBCNd classes

	RandOpBcVLG3d(bool continuityFlag = true,
	int  X_BC = RandOpBcVLG3d::DIRICHLET,
	int  Y_BC = RandOpBcVLG3d::DIRICHLET,
	int  Z_BC = RandOpBcVLG3d::DIRICHLET)
	{
	initialize(continuityFlag,X_BC,Y_BC,Z_BC);
	}

	RandOpBcVLG3d(const RandOpBcVLG3d& R)
	{
		this->continuityFlag = R.continuityFlag;
		this->X_BC           = R.X_BC;
		this->Y_BC           = R.Y_BC;
		this->Z_BC           = R.Z_BC;

		this->randOpBC3d.initialize(X_BC,Y_BC,RandOpBcVLG3d::NONE);
	}

	void initialize(bool continuityFlag = true,
	int  X_BC = RandOpBcVLG3d::DIRICHLET,
	int  Y_BC = RandOpBcVLG3d::DIRICHLET,
	int  Z_BC = RandOpBcVLG3d::DIRICHLET)
	{
		this-> X_BC = X_BC;
		this-> Y_BC = Y_BC;
		this-> Z_BC = Z_BC;

		this->randOpBC3d.initialize(X_BC,Y_BC,RandOpBcVLG3d::NONE);

		this->continuityFlag = continuityFlag;
	}

    virtual ~RandOpBcVLG3d(){}

    virtual void randomize(VLayeredGridFun3d& v)
    {
    	long layerCount = v.getLayerCount();
    	for(long k = 0; k < layerCount; k++)
    	{
    	randOpBC3d.randomize(v.layer[k]);

        if(X_BC == RandOpBcVLG3d::DIRICHLET) {v.layer[k].setXboundaryValues(0.0);}
		if(Y_BC == RandOpBcVLG3d::DIRICHLET) {v.layer[k].setYboundaryValues(0.0);}

		if(X_BC == RandOpBcVLG3d::PERIODIC) {v.layer[k].enforceXperiodicity();}
		if(Y_BC == RandOpBcVLG3d::PERIODIC) {v.layer[k].enforceYperiodicity();}
    	}

    	long zPanels;
    	long xPanels    = v.getXpanelCount();
    	long yPanels    = v.getYpanelCount();

    	if(Z_BC == RandOpBcVLG3d::DIRICHLET)
        {
    		zPanels  = v.layer[layerCount-1].zPanels;
        	for(long i = 0; i <= xPanels; i++)
    		{
    		for(long j = 0; j <= yPanels; j++)
    		{
    			v.layer[0](i,j,0) = 0.0;
    			v.layer[layerCount-1](i,j,zPanels) = 0.0;
    		}}
        }
    	else if(Z_BC == RandOpBcVLG3d::PERIODIC)
        {
    		zPanels  = v.layer[layerCount-1].zPanels;
        	for(long i = 0; i <= xPanels; i++)
    		{
    		for(long j = 0; j <= yPanels; j++)
    		{

    			v.layer[layerCount-1](i,j,zPanels) = v.layer[0](i,j,0);
    		}}
        }

    	// ensure continuity

    	if(continuityFlag)
    	{
    		for(long k = 1; k <= v.getLayerCount()-1; k++)
    		{
    			zPanels             = v.layer[k-1].zPanels;
    			for(long i = 0; i <= xPanels; i++)
    			{
    			for(long j = 0; j <= yPanels; j++)
    			{
    			v.layer[k](i,j,0)   = v.layer[k-1](i,j,zPanels);
    			}}
    		}
    	}
    }

    int  X_BC; int  Y_BC; int  Z_BC;

    RandOpBC3d    randOpBC3d;
    bool      continuityFlag;
};
} // namespace SCC

#endif /*  */
