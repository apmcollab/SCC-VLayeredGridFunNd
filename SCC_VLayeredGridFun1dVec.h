/*
 * SCC_VLayeredGridFun1dVec.h
 *
 *  Created on: Nov 14, 2022
 *      Author: anderson
 */

//
// A class derived from VLayeredGridFun1d
// that provides member functions required for
// operators templated with respect to a mathematical
// vector class.
//

// In this derived class the getDimesion() member function
// returns the number of undetermined function values
// for an instance of VLayeredGridFun1d that is continuous
// across layer boundaries.
//
// To facilitate creating matrix representations of
// linear operators, the member function initializeBasisVector(...)
// is used to initialize the ith unit basis function for
// VLayeredGridFun1d that is continuous across layer boundaries.
//

#include "SCC_VLayeredGridFun1d.h"

#ifndef  VLAYERED_GRID_FUN_1D_VEC_
#define  VLAYERED_GRID_FUN_1D_VEC_

namespace SCC
{
class VLayeredGridFun1dVec : public SCC::VLayeredGridFun1d
{
public:

	VLayeredGridFun1dVec() : VLayeredGridFun1d()
	{}

	VLayeredGridFun1dVec(const VLayeredGridFun1dVec& V) : VLayeredGridFun1d(V)
	{}

	VLayeredGridFun1dVec(long zLayerCount, std::vector<long> zPanels, std::vector<double> zBdrys)
	: VLayeredGridFun1d(zLayerCount, zPanels, zBdrys){}


	long getDimension()
	{
		return getZpanelCountSum() + 1;
	}

	void initializeBasisVector(long i)
	{
		long zIndexMin = 0;
		long zIndexMax = 0;
		long zIndex;

		setToValue(0.0);

		if(i == 0)
		{
			layer[0](0) = 1.0;
			return;
		}

		for(long k = 0; k < layerCount; k++)
		{
			zIndexMin = zIndexMax;
			zIndexMax = zIndexMin + zPanels[k];
			if((zIndexMin + 1 <= i)&&(i <= zIndexMax))
			{
				zIndex = i-zIndexMin;
				layer[k](zIndex) = 1.0;

				// ensure continuity

				if((i == zIndexMax) && (k < layerCount-1))
				{
				layer[k+1](0) = 1.0;
				}
				return;
			}
		}
	}
};

}





#endif /*  VLAYERED_GRID_FUN_1D_VEC_  */
