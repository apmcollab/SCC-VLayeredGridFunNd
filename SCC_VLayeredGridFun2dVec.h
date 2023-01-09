/*
 * SCC_VLayeredGridFun2dVec.h
 *
 *  Created on: Nov 15, 2022
 *      Author: anderson
 */

//
// A class derived from VLayeredGridFun2d
// that provides member functions required for
// operators templated with respect to a mathematical
// vector class.
//

// In this derived class the getDimesion() member function
// returns the number of undetermined function values
// for an instance of VLayeredGridFun2d that is continuous
// across layer boundaries.

// To facilitate creating matrix representations of
// linear operators, the member function initializeBasisVector(...)
// is used to initialize the ith unit basis function for
// VLayeredGridFun1d that is continuous across layer boundaries.
//
#include "SCC_VLayeredGridFun1dVec.h"
#include "SCC_VLayeredGridFun2d.h"


#ifndef  VLAYERED_GRID_FUN_2D_VEC_
#define  VLAYERED_GRID_FUN_2D_VEC_

namespace SCC
{
class VLayeredGridFun2dVec : public SCC::VLayeredGridFun2d
{
public:

	VLayeredGridFun2dVec() : VLayeredGridFun2d()
	{}

	VLayeredGridFun2dVec(const VLayeredGridFun2dVec& V) : VLayeredGridFun2d(V)
	{}

	VLayeredGridFun2dVec(long xPanels, double xMin, double xMax,
	long zLayerCount, std::vector<long> zPanels, std::vector<double> zBdrys)
	: VLayeredGridFun2d(xPanels,xMin,xMax,zLayerCount, zPanels, zBdrys){}


	long getDimension()
	{
		return (getZpanelCountSum() + 1)*(xPanels+1);
	}

	void initializeBasisVector(long bIndex)
	{
		SCC::VLayeredGridFun1dVec vTmp1d(layerCount, zPanels,zBdrys);

		long vDimension = vTmp1d.getDimension();
		long xIndex     = bIndex/vDimension;  // N.B, integer division
		long zIndex     = bIndex - (xIndex)*vDimension;

		vTmp1d.initializeBasisVector(zIndex);
    	setToValue(0.0);
    	setConstantXslice(xIndex,vTmp1d);
	}
};

}





#endif /*  VLAYERED_GRID_FUN_1D_VEC_  */
