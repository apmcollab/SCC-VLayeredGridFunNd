/*
 * SCC_VLayeredGridFun3dVec.h
 *
 *  Created on: Nov 15, 2022
 *      Author: anderson
 */

//
// A class derived from VLayeredGridFun3d
// that provides member functions required for
// operators templated with respect to a mathematical
// vector class.
//

// In this derived class the getDimesion() member function
// returns the number of undetermined function values
// for an instance of VLayeredGridFun3d that is continuous
// across layer boundaries.

// To facilitate creating matrix representations of
// linear operators, the member function initializeBasisVector(...)
// is used to initialize the ith unit basis function for
// VLayeredGridFun1d that is continuous across layer boundaries.
//
#include "SCC_VLayeredGridFun1dVec.h"
#include "SCC_VLayeredGridFun3d.h"


#ifndef  VLAYERED_GRID_FUN_3D_VEC_
#define  VLAYERED_GRID_FUN_3D_VEC_

namespace SCC
{
class VLayeredGridFun3dVec : public SCC::VLayeredGridFun3d
{
public:

	VLayeredGridFun3dVec() : VLayeredGridFun3d()
	{}

	VLayeredGridFun3dVec(const VLayeredGridFun3dVec& V) : VLayeredGridFun3d(V)
	{}

	VLayeredGridFun3dVec(long xPanels, double xMin, double xMax,long yPanels,double yMin, double yMax,
	long zLayerCount, std::vector<long> zPanels, std::vector<double> zBdrys)
	: VLayeredGridFun3d(xPanels,xMin,xMax,yPanels,yMin,yMax,zLayerCount, zPanels, zBdrys){}


	long getDimension()
	{
		return (getZpanelCountSum() + 1)*(xPanels+1)*(yPanels+1);
	}

	void initializeBasisVector(long bIndex)
	{
		SCC::VLayeredGridFun1dVec vTmp1d(layerCount, zPanels,zBdrys);

		long vDimension  = vTmp1d.getDimension();
		long xyIndex     = bIndex/vDimension;  // N.B, integer division
		long zIndex      = bIndex - (xyIndex)*vDimension;

		//xyIndex = yIndex  + xIndex*(yPanels+1);

		long xIndex       = xyIndex/(yPanels+1);
		long yIndex       = xyIndex - (xIndex)*(yPanels+1);

		vTmp1d.initializeBasisVector(zIndex);
    	setToValue(0.0);
    	setConstantXYslice(xIndex,yIndex,vTmp1d);
	}
};

}





#endif /*  VLAYERED_GRID_FUN_3D_VEC_  */
