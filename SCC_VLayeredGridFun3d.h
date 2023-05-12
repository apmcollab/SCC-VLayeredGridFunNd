
//
//####################################################################
//                    VLayeredGridFun3d.h
//####################################################################
/**
   An aggregate linear array of SCC::GridFunnction3d objects that is used
   to represented a layered structure stacked "vertically" with the 
   z-coordinate being the vertical direction  <p>

   Indexing for the layers starts at 0.


   The coordinates for the SCC::GridFunction3d instances that comprise
   the layers are (x,y,z), so the coordinate labels are identical
   to the coordinate labels of the layered grid coordinates.

*/
//####################################################################
// Chris Anderson (C) UCLA                            August 31, 2018
//####################################################################
//

/*
#############################################################################
#
# Copyright 2018 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/
#ifdef _MSC_VER
#include "iso646.h"          // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#undef min
#undef max
#endif


#include <functional>
#include <iostream>
#include <vector>


#include "GridFunctionNd/SCC_GridFunction3d.h"
#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "GridFunctionNd/SCC_GridFunction1d.h"

#include "SCC_VLayeredGridFun1d.h"
#include "SCC_VLayeredGridFun2d.h"


#ifndef SCC_VLAYERED_GRID_FUN3D_
#define SCC_VLAYERED_GRID_FUN3D_


namespace SCC
{

class VLayeredGridFun3d
{

public:

VLayeredGridFun3d()
{
	initialize();
}

VLayeredGridFun3d(long xPanels,    double xMin, double xMax,
		          long yPanels,    double yMin, double yMax,
				  long layerCount, const std::vector<long>& zPanels, const std::vector<double>& zBdrys)
{
    initialize(xPanels,xMin,xMax,yPanels,yMin,yMax,layerCount,zPanels,zBdrys);
}

VLayeredGridFun3d(const VLayeredGridFun3d& M)
{
    if(M.layer.size() == 0){initialize(); return;}
    initialize(M);
}

virtual void initialize()
{
	xWidth   = 1.0;
	xMin     = 0.0;
	xMax     = 1.0;
	xPanels  = 1;

	yWidth   = 1.0;
	yMin     = 0.0;
	yMax     = 1.0;
	yPanels  = 1;

	layerCount = 0;
	layer.clear();
    zWidth.clear();
    zBdrys.clear();
    zPanels.clear();
}

virtual void initialize(const VLayeredGridFun3d& M)
{
	xWidth   = M.xWidth;
    xMin     = M.xMin;
    xMax     = M.xMax;
    xPanels  = M.xPanels;

    yWidth   = M.yWidth;
    yMin     = M.yMin;
    yMax     = M.yMax;
    yPanels  = M.yPanels;

    layerCount = M.layerCount;

    // Can't use std::vector = for layer because the bounds checking
    // of the underlying arrays inhibits initializing with an
    // instance of a different size

    layer.clear();
    layer.resize(layerCount);
    for(long k = 0; k < layerCount; k++)
    {
    	layer[k].initialize(M.layer[k]);
	}

    zWidth     = M.zWidth;
    zBdrys     = M.zBdrys;
    zPanels    = M.zPanels;
}

virtual void initialize(long xPanels, double xMin, double xMax,
		        long yPanels, double yMin, double yMax,
			    long layerCount, const std::vector<long>& zPanels, const std::vector<double>& zBdrys)
{
    initialize();

	this->xWidth   = xMax-xMin;
    this->xMin     = xMin;
    this->xMax     = xMax;
    this->xPanels  = xPanels;

    this->yWidth   = yMax-yMin;
    this->yMin     = yMin;
    this->yMax     = yMax;
    this->yPanels  = yPanels;


    this->layerCount = layerCount;
    this->layer.resize(layerCount);

    this->zBdrys  = zBdrys;
    this->zPanels = zPanels;

    this->zWidth.resize(layerCount);

    double e; double f;


    for(long i = 0; i < layerCount; i++)
    {
    e = zBdrys[i];
    f = zBdrys[i+1];
    layer[i].initialize(xPanels,xMin,xMax,yPanels,yMin,yMax,zPanels[i],e,f);
    this->zWidth[i]  = f-e;
    }
}

virtual ~VLayeredGridFun3d()
{}




// Two structures are equal if they have the same number of
// layers and each layer is the same width and has the same
// number of panels. It is not required that the structures
// have identical x-boundary points.
//
virtual bool isEqualStructure(const VLayeredGridFun3d& V) const
{

    if((layer.size() == 0)&&(V.layer.size() != 0))       return false;
    if(layerCount != V.layerCount)                       return false;

    // transverse structure

    if(xPanels != V.xPanels){return false;}
    if(yPanels != V.yPanels){return false;}
    if(xWidth  != V.xWidth) {return false;}
    if(yWidth  != V.yWidth) {return false;}

    // vertical structure

    for(long i = 0; i < layerCount; i++)
    {
    if(zPanels[i] != V.zPanels[i]){return false;}
    else if(std::abs(zWidth[i] - V.zWidth[i]) > 1.0e-012*std::abs(zWidth[i])){return false;}
    }

    return true;
}


void addToValues(double val)
{
   for(long k = 0; k < layerCount; k++)
   {
    layer[k].addValue(val);
   }
}

void setToValue(double val)
{
   for(long k = 0; k < layerCount; k++)
   {
    layer[k].setToValue(val);
   }
}

void specify(std::function<double(double,double,double)> F)
{

	for(long k = 0; k < layerCount; k++)
	{
	layer[k].specify(F);
	}
}


VLayeredGridFun3d  extractLayers(long begIndex, long endIndex) const
{

	long  lcount = (endIndex - begIndex) + 1; 
	std::vector<long>   zPn(lcount);
	std::vector<double> zBd(lcount+1);

	for(long k = begIndex; k <= endIndex; k++)
	{
		zPn[k-begIndex] = layer[k].getZpanelCount();
	}

	zBd[0] = layer[begIndex].getZmin();
	for(long k = begIndex; k <= endIndex; k++)
	{
		zBd[k-begIndex+1] = layer[k].getZmax();
	}

   VLayeredGridFun3d R(xPanels,xMin,xMax,yPanels,yMin,yMax,lcount,zPn,zBd);

   for(long k = begIndex; k <= endIndex; k++)
   {
	R.layer[k-begIndex] = layer[k];
   }

   return R;
}


void  insertLayers(long begIndex, long endIndex, const VLayeredGridFun3d& V)
{
   for(long k = begIndex; k <= endIndex; k++)
   {
   layer[k] = V.layer[k-begIndex];
   }
}


void createProductFunction(const GridFunction1d& funX, const GridFunction1d& funY, const VLayeredGridFun1d& funZ)
{
    initialize(funX.xPanels, funX.xMin, funX.xMax,
               funY.xPanels, funY.xMin, funY.xMax,
               funZ.layerCount, funZ.zPanels, funZ.zBdrys);

	long i; long j; long k; long p;

	double fX; double fY; double fZ;

	for(p = 0; p < layerCount; p++)
	{

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX(i);
	for(j = 0; j <= yPanels; j++)
	{
	fY = funY(j);
	for(k = 0; k <= zPanels[p]; k++)
	{
			fZ = funZ.layer[p](k);
			layer[p](i,j,k) = fX*fY*fZ;
	}}}

	}
}

void createProductFunction(const GridFunction2d& funXY,  const VLayeredGridFun1d& funZ)
{
    initialize(funXY.xPanels, funXY.xMin, funXY.xMax,
               funXY.yPanels, funXY.yMin, funXY.yMax,
               funZ.layerCount, funZ.zPanels, funZ.zBdrys);

	long i; long j; long k; long p;

	double fXY; double fZ;

	for(p = 0; p < layerCount; p++)
	{

	for(i = 0; i <= xPanels; i++)
	{
    for(j = 0; j <= yPanels; j++)
	{
	fXY = funXY.Values(i,j);
	for(k = 0; k <= zPanels[p]; k++)
	{
	fZ = funZ.layer[p](k);

	layer[p](i,j,k) = fXY*fZ;
	}}}

	}
}

void createProductFunction(const GridFunction1d& funX,const  VLayeredGridFun2d& funYZ)
{
    initialize(funX.xPanels,  funX.xMin,  funX.xMax,
               funYZ.xPanels, funYZ.xMin, funYZ.xMax,
               funYZ.layerCount, funYZ.zPanels, funYZ.zBdrys);

	long i; long j; long k; long p;
	double fX; double fYZ;

	for(p = 0; p < layerCount; p++)
	{
	for(i = 0; i <= xPanels; i++)
	{
	fX = funX(i);
	for(j = 0; j <= yPanels; j++)
	{
	for(k = 0; k <= zPanels[p]; k++)
	{
	fYZ = funYZ.layer[p](j,k);
	layer[p](i,j,k) = fX*fYZ;
	}}}

	}
}
void  incrementLayers(long begIndex, long endIndex, const VLayeredGridFun3d& V)
{
   for(long k = begIndex; k <= endIndex; k++)
   {
   layer[k] += V.layer[k-begIndex];
   }
}

GridFunction2d getConstantZslice(long layerIndex, long zIndex) const //(x-y function)
{
    GridFunction2d R(xPanels,xMin,xMax,yPanels,yMin,yMax);
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    R.Values(i,j) = layer[layerIndex](i,j,zIndex);
    }}
    return R;
}

void getConstantZslice(long layerIndex, long zIndex,GridFunction2d& R) const //(x-y function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    R.Values(i,j) = layer[layerIndex](i,j,zIndex);
    }}
}

void setConstantZslice(long layerIndex, long zIndex, const GridFunction2d& R)  //(x-y function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long j = 0; j <= yPanels; j++)
    {
    layer[layerIndex](i,j,zIndex) = R.Values(i,j);
    }}
}


VLayeredGridFun2d getConstantYslice(long yIndex) const //(x-z function)
{
	VLayeredGridFun2d R(xPanels,xMin,xMax,layerCount,zPanels,zBdrys);

	for(long p = 0; p < layerCount; p++)
	{
		for(long i = 0; i <= xPanels; i++)
		{
			for(long k = 0; k <= zPanels[p]; k++)
			{
				R.layer[p](i,k) = layer[p](i,yIndex,k);
		}   }
    }
    return R;
}

void getConstantYslice(long yIndex,VLayeredGridFun2d& R) const //(x-z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long i = 0; i <= xPanels; i++)
		{
			for(long k = 0; k <= zPanels[p]; k++)
			{
				R.layer[p](i,k) = layer[p](i,yIndex,k);
		}   }
    }
}

void setConstantYslice(long yIndex,const VLayeredGridFun2d& R)  //(x-z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long i = 0; i <= xPanels; i++)
		{
			for(long k = 0; k <= zPanels[p]; k++)
			{
				layer[p](i,yIndex,k) = R.layer[p](i,k);
		}   }
    }
}


VLayeredGridFun2d getConstantXslice(long xIndex) const //(y-z function)
{
	VLayeredGridFun2d R(yPanels,yMin,yMax,layerCount,zPanels,zBdrys);

	for(long p = 0; p < layerCount; p++)
	{
		for(long j = 0; j <= yPanels; j++)
		{
			for(long k = 0; k <= zPanels[p]; k++)
			{
				R.layer[p](j,k) = layer[p](xIndex,j,k);
			}
		}
    }
    return R;
}

void getConstantXslice(long xIndex,VLayeredGridFun2d& R) const //(y-z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long j = 0; j <= yPanels; j++)
		{
			for(long k = 0; k <= zPanels[p]; k++)
			{
				R.layer[p](j,k) = layer[p](xIndex,j,k);
			}
		}
    }
}

void setConstantXslice(long xIndex,const VLayeredGridFun2d& R)  //(y-z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long j = 0; j <= yPanels; j++)
		{
			for(long k = 0; k <= zPanels[p]; k++)
			{
				layer[p](xIndex,j,k) = R.layer[p](j,k);
			}
		}
    }
}


GridFunction1d getConstantYZslice(long yIndex, long layerIndex, long zIndex) const  // (x function)
{
	GridFunction1d R(xPanels,xMin,xMax);

	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = layer[layerIndex](i,yIndex,zIndex);
	}
	return R;
}

void getConstantYZslice(long yIndex, long layerIndex, long zIndex,GridFunction1d& R) const  // (x function)
{
	for(long i = 0; i <= xPanels; i++)
	{
	R.Values(i) = layer[layerIndex](i,yIndex,zIndex);
	}
}

void setConstantYZslice(long yIndex, long layerIndex, long zIndex,const GridFunction1d& R) // (x function)
{
	for(long i = 0; i <= xPanels; i++)
	{
	layer[layerIndex](i,yIndex,zIndex) = R.Values(i);
	}
}


GridFunction1d getConstantXZslice(long xIndex, long layerIndex, long zIndex) const  //( y function)
{
	GridFunction1d R(yPanels,yMin,yMax);

	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = layer[layerIndex](xIndex,j,zIndex);
	}
	return R;
}

void getConstantXZslice(long xIndex, long layerIndex, long zIndex,GridFunction1d& R) const  //( y function)
{
	for(long j = 0; j <= yPanels; j++)
	{
	R.Values(j) = layer[layerIndex](xIndex,j,zIndex);
	}
}

void setConstantXZslice(long xIndex, long layerIndex, long zIndex, const GridFunction1d& R)   //( y function)
{
	for(long j = 0; j <= yPanels; j++)
	{
	layer[layerIndex](xIndex,j,zIndex) = R.Values(j);
	}
}


VLayeredGridFun1d getConstantXYslice(long xIndex, long yIndex) const  //( z function)
{

	VLayeredGridFun1d R(layerCount, zPanels,zBdrys);

	for(long p = 0; p < layerCount; p++)
	{
		for(long k = 0; k <= zPanels[p]; k++)
		{
			R.layer[p](k) = layer[p](xIndex,yIndex,k);
		}
	}
	return R;
}

void getConstantXYslice(long xIndex, long yIndex,VLayeredGridFun1d& R) const  //( z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long k = 0; k <= zPanels[p]; k++)
		{
			R.layer[p](k) = layer[p](xIndex,yIndex,k);
		}
	}
}

void  setConstantXYslice(long xIndex, long yIndex, const VLayeredGridFun1d& R) //( z function)
{

	for(long p = 0; p < layerCount; p++)
	{
		for(long k = 0; k <= zPanels[p]; k++)
		{
			layer[p](xIndex,yIndex,k) = R.layer[p](k);
		}
	}

}


long getLayerCount() const
{
	return layerCount;
}

std::vector<double> getZbdrys() const
{
	return zBdrys;
}

std::vector<long> getZpanels() const
{
	return zPanels;
}

double getXmin() const {return xMin;}
double getXmax() const {return xMax;}

double getYmin() const {return yMin;}
double getYmax() const {return yMax;}

double getHx() const {return xWidth/(double)xPanels;}
double getHy() const {return yWidth/(double)yPanels;}

double getZmin() const
{return layer[0].getZmin();}

double getZmax() const
{return layer[layerCount-1].getZmax();}

long getXpanelCount() const {return xPanels;}

long getYpanelCount() const {return yPanels;}

long getZpanelCountSum() const
{
	long panelSum = 0;
	for(long k = 0; k < layerCount; k++)
	{
		panelSum += layer[k].getZpanelCount();
	}
	return panelSum;
}

void operator=(const VLayeredGridFun3d& M)
{
	initialize(M);
}

void setLayerValues(long k,const SCC::GridFunction3d& L)
{
    layer[k] = (SCC::DoubleVector3d)L;
}

SCC::GridFunction3d  getLayer(long k)
{
    return layer[k];
}

SCC::GridFunction3d* getLayerPtr(long k)
{
    return &layer[k];
}

void  operator*=(const double alpha)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *= alpha;
    }
}

void operator*=(const VLayeredGridFun3d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *= M.layer[i];
    }
}

VLayeredGridFun3d operator-(const VLayeredGridFun3d& M)
{
    VLayeredGridFun3d R(M);
    for(long i = 0; i < layerCount; i++)
    {
    R.layer[i] = layer[i]- M.layer[i];
    }
    return R;
}

void operator-=(const VLayeredGridFun3d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] -= M.layer[i];
    }
}

void operator+=(const VLayeredGridFun3d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] += M.layer[i];
    }
}

VLayeredGridFun3d operator+(const VLayeredGridFun3d& M)
{
    VLayeredGridFun3d R(M);
    for(long i = 0; i < layerCount; i++)
    {
    R.layer[i] = layer[i] + M.layer[i];
    }
    return R;

}


void operator*=(const std::function<double(double,double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *=F;
    }
}

void operator/=(const std::function<double(double,double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] /=F;
    }
}

void operator+=(const std::function<double(double,double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] +=F;
    }
}

void operator-=(const std::function<double(double,double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] -=F;
    }
}


void squareValues()
{
//
//  Square the function values 
//
    for(long k = 0; k < layerCount; k++)
    {
        layer[k].squareValues();
    }
}

//
// dot == scaled dot == sum of scaled dot products
// of each layer SCC::GridFunction1d instance.
//
double dot(const VLayeredGridFun3d& M) const
{
    return scaledDot(M);
}



// Dot product in which the integral is evaluated using the
// Trapezoidal method.

double scaledDot(const VLayeredGridFun3d& M) const
{
    double dotVal = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    dotVal += layer[i].scaledDot(M.layer[i]);
    }
    return dotVal;
}

double normInf() const
{
    double infVal = 0.0;
    double infTmp;
    for(long i = 0; i < layerCount; i++)
    {
    infTmp = layer[i].normInf();
    infVal = (infTmp > infVal) ? infTmp : infVal;
    }
    return infVal;
}

// The required integral is approximated using the Trapezoidal method

double norm1() const
{
    double norm1val = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    norm1val += layer[i].norm1();
    }
    return norm1val;
}

// The required integral is approximated using the Trapezoidal method

double norm2() const
{
    double norm2val = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    norm2val += layer[i].norm2squared();
    }
    return sqrt(norm2val);
}

// The required integral is approximated using the Trapezoidal method

double integral() const
{
    double intVal = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    intVal += layer[i].integralTrapezoidal();
    }
    return intVal;
}

//
// Integral of F( grid function values )
//
double integral(std::function<double(double)> F) const
{
	double intVal = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    intVal += layer[i].integralTrapezoidal(F);
    }
    return intVal;
}

double min() const
{
    double minVal = layer[0].min();
    double minTmp;

    for(long i = 1; i < layerCount; i++)
    {
    minTmp = layer[i].min();
    minVal = (minTmp < minVal) ? minTmp : minVal;
    }
    return minVal;
}

double max() const
{
    double maxVal = layer[0].max();
    double maxTmp;
    for(long i = 1; i < layerCount; i++)
    {
    maxTmp = layer[i].max();
    maxVal = (maxTmp > maxVal) ? maxTmp : maxVal;
    }
    return maxVal;
}


void axpy(double alpha, const VLayeredGridFun3d& x)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].axpy(alpha,(x.layer[i]));
    }
} 

void axpby(double alpha,const VLayeredGridFun3d& x,double beta)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].axpby(alpha,(x.layer[i]),beta);
    }
}

void copy(const VLayeredGridFun3d& x)
{
	this->operator=(x);
}

void scal(double alpha)
{
    this->operator*=(alpha);
}

//
// The square root of the sum of the squares of all the grid node
// values. This value is not scaled by the mesh size. The interior
// interface values are averages of the interface values associated
// with the layers on either side of the interface.
//

double nrm2() const
{
	return std::sqrt(std::abs(this->dot(*this)));
}


void enforceXYperiodicity()
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].enforceXYperiodicity();
    }
}

void enforceXYperiodicityAveraging()
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].enforceXYperiodicityAveraging();
    }
}


void setBoundaryValues(double val)
{
	// Sides

	for(long Lindex = 0; Lindex < layerCount; Lindex++)
	{
		for(long k = 0; k <= zPanels[Lindex]; k++)
		{
			for(long i = 0; i <= xPanels; i++)
			{
			layer[Lindex](i,0,k)       = val;
			layer[Lindex](i,yPanels,k) = val;
			}

			for(long j = 0; j <= yPanels; j++)
			{
			layer[Lindex](0,j,k)         = val;
			layer[Lindex](xPanels,j,k)   = val;
			}

		}
	}

	// Top and bottom


	for(long i = 0; i <= xPanels; i++)
	{
	for(long j = 0; j <= yPanels; j++)
	{
		layer[0](i,j,0)                                = val;
		layer[layerCount-1](i,j,zPanels[layerCount-1]) = val;
	}}
}

bool isNull() const
{
	if(layerCount == 0) {return true;}
	return false;
}
//
//  Data Members
//

    double                      xWidth;
    double                        xMin;
    double                        xMax;
    long                       xPanels;


    double                       yWidth;
    double                         yMin;
    double                         yMax;
    long                        yPanels;

    std::vector< GridFunction3d >      layer;
    long                            layerCount;
    std::vector<double> 			     zWidth;
    std::vector<double> 			     zBdrys;
    std::vector<long>                   zPanels;

};
}

#endif
