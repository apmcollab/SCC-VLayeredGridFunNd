
//
//####################################################################
//                    VLayeredGridFun2d.h
//####################################################################
/**
   An aggregate linear array of SCC::GridFunnction2d objects that is used
   to represented a layered structure stacked "vertically" with the 
   z-coordinate being the vertical direction  <p>

   Indexing for the layers starts at 0.

   The coordinates for the SCC::GridFunction2d instances that comprise
   the layers are (x,y), so the coordinate labels used when accessing the
   layers individually has to accommodate this difference. Specifically,
   each layer consists of an SCC::GridFunction2d instance whose
   y-coordinate then becomes the z-coordinate of the layer in
   the layered structure.

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


#ifndef SCC_VLAYERED_GRID_FUN2D_
#define SCC_VLAYERED_GRID_FUN2D_

#include <functional>
#include <iostream>
#include <vector>

#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "SCC_VLayeredGridFun1d.h"

namespace SCC
{

class VLayeredGridFun2d
{

public:

VLayeredGridFun2d()
{
	initialize();
}

VLayeredGridFun2d(long xPanels, double xMin, double xMax,
				  long layerCount, const std::vector<long>& zPanels, const std::vector<double>& zBdrys)
{
    initialize(xPanels,xMin,xMax,layerCount,zPanels,zBdrys);
}

VLayeredGridFun2d(const VLayeredGridFun2d& M)
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

	layerCount = 0;
	layer.clear();
    zWidth.clear();
    zBdrys.clear();
    zPanels.clear();
}

virtual void initialize(const VLayeredGridFun2d& M)
{
	xWidth   = M.xWidth;
    xMin     = M.xMin;
    xMax     = M.xMax;
    xPanels  = M.xPanels;

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
			    long layerCount, const std::vector<long>& zPanels, const std::vector<double>& zBdrys)
{
    initialize();

	this->xWidth   = xMax-xMin;
    this->xMin     = xMin;
    this->xMax     = xMax;
    this->xPanels  = xPanels;

    this->layerCount = layerCount;
    this->layer.resize(layerCount);

    this->zBdrys  = zBdrys;
    this->zPanels = zPanels;

    this->zWidth.resize(layerCount);

    double c; double d;


    for(long i = 0; i < layerCount; i++)
    {
    c = zBdrys[i];
    d = zBdrys[i+1];
    layer[i].initialize(xPanels,xMin,xMax,zPanels[i],c,d);
    this->zWidth[i]  = d-c;
    }
}

virtual ~VLayeredGridFun2d()
{}




// Two structures are equal if they have the same number of
// layers and each layer is the same width and has the same
// number of panels. It is not required that the structures
// have identical x-boundary points.
//
virtual bool isEqualStructure(const VLayeredGridFun2d& V) const
{

    if((layer.size() == 0)&&(V.layer.size() != 0))       return false;
    if(layerCount != V.layerCount)                       return false;

    // transverse structure

    if(xPanels != V.xPanels){return false;}
    if(xWidth  != V.xWidth) {return false;}

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

void specify(std::function<double(double,double)> F)
{

	for(long k = 0; k < layerCount; k++)
	{
	layer[k].specify(F);
	}
}


VLayeredGridFun2d  extractLayers(long begIndex, long endIndex) const
{

	long  lcount = (endIndex - begIndex) + 1; 
	std::vector<long>   zPn(lcount);
	std::vector<double> zBd(lcount+1);

	for(long k = begIndex; k <= endIndex; k++)
	{
		zPn[k-begIndex] = layer[k].getYpanelCount();
	}

	zBd[0] = layer[begIndex].getYmin();
	for(long k = begIndex; k <= endIndex; k++)
	{
		zBd[k-begIndex+1] = layer[k].getYmax();
	}

   VLayeredGridFun2d R(xPanels,xMin,xMax,lcount,zPn,zBd);

   for(long k = begIndex; k <= endIndex; k++)
   {
	R.layer[k-begIndex] = layer[k];
   }

   return R;
}


void  insertLayers(long begIndex, long endIndex, const VLayeredGridFun2d& V)
{
   for(long k = begIndex; k <= endIndex; k++)
   {
   layer[k] = V.layer[k-begIndex];
   }
}

void  incrementLayers(long begIndex, long endIndex, const VLayeredGridFun2d& V)
{
   for(long k = begIndex; k <= endIndex; k++)
   {
   layer[k] += V.layer[k-begIndex];
   }
}



VLayeredGridFun1d getConstantXslice(long xIndex) const //( z function)
{
	VLayeredGridFun1d R(layerCount,zPanels,zBdrys);

	for(long p = 0; p < layerCount; p++)
	{
		for(long k = 0; k <= zPanels[p]; k++)
		{
				R.layer[p](k) = layer[p](xIndex,k);
		}
    }
    return R;
}

// For the slice operations in which a reference object is
// passed in, the routines are not guaranteed to be thread
// safe unless the input argument is not shared among
// threads.
//

void getConstantXslice(long xIndex,VLayeredGridFun1d& R) const //( z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long k = 0; k <= zPanels[p]; k++)
		{
			R.layer[p](k) = layer[p](xIndex,k);
		}
    }
}

void setConstantXslice(long xIndex,const VLayeredGridFun1d& R)   //( z function)
{
	for(long p = 0; p < layerCount; p++)
	{
		for(long k = 0; k <= zPanels[p]; k++)
		{
		 layer[p](xIndex,k) = R.layer[p](k);
		}
    }
}

GridFunction1d getConstantZslice(long layerIndex, long zIndex) const //(x function)
{
    GridFunction1d R(xPanels,xMin,xMax);
    for(long i = 0; i <= xPanels; i++)
    {
    R.Values(i) = layer[layerIndex](i,zIndex);
    }
    return R;
}

void getConstantZslice(long layerIndex, long zIndex,GridFunction1d& R) const //(x function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    R.Values(i) = layer[layerIndex](i,zIndex);
    }
}

void setConstantZslice(long layerIndex, long zIndex, const GridFunction1d& R)  //(x function)
{
    for(long i = 0; i <= xPanels; i++)
    {
    layer[layerIndex](i,zIndex) = R.Values(i);
    }
}


void createProductFunction(const GridFunction1d& funX, const VLayeredGridFun1d& funZ)
{
    initialize(funX.xPanels, funX.xMin, funX.xMax,
               funZ.layerCount, funZ.zPanels, funZ.zBdrys);

	long i; long j; long p;

	double fX; double fZ;

	for(p = 0; p < layerCount; p++)
	{

	for(i = 0; i <= xPanels; i++)
	{
	fX = funX(i);
	for(j = 0; j <= zPanels[p]; j++)
	{
	fZ = funZ.layer[p](j);
	layer[p](i,j) = fX*fZ;
	}}

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

double getXmin()      const {return xMin;}
double getXmax()      const {return xMax;}
double getHx()        const {return xWidth/(double)xPanels;}
long getXpanelCount() const {return xPanels;}


double getZmin() const
{return layer[0].getYmin();}

double getZmax() const
{return layer[layerCount-1].getYmax();}


long getZpanelCountSum() const
{
	long panelSum = 0;
	for(long k = 0; k < layerCount; k++)
	{
		panelSum += layer[k].getYpanelCount();
	}
	return panelSum;
}

void operator=(const VLayeredGridFun2d& M)
{
	initialize(M);
}

void setLayerValues(long k,const SCC::GridFunction2d& L)
{
    layer[k] = (SCC::DoubleVector2d)L;
}

SCC::GridFunction2d  getLayer(long k)
{
    return layer[k];
}

SCC::GridFunction2d* getLayerPtr(long k)
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

VLayeredGridFun2d operator-(const VLayeredGridFun2d& M)
{
    VLayeredGridFun2d R(M);
    for(long i = 0; i < layerCount; i++)
    {
    R.layer[i] = layer[i]- M.layer[i];
    }
    return R;
}

void operator-=(const VLayeredGridFun2d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] -= M.layer[i];
    }
}

void operator+=(const VLayeredGridFun2d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] += M.layer[i];
    }
}

VLayeredGridFun2d operator+(const VLayeredGridFun2d& M)
{
    VLayeredGridFun2d R(M);
    for(long i = 0; i < layerCount; i++)
    {
    R.layer[i] = layer[i] + M.layer[i];
    }
    return R;

}

void operator*=(const VLayeredGridFun2d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *= M.layer[i];
    }
}


void operator*=(const std::function<double(double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *=F;
    }
}


void operator/=(const std::function<double(double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] /=F;
    }
}

void operator+=(const std::function<double(double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] +=F;
    }
}

void operator-=(const std::function<double(double,double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] -=F;
    }
}

void enforceXperiodicity()
{
    for(long k = 0; k < layerCount; k++)
    {
        layer[k].enforceXperiodicity();
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
double dot(const VLayeredGridFun2d& M) const
{
    return scaledDot(M);
}



// Dot product in which the integral is evaluated using the
// Trapezoidal method.

double scaledDot(const VLayeredGridFun2d& M) const
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
    return std::sqrt(norm2val);
}

// The required integral is approximated using the Trapezoidal method

double  integral() const
{
    double intVal = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    intVal += layer[i].integralTrapezoidal();
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


void axpy(double alpha, const VLayeredGridFun2d& x)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].axpy(alpha,(x.layer[i]));
    }
} 

void axpby(double alpha,const VLayeredGridFun2d& x,double beta)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].axpby(alpha,(x.layer[i]),beta);
    }
}

void copy(const VLayeredGridFun2d& x)
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

void setBoundaryValues(double val)
{
	// Sides

	for(long k = 0; k < layerCount; k++)
	{
		for(long j = 0; j <= zPanels[k]; j++)
		{
			layer[k](0,j)       = val;
			layer[k](xPanels,j) = val;
		}
	}

	// Top and bottom

	for(long i = 0; i <= xPanels; i++)
	{
		layer[0](i,0)                                = val;
		layer[layerCount-1](i,zPanels[layerCount-1]) = val;
	}
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


    std::vector< GridFunction2d >         layer;
    long                             layerCount;
    std::vector<double> 			     zWidth;
    std::vector<double> 			     zBdrys;
    std::vector<long>                   zPanels;

};
}

#endif
