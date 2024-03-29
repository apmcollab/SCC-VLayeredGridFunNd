
//
//####################################################################
//                    VLayeredGridFun1d.h
//####################################################################
/**
   An aggregate linear array of SCC::GridFunnction1d objects that is used
   to represented a layered structure stacked "vertically" with the 
   z-coordinate being the vertical direction  <p>
   Indexing for the layers starts at 0.


   The coordinate for the SCC::GridFunction1d instances that comprise
   the layers is x, so the coordinate labels used when accessing the
   layers individually has to accommdodate this difference.

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


#ifndef SCC_VLAYERED_GRID_FUN1D_
#define SCC_VLAYERED_GRID_FUN1D_

#include <functional>
#include <iostream>
#include <vector>

#include "GridFunctionNd/SCC_GridFunction1d.h"

namespace SCC
{

class VLayeredGridFun1d
{

public:

VLayeredGridFun1d()
{
	initialize();
}

VLayeredGridFun1d(long layerCount, const std::vector<long>& zPanels, const std::vector<double>& zBdrys)
{
    initialize(layerCount,zPanels,zBdrys);
}

VLayeredGridFun1d(const VLayeredGridFun1d& M)
{
    if(M.layer.size() == 0){initialize(); return;}
    initialize(M);
}

virtual void initialize()
{
	layerCount = 0;
	layer.clear();
    zWidth.clear();
    zBdrys.clear();
    zPanels.clear();
}

virtual void initialize(const VLayeredGridFun1d& M)
{
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

virtual void initialize(long layerCount, const std::vector<long>& zPanels, const std::vector<double>& zBdrys)
{
    initialize();

    this->layerCount = layerCount;
    this->layer.resize(layerCount);
    this->zBdrys.resize(layerCount+1);
    this->zWidth.resize(layerCount);
    this->zPanels.resize(layerCount);

    double a; double b;

    for(long i = 0; i < layerCount+1; i++)
    {
    this->zBdrys[i] = zBdrys[i];
    }

    for(long i = 0; i < layerCount; i++)
    {
    a = zBdrys[i];
    b = zBdrys[i+1];
    layer[i].initialize(zPanels[i],a,b);
    this->zWidth[i] = b-a;
    this->zPanels[i] = zPanels[i];
    }
}

virtual ~VLayeredGridFun1d()
{
}

long getLayerCount() const
{
	return layerCount;
}

std::vector<double> getZbdrys() const
{
	return zBdrys;
}

double getZmin() const
{return layer[0].getXmin();}

double getZmax() const
{return layer[layerCount-1].getXmax();}


std::vector<long> getZpanels() const
{
	return zPanels;
}

long getZpanelCountSum() const
{
	long panelSum = 0;
	for(long k = 0; k < layerCount; k++)
	{
		panelSum += layer[k].getXpanelCount();
	}
	return panelSum;
}



// Two structures are equal if they have the same number of
// layers and each layer is the same width and has the same
// number of panels. It is not required that the structures
// have identical x-boundary points.
//
virtual bool isEqualStructure(const VLayeredGridFun1d& V) const
{

    if((layer.size() == 0)&&(V.layer.size() != 0))       return false;
    if(layerCount != V.layerCount)                       return false;

    for(long i = 0; i < layerCount; i++)
    {
    if(zPanels[i] != V.zPanels[i]){return false;}
    else if(std::abs(zWidth[i] - V.zWidth[i]) > 1.0e-012*std::abs(zWidth[i])){return false;}
    }

    return true;
}


void setBoundaryValues(double val)
{
	setAboundaryValue(val);
	setBboundaryValue(val);
}

inline double getAboundaryValue() const
{
    return layer[0].operator()(0);
}

inline double getBboundaryValue() const
{
    return layer[layerCount-1].operator()(layer[layerCount-1].xPanels);
}

inline void setAboundaryValue(double v)
{
    layer[0].operator()(0) = v;
}
    
inline void setBboundaryValue(double v)
{
    layer[layerCount-1].operator()(layer[layerCount-1].xPanels) = v;
} 

inline void zeroAboundaryValue()
{
    layer[0].operator()(0) = 0.0;
}

inline void zeroBboundaryValue()
{
    layer[layerCount-1].operator()(layer[layerCount-1].xPanels) = 0.0;
}

inline void addToAboundaryValue(double v)
{
    layer[0].operator()(0) += v;
}

inline void addToBboundaryValue(double v)
{
    layer[layerCount-1].operator()(layer[layerCount-1].xPanels) += v;
}

inline void scaleAboundaryValue(double s)
{
    layer[0].operator()(0) *= s;
}

inline void scaleBboundaryValue(double s)
{
   layer[layerCount-1].operator()(layer[layerCount-1].xPanels) *= s;
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

void specify(std::function<double(double)> F)
{

	for(long k = 0; k < layerCount; k++)
	{
	layer[k].specify(F);
	}
}


VLayeredGridFun1d  extractLayers(long begIndex, long endIndex) const
{

	long  lcount = (endIndex - begIndex) + 1; 
	std::vector<long>   zPn(lcount);
	std::vector<double> zBd(lcount+1);

   for(long k = begIndex; k <= endIndex; k++)
   {
   zPn[k-begIndex] = layer[k].getXpanelCount();
   }

   zBd[0] = layer[begIndex].getXmin();
   for(long k = begIndex; k <= endIndex; k++)
   {
   zBd[k-begIndex+1] = layer[k].getXmax();
   }

   VLayeredGridFun1d R(lcount,zPn,zBd);

   for(long k = begIndex; k <= endIndex; k++)
   {
   R.layer[k-begIndex] = layer[k];
   }

   return R;
}


void  insertLayers(long begIndex, long endIndex,
const VLayeredGridFun1d& V)
{
   for(long k = begIndex; k <= endIndex; k++)
   {
   layer[k] = V.layer[k-begIndex];
   }
}

double getXmin() const
{
    return layer[0].getXmin();
}

double getXmax() const
{
    return layer[layerCount-1].getXmax();
}


void operator=(const VLayeredGridFun1d& M)
{
	initialize(M);
}

void setLayerValues(long k,const GridFunction1d& L)
{
    layer[k] = (DoubleVector1d) L;
}

GridFunction1d  getLayer(long k)
{
    return layer[k];
}

GridFunction1d* getLayerPtr(long k)
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

VLayeredGridFun1d operator-(const VLayeredGridFun1d& M)
{
    VLayeredGridFun1d R(M);
    for(long i = 0; i < layerCount; i++)
    {
    R.layer[i] = layer[i]- M.layer[i];
    }
    return R;
}

void operator-=(const VLayeredGridFun1d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] -= M.layer[i];
    }
}

void operator+=(const VLayeredGridFun1d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] += M.layer[i];
    }
}

void operator*=(const VLayeredGridFun1d& M)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *= M.layer[i];
    }
}

VLayeredGridFun1d operator+(const VLayeredGridFun1d& M)
{
    VLayeredGridFun1d R(M);
    for(long i = 0; i < layerCount; i++)
    {
    R.layer[i] = layer[i] + M.layer[i];
    }
    return R;

}

void operator*=(const std::function<double(double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] *=F;
    }
}

void operator/=(const std::function<double(double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] /=F;
    }
}

void operator+=(const std::function<double(double)>& F)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i] +=F;
    }
}

void operator-=(const std::function<double(double)>& F)
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
double dot(const VLayeredGridFun1d& M) const
{
    return scaledDot(M);
}

double scaledDot(const VLayeredGridFun1d& M) const
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

double norm1() const
{
    double norm1val = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    norm1val += layer[i].norm1();
    }
    return norm1val;
}

double norm2() const
{
    double norm2val = 0.0;
    for(long i = 0; i < layerCount; i++)
    {
    norm2val += layer[i].norm2squared();
    }
    return std::sqrt(norm2val);
}


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


void axpy(double alpha, const VLayeredGridFun1d& x)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].axpy(alpha,(x.layer[i]));
    }
} 

void axpby(double alpha,const VLayeredGridFun1d& x,double beta)
{
    for(long i = 0; i < layerCount; i++)
    {
    layer[i].axpby(alpha,(x.layer[i]),beta);
    }
}

void copy(const VLayeredGridFun1d& x)
{
	this->operator=(x);
}

void scal(double alpha)
{
    this->operator*=(alpha);
}

double nrm2() const
{
	return std::sqrt(std::abs(this->dot(*this)));
}

bool isNull() const
{
	if(layerCount == 0) {return true;}
	return false;
}


//
//  Data Members
//
    std::vector< GridFunction1d >         layer;
    long                             layerCount;
    std::vector<double> 			     zWidth;
    std::vector<double> 			     zBdrys;
    std::vector<long>                    zPanels;

};
}

#endif
