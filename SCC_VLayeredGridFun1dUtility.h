//
// VLayeredGridFun1dUtility.cpp
//
// Author: Chris Anderson  
// (C) UCLA 2018
//
// Aug. 20, 2018
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;


#include "GridFunctionNd/SCC_GridFunction1d.h"
#include "GridFunctionNd/SCC_GridFunction1dUtility.h"
#include "SCC_VLayeredGridFun1d.h"


#ifndef _SCC_VLayeredGridFun1dUtility_
#define _SCC_VLayeredGridFun1dUtility_

// MS compilers generate warnings if fopen is used instead of fopen_s (a Microsoft specific language
// extension, so this macro implements the appropriate equivalent to fopen that MS wants when
// MS compilers are being used. In both versions, the
// macro returns a non-zero value if the open fails (e.g. a non-zero error code).
//
#ifndef _MSC_VER
#ifndef OPENFILE
#define OPENFILE(dataFile,filename,mode) ((dataFile = fopen(filename,  mode)) == NULL)
#endif
#else
#ifndef OPENFILE
#define OPENFILE(dataFile,fileName,mode) ((fopen_s(&dataFile,fileName, mode)) != 0)
#pragma warning(push)
#pragma warning(disable: 4996)
#endif
#endif


namespace SCC
{

class VLayeredGridFun1dUtility
{

public:

void outputToGNUplot(const VLayeredGridFun1d& gridFun,const string& fileName,
					 const string& formatString = "%20.15e")
{
    ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << " \n";

    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    double   xp;
    double   hx;
    long mPanel;
    double xMin;

    for(long k = 0; k < gridFun.getLayerCount(); k++)
    {
    hx     = gridFun.layer[k].getHx();
    mPanel = gridFun.layer[k].getXpanelCount();
    xMin   = gridFun.layer[k].getXmin();

    if(k == 0){fprintf(dataFile,"\n");}

    for(long i = 0;  i <= mPanel; i++)
    {
    xp = xMin + i*hx;
    fprintf(dataFile,(s.str()).c_str(),xp, gridFun.layer[k](i));
    }
    }
    
    fclose(dataFile);
}


void appendToGNUplot(const VLayeredGridFun1d& gridFun,const string& fileName,
					 const string& formatString = "%20.15e")
{
    ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << " \n"; 


    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "a+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    double   xp;
    double   hx;
    long mPanel;
    double xMin;

    for(long k = 0; k < gridFun.getLayerCount(); k++)
    {
    hx     = gridFun.layer[k].getHx();
    mPanel = gridFun.layer[k].getXpanelCount();
    xMin   = gridFun.layer[k].getXmin();

    for(long i = 0;  i <= mPanel; i++)
    {
    xp = xMin + i*hx;
    fprintf(dataFile,(s.str()).c_str(),xp, gridFun.layer[k](i));
    }
    }

    fclose(dataFile);
} 



//###################################################
//   VLayeredGridFun1d Input/Output data structure
//###################################################
//
// layerCount
// # of z panels in layer 0
// # of z panels in layer 1
//   ***
// # of z panels in layerCount-1
// zBdry[0]
// zBdry[1]
//   *
//   *
// zBdry[layerCount]
//
// Grid values for layer 0
// Grid values for layer 1
//
// Grid values for layer layerCount-1
//
void outputToDataFile(const SCC::VLayeredGridFun1d& gF, const string& fileName,
const string& formatString = "%20.15e")
{
//
//  Create format string
//
    ostringstream s;
    s.str("");
    s << formatString << " ";
//
//  Open and then write to a file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long layerCount       = gF.getLayerCount();

    vector<long> zPanels  = gF.getZpanels();
    vector<double> zBdrys = gF.getZbdrys();


	fprintf(dataFile,"%ld \n", layerCount);
	for(long i = 0; i < layerCount; i++)
	{
	fprintf(dataFile,"%ld \n", zPanels[i]);
	}

	for(long i = 0; i <= layerCount; i++)
	{
	fprintf(dataFile,"%20.15e \n", zBdrys[i]);
	}

	for(long n = 0; n < layerCount; n++)
	{
		for(long k = 0; k <= zPanels[n]; k++)
		{
			fprintf(dataFile, s.str().c_str(),gF.layer[n](k));
		}
    }


    fclose(dataFile);
}


// Extracts the :VLayeredGridFun2d from the dataFile. The optional fileName argument
// provides the name of the originating file and is used to create a meaningful output
// message if the read fails, otherwise it is ignored.



void inputFromDataFile(SCC::VLayeredGridFun1d& gF, FILE* dataFile, string fileName = "")
{
    long layerCount;

	vector<long>  zPanels;
    vector<double> zBdrys;

    int rValue = 0;

	rValue = fscanf(dataFile,"%ld", &layerCount) != 1 ? 1 : rValue;

	zPanels.resize(layerCount);
	zBdrys.resize(layerCount+1);

	for(long i = 0; i < layerCount; i++)
	{
		rValue = fscanf(dataFile,"%ld", &zPanels[i]) != 1 ? 1 : rValue;
	}

	for(long i = 0; i <= layerCount; i++)
	{
		rValue =fscanf(dataFile,"%lf",&zBdrys[i]) != 1 ? 1 : rValue;
	}


    gF.initialize(layerCount,zPanels,zBdrys);

    for(long n = 0; n < layerCount; n++)
    {
    	for(long k = 0; k <= zPanels[n]; k++)
    	{
    		rValue = fscanf(dataFile,"%lf",&gF.layer[n](k)) != 1 ? 1 : rValue;
    	}
    }


    if(rValue == 1)
    {
    throw std::runtime_error("\nVLayeredGridFun1d could not be initialized from file " + fileName + " \n");
    }

}

void inputFromDataFile(SCC::VLayeredGridFun1d& gF, const string& fileName)
{
//
//  Open input file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "r" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    inputFromDataFile(gF, dataFile,fileName);

	fclose(dataFile);
}



//
// Using std::int64 to insure that the input and output of integers always
// uses 64 bit integer representation.
//
void outputToBinaryDataFile(const SCC::VLayeredGridFun1d& gF, FILE* dataFile)
{
    long dataSize;


    long layerCount                  = gF.getLayerCount();
    std::int64_t layerCount64        = gF.getLayerCount();

    vector<long>         zPanels     = gF.getZpanels();

    vector<std::int64_t> zPanels64(layerCount);

    for(long i = 0; i < layerCount; i++)
    {
    	zPanels64[i] = (std::int64_t)zPanels[i];
    }


    vector<double> zBdrys = gF.getZbdrys();

	//
	//  Write out the grid structure information. Using std:int64
	//  for integer values to avoid problems with machines with
	//  alternate storage sizes for int's and long's
	//


	fwrite(&layerCount64,  sizeof(std::int64_t), 1, dataFile);
	fwrite(&zPanels64[0],  sizeof(std::int64_t), layerCount, dataFile);
	fwrite(&zBdrys[0],  sizeof(double), layerCount+1, dataFile);

//  Write out the function values

	for(long n = 0; n < layerCount; n++)
	{
		dataSize = zPanels[n] + 1;
		fwrite(gF.layer[n].getDataPointer(),  sizeof(double), dataSize, dataFile);
	}
}

void outputToBinaryDataFile(const SCC::VLayeredGridFun1d& gF, const string& fileName)
{
//
//  Open and then write to a file
//
    FILE* dataFile;

    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    outputToBinaryDataFile(gF,dataFile);
    fclose(dataFile);
}

void inputFromBinaryDataFile(SCC::VLayeredGridFun1d& gF, FILE* dataFile, string fileName = "")
{
	size_t rValue;
    long dataSize;

    long            layerCount;
    vector<long>       zPanels;
    vector<double>      zBdrys;

	std::int64_t layerCount64;

	vector<std::int64_t> zPanels64;

    rValue = 0;

	rValue = fread(&layerCount64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;

	layerCount = (long)layerCount64;

	zPanels64.resize(layerCount);

	rValue = fread(&zPanels64[0],  sizeof(std::int64_t), layerCount, dataFile) != (uint) layerCount ? 0 : rValue;

	zPanels.resize(layerCount);

	for(long i = 0; i < layerCount; i++)
	{
		zPanels[i] = (long)zPanels64[i];
	}

	zBdrys.resize(layerCount+1);
    rValue = fread(&zBdrys[0], sizeof(double), layerCount+1, dataFile) != (uint)(layerCount+1) ? 1 : rValue;

	gF.initialize(layerCount,zPanels,zBdrys);

	for(long n = 0; n < layerCount; n++)
	{
	dataSize = (zPanels[n] + 1);
	rValue = fread(gF.layer[n].getDataPointer(),  sizeof(double), dataSize, dataFile) != (uint) dataSize ? 1 : rValue;
	}


    if(rValue == 1)
    {
    throw std::runtime_error("\nVLayeredGridFun1d could not be initialized from file " + fileName + " \n");
    }
}

void inputFromBinaryDataFile(SCC::VLayeredGridFun1d& gF, const string& fileName)
{
	//
	//  Open input file (remember to use the b mode to specify binary!!!!)
	//
	FILE* dataFile = 0;

	if(OPENFILE(dataFile,fileName.c_str(), "rb" ))
    {
    throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

	string dataFileName = fileName;

	inputFromBinaryDataFile(gF,dataFile,dataFileName);

    if(ferror(dataFile))
	{
    throw std::runtime_error("\nVLayeredGridFun1d could not be initialized from " + fileName + " \n");
	}

    fclose(dataFile);
}


//
// nVal = number of function values <= 5
// h    = grid size array of size nVal-1
// f    = function value array of size nVal
//
// returns
//
// dfA  = derivative at left  edge
//
void getOneSidedDerivativeApproximation(long nVal, double* h, double* f, double& dfA)
{
    if(nVal > 5) nVal = 5;

    double   x[5];
    double d1f[4];
    double d2f[3];
    double d3f[2];
    double d4f[1];

// x[0]        x[1]         x[2]         x[3]         x[4]
// |--- h[0] ---|--- h[1] ---|--- h[2] ---|--- h[3] ---|
// f[0]        f[1]         f[2]         f[3]         f[4]
//
// Divided difference construction
//
// x[0]  f0
// x[1]  f1 df0
// x[2]  f2 df1 ddf0
// x[3]  f3 df2 ddf1 dddf0
// x[4]  f4 df3 ddf2 dddf1 ddddf0
//
//  f = d0 + (x-x[0])*d1f0 + (x-x[0])(x-x[1])d2f0 + (x-x[0])(x-x[1])(x-x[2])dd3f0 + (x-x[0])(x-x[1])(x-x[2])(x-x[3])d4f0

    x[0] = 0.0;
    for(long i = 1; i < nVal; i++)
    {
    x[i] = x[i-1] + h[i-1];
    }

    for(long i = 0; i < nVal-1; i++)
    {
    d1f[i] = (f[i+1] - f[i])/(x[i+1] - x[i]);
    }

    for(long i = 0; i < nVal-2; i++)
    {
    d2f[i] = (d1f[i+1] - d1f[i])/(x[i+2] - x[i]);
    }

    for(long i = 0; i < nVal-3; i++)
    {
    d3f[i] = (d2f[i+1] - d2f[i])/(x[i+3] - x[i]);
    }

    /*
    for(long i = 0; i < nVal-4; i++)
    {
    d4f[i] = (d3f[i+1] - d3f[i])/(x[i+4] - x[i]);
    }
    */
    if(nVal == 5)
    {
    d4f[0] = (d3f[1] - d3f[0])/(x[4] - x[0]);
    }

    dfA = d1f[0];
    if(nVal >= 3)
    {
    dfA += (x[0]-x[1])*d2f[0];
    if(nVal >= 4)
    {
    dfA += (x[0]-x[1])*(x[0]-x[2])*d3f[0];
    if(nVal >= 5)
    {
    dfA += (x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3])*d4f[0];
    }
    }}
}
//
// This routine returns the coefficients of up to a fourth order one sided derivative
// approximation based on the Lagrange interpolating polynomial associated with the
// following grid structure
//
// x[0]        x[1]         x[2]         x[3]         x[4]
// |--- h[0] ---|--- h[1] ---|--- h[2] ---|--- h[3] ---|
// f[0]        f[1]         f[2]         f[3]         f[4]
//
//
//  This routine returns two arrays of size nVal so that the positive and negative sums can
//  be accumulated separately to reduce the effects of roundoff when the mesh size is small.
//
//  Code snipet to evaluate the derivative
//
//  getOneSidedDerivativeCoeffients(nVal, h, dAcoeffPos,dAcoeffNeg);
//
// 	dfA   = 0.0;
// 	dfPos = 0.0;
// 	dfNeg = 0.0;
// 	for(long i = 0; i < nVal; i++)
// 	{
// 	dfPos += dAcoeffPos[i]*f[i];
// 	dfNeg += dAcoeffNeg[i]*f[i];
// 	}
// 	dfA    = dfPos + dfNeg;
//

void getOneSidedDerivativeCoeffients(long nVal, double* h, double* dfAcoeffPos, double*dfAcoeffNeg)
{
	if(nVal > 5) nVal = 5;

	double f[5];
	double dfA;

	for(long k = 0; k < nVal; k++)
	{
	for(long j = 0; j < nVal; j++) {f[j]=0.0;}
	f[k] = 1.0;
	dfA = 0.0;
	getOneSidedDerivativeApproximation(nVal,h,f,dfA);
	if(dfA > 0)
	{
	 dfAcoeffPos[k] = dfA;
	 dfAcoeffNeg[k] = 0.0;
	}
    {
     dfAcoeffPos[k] = 0.0;
	 dfAcoeffNeg[k] = dfA;
    }
	}
}
//
// To obtain the one sided derivative on the left side of a layer use interfaceIndex = layerIndex
//
// To obtain the one sided derivative on the right side of a layer use interfaceIndex = layerIndex+1
//
// layer indices    start at 0
// layer interface  indices start at 0
// layer interfaces include first and last interval baoundaries
//
//                                  interface index k
//                                        |
//                            layer k-1   |  Layer k   |
//                                        |
//                                  ----> |  --->
//                   layerIndex     = k-1   layer Index =    k
//                   interfaceIndex = k     interfaceIndex = k
//
//
double getInterfaceDerivative(SCC::VLayeredGridFun1d& F,
long layerIndex, long interfaceIndex)
{
	double   dF;

    long   nVal;

    double h[4];
    double f[5];
    double hx;


    double dAcoeffPos[5];
	double dAcoeffNeg[5];
	double dfPos;
	double dfNeg;


	long layerPanelCount = F.layer[layerIndex].getXpanelCount();
    nVal = layerPanelCount + 1;

    if(nVal > 5) nVal = 5;

	hx = F.layer[layerIndex].getHx();

	for(long i = 0; i < nVal; i++)
	{
		h[i] = hx;
	}

	getOneSidedDerivativeCoeffients(nVal, h, dAcoeffPos,dAcoeffNeg);
    //
    // One sided derivative on the left
    //
	if(layerIndex == interfaceIndex)
	{
    	for(long j = 0; j < nVal; j++)
    	{
    		f[j] = F.layer[layerIndex].Values(j);
    	}

    	dfPos = 0.0; dfNeg = 0.0;
    	for(long k = 0; k < nVal; k++)
    	{
    		dfPos += dAcoeffPos[k]*f[k];
    		dfNeg += dAcoeffNeg[k]*f[k];
    	}

    	dF = dfPos + dfNeg;
	}
	else
	// Assuming (layerIndex == interfaceIndex-1)
	//
	// One sided derivative on the right
	//
	{

    	for(long j = 0; j < nVal; j++)
    	{
    	    f[j] = F.layer[layerIndex].Values(layerPanelCount - j);
    	}

    	dfPos = 0.0; dfNeg = 0.0;
    	for(long k = 0; k < nVal; k++)
    	{
    		dfPos += dAcoeffPos[k]*f[k];
    		dfNeg += dAcoeffNeg[k]*f[k];
    	}

    	dF = -(dfPos + dfNeg); // (-1) since reflecting difference formula
	}

	return dF;

}





};






}
#endif
