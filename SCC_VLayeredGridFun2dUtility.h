//
// VLayeredGridFun2dUtility.cpp
//
// Author: Chris Anderson
// (C) UCLA 2018
//
// Sept. 1, 2018
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


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>

#include "SCC_VLayeredGridFun2d.h"

#ifndef  SCC_VLAYERED_GRID_FUN2D_UTILITY_
#define  SCC_VLAYERED_GRID_FUN2D_UTILITY_

namespace SCC
{

class VLayeredGridFun2dUtility
{

public:


void outputToGNUplot(const VLayeredGridFun2d& gF, const std::string& fileName, const std::string& formatString = "%20.15e")
{
    std::ostringstream s;
    s.str("");
    s << formatString << "  " << formatString << "  " << formatString << " \n";
//
//  Open and then write to a file
//
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    long xPanels = gF.getXpanelCount();

    long layerCount          = gF.getLayerCount();
    std::vector<long>   zPanels   = gF.getZpanels();
    std::vector<double> zBdrys    = gF.getZbdrys();

    double xMin  = gF.getXmin();
    double xMax  = gF.getXmax();

    double xPos = 0; double zPos;

    double hx   = (xMax-xMin)/(double)(xPanels);
    double hz;

    for(long i = 0; i <= xPanels; i++)
    {
    	xPos = xMin + i*hx;
    	for(long n = 0; n < layerCount; n++)
    	{
    		hz = (zBdrys[n+1] - zBdrys[n])/(double)zPanels[n];
    		for(long k = 0; k <= gF.layer[n].getYpanelCount(); k++)
    		{
    			zPos = zBdrys[n] + k*hz;
    			fprintf(dataFile,(s.str()).c_str(),xPos,zPos,gF.layer[n](i,k));
    		}
    }

    fprintf(dataFile,"\n");
    }


    fclose(dataFile);
}



// Outputs the 2D data as 3D data with Y being a dummy coordinate.
//
// The data values output on the interfaces between layers are the average of the
// values of the interface values of the layers on either side of the interface.
//
// dataLabel must be less than 256 characters in length
//

void outputDataToVTKfile(const SCC::VLayeredGridFun2d& gridFun, const std::string& fileName, const std::string& dataLabel)
{
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }
    //
    // Output 3D data for a 2D function by using 2 points in the Y direction
    //
    double a  = gridFun.getXmin();

	double e  = gridFun.getZmin();

    double hx = gridFun.getHx();

    double c  =  gridFun.getXmin();
    double hy = (gridFun.getXmax()-c);

    double hz;

    long mPt = gridFun.getXpanelCount()    + 1;
    long nPt = 2;
    long pPt = gridFun.getZpanelCountSum() + 1;

	std::vector<long>   zPanels   = gridFun.getZpanels();
	std::vector<double> zBdrys    = gridFun.getZbdrys();

	long nLayers             = gridFun.getLayerCount();

    long dataCount           = mPt*nPt*pPt;


    long i; long j; long k;  long  n;
    double xPos; double yPos; double zPos;

	fprintf(dataFile, "# vtk DataFile Version 2.0\n");
    fprintf(dataFile, "%s \n",dataLabel.c_str());
    fprintf(dataFile, "ASCII\n");

    fprintf(dataFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(dataFile, "DIMENSIONS %ld %ld %ld \n",mPt,nPt,pPt);
    fprintf(dataFile, "X_COORDINATES %ld float \n",mPt);
    for(i = 0; i < mPt; i++)
    {
    xPos = i*hx + a;
    fprintf(dataFile, "%20.15e ",xPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Y_COORDINATES %ld float \n",nPt);
    for(j = 0; j < nPt; j++)
    {
    yPos = j*hy + c;
    fprintf(dataFile, "%20.15e ",yPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Z_COORDINATES %ld float \n",pPt);

    for(n = 0; n < nLayers; n++)
    {
    	e  = zBdrys[n];
    	hz = gridFun.layer[n].getHy();
    	for(k = 0; k < zPanels[n]; k++)
    	{
    		zPos = k*hz + e;
    		fprintf(dataFile, "%20.15e ",zPos);
    	}
    }
    zPos = zBdrys[nLayers];
    fprintf(dataFile, "%20.15e ",zPos);

    fprintf(dataFile, "\n");

    fprintf(dataFile, "POINT_DATA %ld\n",dataCount);
    fprintf(dataFile, "SCALARS %s float\n",dataLabel.c_str());
    fprintf(dataFile, "LOOKUP_TABLE default\n");
    double avgVal;

    // The complication in the output is the fact that the values output along the
    // layer interface are the averages of values associated with the layers
    // on either side of the interface.

    for(n = 0; n < nLayers; n++)
    {

        if(n == 0)
        {
        	for(k = 0; k < zPanels[n]; k++)
        	{
        		for(j = 0; j < nPt; j++)
        		{
        			for(i = 0; i < mPt; i++)
        			{
        				fprintf(dataFile, "%20.15e ",gridFun.layer[n](i,k));
        			}
        			fprintf(dataFile, "\n");
        		}
        	}
        }
        else
        {
        	k = 0;
        	for(j = 0; j < nPt; j++)
            {
        	for(i = 0; i < mPt; i++)
            {
        		avgVal = 0.5*(gridFun.layer[n-1](i,zPanels[n-1])+gridFun.layer[n](i,0));
        		fprintf(dataFile, "%20.15e ",avgVal);
            }
        	   fprintf(dataFile, "\n");
            }


        	for(k = 1;  k < zPanels[n]; k++)
        	{
        		for(j = 0; j < nPt; j++)
        		{
        			for(i = 0; i < mPt; i++)
        			{
        				fprintf(dataFile, "%20.15e ",gridFun.layer[n](i,k));
        			}
        			fprintf(dataFile, "\n");
        		}
        	}
        }
    }

    for(j = 0; j < nPt; j++)
    {
    	for(i = 0; i < mPt; i++)
    	{
    		fprintf(dataFile, "%20.15e ",gridFun.layer[nLayers-1](i,zPanels[nLayers-1]));
        }
    	fprintf(dataFile, "\n");
    }

    fclose(dataFile);
}


// Outputs the 2D data as 3D data with Y being a dummy coordinate.
//
//
// The data values output on the interfaces between layers are the average of the
// values of the interface values of the layers on either side of the interface.
//
// This routine scales the coordinate values of the scaling coordinate specified (one of "x" or "z") so that
// the coordinate extent is the scaling value times the largest size of the domain in the other coordinate
// directions. This utility is added so that rectangular regions in which the specified coordinate is thin
// with respect to the others can be viewed with the thin region expanded.
//
// dataLabel must be less than 256 characters in length
//
void outputDataToVTKfile(const VLayeredGridFun2d& gridFun, const std::string& fileName, const std::string& dataLabel,
std::string scalingCoord, double scalingValue)
{
    FILE* dataFile = 0;
    if(OPENFILE(dataFile,fileName.c_str(), "w+" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    double a  = gridFun.getXmin();  double b  = gridFun.getXmax();
	double e  = gridFun.getZmin();  double f  = gridFun.getZmax();

    double hx = gridFun.getHx();

    // Dummy dimension

    double hy = hx;
    double c  = 0.0;

    double hz;

    long mPt = gridFun.getXpanelCount()    + 1;
    long nPt = 2;
    long pPt = gridFun.getZpanelCountSum() + 1;

	std::vector<long>   zPanels   = gridFun.getZpanels();
	std::vector<double> zBdrys    = gridFun.getZbdrys();

	long nLayers             = gridFun.getLayerCount();

    long dataCount           = mPt*nPt*pPt;

    double xScalingFactor = -1.0;
    double yScalingFactor =  1.0;
    double zScalingFactor = -1.0;

    if((scalingCoord == "x") || (scalingCoord == "X")  ){xScalingFactor = scalingValue;}
    if((scalingCoord == "z") || (scalingCoord == "Z")  ){zScalingFactor = scalingValue;}

    double transverseSizeMax;

    if(xScalingFactor > 0.0)
    {
    transverseSizeMax = f-c;
    xScalingFactor = (xScalingFactor*transverseSizeMax)/(b-a);
    }
    else {xScalingFactor = 1.0;}

    if(zScalingFactor > 0.0)
    {
    transverseSizeMax = b-a;
    zScalingFactor = (zScalingFactor*transverseSizeMax)/(f-e);
    }
    else {zScalingFactor = 1.0;}

    long i; long j; long k;  long  n;
    double xPos; double yPos; double zPos;

	fprintf(dataFile, "# vtk DataFile Version 2.0\n");
    fprintf(dataFile, "%s \n",dataLabel.c_str());
    fprintf(dataFile, "ASCII\n");

    fprintf(dataFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(dataFile, "DIMENSIONS %ld %ld %ld \n",mPt,nPt,pPt);
    fprintf(dataFile, "X_COORDINATES %ld float \n",mPt);
    for(i = 0; i < mPt; i++)
    {
    xPos = i*hx*xScalingFactor  + a*xScalingFactor ;
    fprintf(dataFile, "%20.15e ",xPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Y_COORDINATES %ld float \n",nPt);
    for(j = 0; j < nPt; j++)
    {
    yPos = j*hy*yScalingFactor  + c*yScalingFactor ;
    fprintf(dataFile, "%20.15e ",yPos);
    }
    fprintf(dataFile, "\n");

    fprintf(dataFile, "Z_COORDINATES %ld float \n",pPt);

    for(n = 0; n < nLayers; n++)
    {
    	e  = zBdrys[n];
    	hz = gridFun.layer[n].getHy();
    	for(i = 0; i < zPanels[n]; i++)
    	{
    		zPos = i*hz*zScalingFactor  + e*zScalingFactor ;
    		fprintf(dataFile, "%20.15e ",zPos);
    	}
    }
    zPos = zBdrys[nLayers]*zScalingFactor;
    fprintf(dataFile, "%20.15e ",zPos);

    fprintf(dataFile, "\n");

    fprintf(dataFile, "POINT_DATA %ld\n",dataCount);
    fprintf(dataFile, "SCALARS %s float\n",dataLabel.c_str());
    fprintf(dataFile, "LOOKUP_TABLE default\n");
    double avgVal;

    // The complication in the output is the fact that the values output along the
    // layer interface are the averages of values associated with the layers
    // on either side of the interface.

    for(n = 0; n < nLayers; n++)
    {

        if(n == 0)
        {
        	for(k = 0; k < zPanels[n]; k++)
        	{
        		for(j = 0; j < nPt; j++)
        		{
        			for(i = 0; i < mPt; i++)
        			{
        				fprintf(dataFile, "%20.15e ",gridFun.layer[n](i,k));
        			}
        			fprintf(dataFile, "\n");
        		}
        	}
        }
        else
        {
        	k = 0;
        	for(j = 0; j < nPt; j++)
            {
        	for(i = 0; i < mPt; i++)
            {
        		avgVal = 0.5*(gridFun.layer[n-1](i,zPanels[n-1])+gridFun.layer[n](i,0));
        		fprintf(dataFile, "%20.15e ",avgVal);
            }
        	   fprintf(dataFile, "\n");
            }


        	for(k = 1;  k < zPanels[n]; k++)
        	{
        		for(j = 0; j < nPt; j++)
        		{
        			for(i = 0; i < mPt; i++)
        			{
        				fprintf(dataFile, "%20.15e ",gridFun.layer[n](i,k));
        			}
        			fprintf(dataFile, "\n");
        		}
        	}
        }
    }

    for(j = 0; j < nPt; j++)
    {
    	for(i = 0; i < mPt; i++)
    	{
    		fprintf(dataFile, "%20.15e ",gridFun.layer[nLayers-1](i,zPanels[nLayers-1]));
        }
    	fprintf(dataFile, "\n");
    }

    fclose(dataFile);
}



//###################################################
//   VLayeredGridFun2d Input/Output data structure
//###################################################
//
// # of x panels
// layerCount
// # of z panels in layer 0
// # of z panels in layer 1
//   ***
// # of z points in layerCount-1
// xMin
// xMax
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
void outputToDataFile(const SCC::VLayeredGridFun2d& gF, const std::string& fileName,
const std::string& formatString = "%20.15e")
{
//
//  Create format std::string
//
    std::ostringstream s;
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


    long xPanels          = gF.getXpanelCount();
    long layerCount       = gF.getLayerCount();
    std::vector<long> zPanels  = gF.getZpanels();
    std::vector<double> zBdrys = gF.getZbdrys();


    double xMin  = gF.getXmin();
    double xMax  = gF.getXmax();

    fprintf(dataFile,"%ld \n", xPanels);
	fprintf(dataFile,"%ld \n", layerCount);
	for(long i = 0; i < layerCount; i++)
	{
	fprintf(dataFile,"%ld \n", zPanels[i]);
	}

    fprintf(dataFile,"%20.15e \n",xMin);
	fprintf(dataFile,"%20.15e \n",xMax);

	for(long i = 0; i <= layerCount; i++)
	{
	fprintf(dataFile,"%20.15e \n", zBdrys[i]);
	}


	for(long n = 0; n < layerCount; n++)
	{
    for(long i = 0; i <= xPanels; i++)
    {
    for(long k = 0; k <= zPanels[n]; k++)
    {
    fprintf(dataFile, s.str().c_str(),gF.layer[n](i,k));
    }
    fprintf(dataFile, "\n");
    }

	}

    fclose(dataFile);
}


// Extracts the :VLayeredGridFun2d from the dataFile. The optional fileName argument
// provides the name of the originating file and is used to create a meaningful output
// message if the read fails, otherwise it is ignored.



void inputFromDataFile(SCC::VLayeredGridFun2d& gF, FILE* dataFile, std::string fileName = "")
{
    long xPanels;

    double xMin; double xMax;

    long       layerCount;
	std::vector<long>  zPanels;
    std::vector<double> zBdrys;

    int rValue = 0;

    rValue = fscanf(dataFile,"%ld", &xPanels)    != 1 ? 1 : rValue;
	rValue = fscanf(dataFile,"%ld", &layerCount) != 1 ? 1 : rValue;

	zPanels.resize(layerCount);
	zBdrys.resize(layerCount+1);

	for(long i = 0; i < layerCount; i++)
	{
		rValue = fscanf(dataFile,"%ld", &zPanels[i]) != 1 ? 1 : rValue;
	}

    rValue =fscanf(dataFile,"%lf",&xMin) != 1 ? 1 : rValue;
	rValue =fscanf(dataFile,"%lf",&xMax) != 1 ? 1 : rValue;

	for(long i = 0; i <= layerCount; i++)
	{
		rValue =fscanf(dataFile,"%lf",&zBdrys[i]) != 1 ? 1 : rValue;
	}


    gF.initialize(xPanels,xMin,xMax,layerCount,zPanels,zBdrys);

    for(long n = 0; n < layerCount; n++)
    {
    for(long i = 0; i <= xPanels; i++)
    {
	for(long k = 0; k <= zPanels[n]; k++)
    {
    rValue = fscanf(dataFile,"%lf",&gF.layer[n](i,k)) != 1 ? 1 : rValue;
    }
    }

    }


    if(rValue == 1)
    {
    throw std::runtime_error("\nVLayeredGridFun2d could not be initialized from file " + fileName + " \n");
    }

}

void inputFromDataFile(SCC::VLayeredGridFun2d& gF, const std::string& fileName)
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
void outputToBinaryDataFile(const SCC::VLayeredGridFun2d& gF, FILE* dataFile)
{
    long dataSize;

    long xPanels                    = gF.getXpanelCount();
    std::int64_t xPanels64          = (std::int64_t)xPanels;

    long layerCount                 = gF.getLayerCount();
    std::int64_t layerCount64       = gF.getLayerCount();


    std::vector<long>         zPanels     = gF.getZpanels();
    std::vector<std::int64_t> zPanels64(layerCount);

    for(long i = 0; i < layerCount; i++)
    {
    	zPanels64[i] = (std::int64_t)zPanels[i];
    }


    double xMin  = gF.getXmin();
    double xMax  = gF.getXmax();

    std::vector<double> zBdrys = gF.getZbdrys();

	//
	//  Write out the grid structure information. Using std:int64
	//  for integer values to avoid problems with machines with
	//  alternate storage sizes for int's and long's
	//

    fwrite(&xPanels64,       sizeof(std::int64_t), 1, dataFile);
	fwrite(&layerCount64,    sizeof(std::int64_t), 1, dataFile);
	fwrite(&zPanels64[0],    sizeof(std::int64_t), layerCount, dataFile);

    fwrite(&xMin,       sizeof(double), 1, dataFile);
	fwrite(&xMax,       sizeof(double), 1, dataFile);
	fwrite(&zBdrys[0],  sizeof(double), layerCount+1, dataFile);

//  Write out the function values

	for(long n = 0; n < layerCount; n++)
	{
		dataSize = (xPanels+1)*(zPanels[n] + 1);
		fwrite(gF.layer[n].getDataPointer(),  sizeof(double), dataSize, dataFile);
	}
}

void outputToBinaryDataFile(const SCC::VLayeredGridFun2d& gF, const std::string& fileName)
{
//
//  Open and then write to a file
//
    FILE* dataFile;

    if(OPENFILE(dataFile,fileName.c_str(), "w+b" ))
    {
      throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

    outputToBinaryDataFile(gF,dataFile);
    fclose(dataFile);
}

void inputFromBinaryDataFile(SCC::VLayeredGridFun2d& gF, FILE* dataFile, std::string fileName = "")
{
	size_t rValue = 0;
    long dataSize;

    long               xPanels;
    double xMin;   double xMax;

    long            layerCount;
    std::vector<long>       zPanels;
    std::vector<double>      zBdrys;

	std::int64_t         xPanels64;
	std::int64_t      layerCount64;
	std::vector<std::int64_t> zPanels64;


	rValue = fread(&xPanels64,     sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&layerCount64,  sizeof(std::int64_t), 1, dataFile) != 1 ? 1 : rValue;

	layerCount = (long)layerCount64;
    xPanels    = (long)xPanels64;

	zPanels64.resize(layerCount);
    zPanels.resize(layerCount);

	rValue = fread(&zPanels64[0],  sizeof(std::int64_t), layerCount, dataFile) != (unsigned int) layerCount ? 0 : rValue;

	for(long i = 0; i < layerCount; i++)
	{
		zPanels[i] = (long)zPanels64[i];
	}

	zBdrys.resize(layerCount+1);

	rValue = fread(&xMin,      sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
	rValue = fread(&xMax,      sizeof(double), 1, dataFile) != 1 ? 1 : rValue;
    rValue = fread(&zBdrys[0], sizeof(double), layerCount+1, dataFile) != (unsigned int)(layerCount+1) ? 1 : rValue;

	gF.initialize(xPanels,xMin,xMax,layerCount,zPanels,zBdrys);

	for(long n = 0; n < layerCount; n++)
	{
	dataSize = (xPanels+1)*(zPanels[n] + 1);
	rValue = fread(gF.layer[n].getDataPointer(),  sizeof(double), dataSize, dataFile) != (unsigned int) dataSize ? 1 : rValue;
	}


    if(rValue == 1)
    {
    throw std::runtime_error("\nVLayeredGridFun2d could not be initialized from file " + fileName + " \n");
    }
}

void inputFromBinaryDataFile(SCC::VLayeredGridFun2d& gF, const std::string& fileName)
{
	//
	//  Open input file (remember to use the b mode to specify binary!!!!)
	//
	FILE* dataFile = 0;

	if(OPENFILE(dataFile,fileName.c_str(), "rb" ))
    {
    throw std::runtime_error("\nCannot open " + fileName + " \nFile not found.\n");
    }

	std::string dataFileName = fileName;

	inputFromBinaryDataFile(gF,dataFile,dataFileName);

    if(ferror(dataFile))
	{
    throw std::runtime_error("\nVLayeredGridFun2d could not be initialized from " + fileName + " \n");
	}

    fclose(dataFile);
}

};



}
#endif /*_VlayeredGridFun2dUtility_ */
