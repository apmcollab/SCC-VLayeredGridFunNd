
//
//####################################################################
//                    VLayeredGridFunCnvrsn1d.h
//####################################################################
/**
   A class whose member functions are used to convert between a 
   layered grid function instance and a uniform grid. 

*/
//####################################################################
// Chris Anderson                                    Sept. 11, 2018
//####################################################################
//

/*
######################################################################
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
###################################################################
*/
#include <string>
#include <exception>

#include "VLayeredGridFunNd/SCC_VLayeredGridFun1d.h"
#include "VLayeredGridFunNd/SCC_VLayeredGridFun2d.h"
#include "VLayeredGridFunNd/SCC_VLayeredGridFun3d.h"

#include "GridFunctionNd/SCC_GridFunction1d.h"
#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "GridFunctionNd/SCC_GridFunction3d.h"

#ifndef SCC_VLAYERED_GRID_CNVRSN_ND_
#define SCC_VLAYERED_GRID_CNVRSN_ND_

namespace SCC
{

class VLayeredGridFunCnvrsnNd
{
public: 


bool identicalHz(const SCC::VLayeredGridFun1d& V)
{
    double hz = V.layer[0].getHx();
    for(long i = 0; i < V.getLayerCount(); i++)
    {
    	if(std::abs(hz - V.layer[i].getHx()) > 1.0e-14*hz){return false;}
    }
    return true;
}

bool identicalHz(const SCC::VLayeredGridFun2d& V)
{
    double hz = V.layer[0].getHy();
    for(long i = 0; i < V.getLayerCount(); i++)
    {
    	if(std::abs(hz - V.layer[i].getHy()) > 1.0e-14*hz){return false;}
    }
    return true;
}

bool identicalHz(const SCC::VLayeredGridFun3d& V)
{
    double hz = V.layer[0].getHz();
    for(long i = 0; i < V.getLayerCount(); i++)
    {
    	if(std::abs(hz - V.layer[i].getHz()) > 1.0e-14*hz){return false;}
    }
    return true;
}

//
// These routines check that the SCC::GridFunctionXd and SCC::VLayeredGridFunXd arguments have
// a consistent structure to support either extracting or inserting data from one to the other.
//
// This routine only checks that the values on the associated grids map from one to another,
// and does not check for consistency with respect uniform translation, i.e. consistency
// with the domain boundary specification.
//
// Returns:
//
// true  :  The arguments have a consistent structure
// false :  otherwise
//
bool checkForConsistentGridStructure(const SCC::GridFunction1d& uniformG, const SCC::VLayeredGridFun1d& G)
{
    if(not identicalHz(G))
    {
    	throw std::runtime_error("\n Layered grid Z mesh widths unequal : conversion to uniform grid function unsupported \n");
    }

    std::string eMessage;

    bool checkVal = true;
    if(uniformG.getXpanelCount() != G.getZpanelCountSum())    checkVal = false;

    if(not checkVal)
    {
    	eMessage =  "\nXXX VLayeredGridFunCnvrsnNd  XXX\n";
    	eMessage += "Layered grid function conversion to/from uniform grid function failed. \n";
    	throw std::runtime_error(eMessage);
    }

    return checkVal;
}
bool checkForConsistentGridStructure(const SCC::GridFunction2d& uniformG, const SCC::VLayeredGridFun2d& G)
{
    if(not identicalHz(G))
    {
    	throw std::runtime_error("\n Layered grid Z mesh widths unequal : conversion to uniform grid function unsupported \n");
    }

    std::string eMessage;

    bool checkVal = true;
    if(uniformG.getXpanelCount() != G.getXpanelCount())    checkVal = false;
    if(uniformG.getYpanelCount() != G.getZpanelCountSum()) checkVal = false;


    if(not checkVal)
    {
    	eMessage =  "\nXXX VLayeredGridFunCnvrsnNd  XXX\n";
    	eMessage += "Layered grid function conversion to/from uniform grid function failed. \n";
    	throw std::runtime_error(eMessage);
    }

    return checkVal;
}

bool checkForConsistentGridStructure(const SCC::GridFunction3d& uniformG, const SCC::VLayeredGridFun3d& G)
{
    if(not identicalHz(G))
    {
    	throw std::runtime_error("\n Layered grid Z mesh widths unequal : conversion to uniform grid function unsupported \n");
    }

    std::string eMessage;

    bool checkVal = true;
    if(uniformG.getXpanelCount() != G.getXpanelCount())    checkVal = false;
    if(uniformG.getYpanelCount() != G.getYpanelCount())    checkVal = false;
    if(uniformG.getZpanelCount() != G.getZpanelCountSum()) checkVal = false;


    if(not checkVal)
    {
    	eMessage =  "\nXXX VLayeredGridFunCnvrsnNd  XXX\n";
    	eMessage += "Layered grid function conversion to/from uniform grid function failed. \n";
    	throw std::runtime_error(eMessage);
    }

    return checkVal;
}

SCC::GridFunction1d createUniformGridFunction(const SCC::VLayeredGridFun1d& V)
{
    if(not identicalHz(V))
    {
    	throw std::runtime_error("\n Layered grid Z mesh widths unequal : conversion to uniform grid function unsupported \n");
    }

    long zPanels = V.getZpanelCountSum();
    double zMin  = V.getZmin();
    double zMax  = V.getZmax();

    SCC::GridFunction1d uniformGrid(zPanels,zMin,zMax);
    uniformGrid.setToValue(0.0);

    long zIndex = 0;
    long k;  long p;

   // First layer
    
    for(k = 0; k <=  V.zPanels[0]; k++)
    {
    	uniformGrid(zIndex) = V.layer[0](k);
    	zIndex += 1;
    }
    
    zIndex -= 1;
    
    // Remaining layers
    // Insert average interface values at interior interfaces
    
    for(p = 1; p < V.getLayerCount(); p++)
    {
    		uniformGrid(zIndex) = 0.5*(V.layer[p-1](V.zPanels[p-1]) + V.layer[p](0));
                zIndex += 1;
    
    		for(k = 1; k <= V.zPanels[p]; k++)
    		{
    			uniformGrid(zIndex) = V.layer[p](k);
    			zIndex += 1;
    		}
    
    		zIndex -= 1;
        }
    
    return uniformGrid;
    }
    
    
    SCC::GridFunction2d createUniformGridFunction(const SCC::VLayeredGridFun2d& V)
    {
    if(not identicalHz(V))
    {
    	throw std::runtime_error("\n Layered grid Z mesh widths unequal : conversion to uniform grid function unsupported \n");
    }
    
    long zPanels = V.getZpanelCountSum();
    double zMin  = V.getZmin();
    double zMax  = V.getZmax();
    
    SCC::GridFunction2d uniformGrid(V.xPanels,V.xMin,V.xMax,zPanels,zMin,zMax);
    
        long zIndex = 0;
        long zIndexStart;
    
        long i; long k; long p;
    
    // First layer
    
    for(i = 0; i <= V.xPanels; i++)
    {
    	zIndex = 0;
    	for(k = 0; k <=  V.zPanels[0]; k++)
    	{
    	uniformGrid(i,zIndex) = V.layer[0](i,k);
    	zIndex++;
    	}
    }
    
    zIndexStart = zIndex - 1;
    
    // Remaining layers
    // Insert average interface values at interior interfaces
    
    for(p = 1; p < V.getLayerCount(); p++)
    {
    	for(i = 0; i <= V.xPanels; i++)
    	{
    		zIndex = zIndexStart;
    		uniformGrid(i,zIndex) = 0.5*(V.layer[p-1](i,V.zPanels[p-1]) + V.layer[p](i,0));
    		zIndex++;
    
    		for(k = 1; k <= V.zPanels[p]; k++)
    		{
    			uniformGrid(i,zIndex) = V.layer[p](i,k);
    			zIndex++;
    		}
    	}
    	zIndexStart = zIndex - 1;
    }
    
    
    return uniformGrid;
}
    
    
SCC::GridFunction3d createUniformGridFunction(const SCC::VLayeredGridFun3d& V)
{
    if(not identicalHz(V))
    {
    	throw std::runtime_error("\n Layered grid Z mesh widths unequal : conversion to uniform grid function unsupported \n");
    }
    
    long zPanels = V.getZpanelCountSum();
    double zMin  = V.getZmin();
    double zMax  = V.getZmax();
    
	SCC::GridFunction3d uniformGrid(V.xPanels,V.xMin,V.xMax,V.yPanels,V.yMin,V.yMax,zPanels,zMin,zMax);
    
	long zIndex = 0;
	long zIndexStart;
    
	long i; long j; long k; long p;
    
    // First layer
    
    for(i = 0; i <= V.xPanels; i++)
    {
    for(j = 0; j <= V.yPanels; j++)
    {
    	zIndex = 0;
    	for(k = 0; k <=  V.zPanels[0]; k++)
    	{
    	uniformGrid(i,j,zIndex) = V.layer[0](i,j,k);
    	zIndex++;
    	}
    }}
    
    zIndexStart = zIndex - 1;
    
    
    // Remaining layers
    // Insert average interface values at interior interfaces
    
    for(p = 1; p < V.getLayerCount(); p++)
    {
    	for(i = 0; i <= V.xPanels; i++)
    	{
    	for(j = 0; j <= V.yPanels; j++)
    	{
    		zIndex = zIndexStart;
    		uniformGrid(i,j,zIndex) = 0.5*(V.layer[p-1](i,j,V.zPanels[p-1]) + V.layer[p](i,j,0));
    		zIndex++;
    
    		for(k = 1; k <= V.zPanels[p]; k++)
    		{
    			uniformGrid(i,j,zIndex) = V.layer[p](i,j,k);
    			zIndex++;
    		}
    	}}
    	zIndexStart = zIndex - 1;
    }
    
    return uniformGrid;
}


void insertUniformGridFunction3d(const SCC::GridFunction3d& uniformGrid, SCC::VLayeredGridFun3d& V)
{
    if(not checkForConsistentGridStructure(uniformGrid,V))
    {
    	throw std::runtime_error("\n Insertion of uniform grid function into layered grid function failed \n");
    }


	long zIndex = 0;
	long zIndexStart;

	long i; long j; long k; long p;

	// First layer

	for(i = 0; i <= V.xPanels; i++)
	{
	for(j = 0; j <= V.yPanels; j++)
	{
		zIndex = 0;
		for(k = 0; k <=  V.zPanels[0]; k++)
		{
        	V.layer[0](i,j,k) = uniformGrid(i,j,zIndex);
        	zIndex++;
		}
	}}

	zIndexStart = zIndex - 1;


	// Remaining layers
	//
    // We assume continuity so duplicate interface values at interior interfaces

	for(p = 1; p < V.getLayerCount(); p++)
	{
		for(i = 0; i <= V.xPanels; i++)
		{
	    for(j = 0; j <= V.yPanels; j++)
		{
	    	zIndex = zIndexStart;
	    	V.layer[p-1](i,j,V.zPanels[p-1]) = uniformGrid(i,j,zIndex);
	    	V.layer[p](i,j,0)                = uniformGrid(i,j,zIndex);

        	zIndex++;

        	for(k = 1; k <= V.zPanels[p]; k++)
        	{
        		V.layer[p](i,j,k) = uniformGrid(i,j,zIndex);
        		zIndex++;
        	}
        }}

		zIndexStart = zIndex - 1;
      }

}


};

}
#endif
    
 
