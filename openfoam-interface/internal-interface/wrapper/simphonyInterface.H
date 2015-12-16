/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2015 SimPhoNy -project
   \\/     M anipulation  |
-------------------------------------------------------------------------------


License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef simphonyInterface_H
#define simphonyInterface_H

#include "fvSchemes/fvSchemes.H"
#include "solution/solution.H"
#include "fvCFD.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "time/TimeMod.H"

#include<vector>
#include<string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// init dictionary for II wrapper
void foam_init(std::string name, std::string cD);
// init Foam mesh name (case name)  and rootPath for IO wrapper
void foam_init_IO(std::string name, std::string rootPath);
// read Foam mesh from case named name to wrapper
void foam_readMesh(std::string name);
// get mesh single point coordinates 
std::vector<double> foam_getPointCoordinates(std::string name, int label);
// get mesh  points coordinates
std::vector<double> foam_getPointCoordinates(std::string name);
// get data names associated (interpolated) to mesh points
std::vector<std::string> foam_getPointDataNames(std::string name);
// get mesh single point data named dataname
std::vector<double> foam_getPointData(std::string name,int label, std::string dataname);
// get mesh single face points
std::vector<int> foam_getFacePoints(std::string name, int label);
// get mesh faces points
std::vector<int> foam_getFacePoints(std::string name);

// get mesh single face data named dataname
std::vector<double> foam_getFaceData(std::string name,int label, std::string dataname);
// get mesh single cell points
std::vector<int> foam_getCellPoints(std::string name, int label);
// get mesh cells points
std::vector<int> foam_getCellPoints(std::string name);
// get data names associated to mesh cells
std::vector<std::string> foam_getCellDataNames(std::string name);
// get vector data names associated to mesh cells
std::vector<std::string> foam_getCellVectorDataNames(std::string name);
// get mesh single cell scalar data named dataname
double foam_getCellData(std::string name, int label,std::string dataname);
// get mesh cells scalar data named dataname
std::vector<double> foam_getCellData(std::string name, std::string dataname);
// get mesh single cell vector data named dataname
std::vector<double> foam_getCellVectorData(std::string name, int label,std::string dataname);
// get mesh cells vector data named dataname
std::vector<double> foam_getCellVectorData(std::string name, std::string dataname);
// get mesh boundary patch names
std::vector<std::string> foam_getBoundaryPatchNames(std::string name);
// get mesh boundary patches faces
std::vector<int> foam_getBoundaryPatchFaces(std::string name);
// update cell data (read from file to runTime object registry)
void foam_updateCellData(std::string name, std::string dataname);
// update cell vector data (read from file to runTime object registry)
void foam_updateCellVectorData(std::string name, std::string dataname);
// add mesh to wrapper registry
void foam_addMesh(std::string name,std::vector<double> points,  std::vector<int> cellpoints, std::vector<int> facepoints, std::vector<std::string> patchnames, std::vector<int> patchfaces, std::vector<std::string> patchtypes);
// write mesh
void foam_writeMesh(std::string name);
// deletes mesh from wrapper registry
void foam_deleteMesh(std::string name);
// set mesh single cell data (data must exists in registry)
void foam_setCellData(std::string name, int label,std::string dataname,double value);
// set mesh single cell vector data (data must exists in registry)
void foam_setCellVectorData(std::string name, int label,std::string dataname,std::vector<double> values);
// set mesh cell data
void foam_setCellData(std::string name, std::string dataname, std::vector<double> values);
// set mesh cell vector data
void foam_setCellVectorData(std::string name, std::string dataname, std::vector<double> values);
// write cell data 
void foam_writeCellData(std::string name,std::string dataname);
// write cell vector data
void foam_writeCellVectorData(std::string name,std::string dataname);
// get mesh faces count
int foam_getFaceCount(std::string name);
// get mesh points count
int foam_getPointCount(std::string name);
// get mesh cells count
int foam_getCellCount(std::string name);
// set active mesh
void foam_setMesh(std::string name);


// create default files
void foam_createDefaultFields(std::string name, std::string solver);
// modifying numerics fvSolution and fvSchemes
void foam_modifyNumerics(std::string name, std::string fvSch, std::string fvSol, std::string cD);
//modifying boundary conditions
void foam_setBC(std::string name, std::string fieldname, std::string dict);
//running the time steps
double foam_run(std::string name, int nproc, std::string solver);

template<class Type> 
GeometricField<Type, fvPatchField, volMesh> & find_Data(fvMesh& mesh,std::string& dataname)
{
  return const_cast< GeometricField < Type, fvPatchField, volMesh > &>(mesh.lookupObject< GeometricField <Type, fvPatchField, volMesh> >(word(dataname))); 
}


#endif

// ************************************************************************* //
