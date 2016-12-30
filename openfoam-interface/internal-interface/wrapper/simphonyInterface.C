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

#include "simphonyInterface.H"

#include "solution/solution.H"
#include "time/TimeMod.H"

#include "pointFields.H"
#include "volPointInterpolation.H"
#include "error.H"
#include "face.H"
#include "cell.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "preservePatchTypes.H"
#include "OFstream.H"

#include "RASModel.H"
#include "kEpsilon.H"
#include "laminar.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.H"
#include "relativeVelocityModel.H"
#include "dictionaryEntry.H"
#include "UList.H"
#include "mpi.h"
#include "CMULES.H"
#include "subCycle.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "hashedWordList.H"

#include <map>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::map<std::string,Foam::Time *> runTimes;


void foam_init(std::string caseName,std::string cD)
{

  
    dictionary controlDict_(dictionary::null,IStringStream
    (
        cD
    )());

    Foam::Time *runTime = new Foam::TimeMod(controlDict_);
    runTimes[caseName] = runTime;
}



void foam_init_IO(std::string caseName, std::string rootPath,std::string cD)
{
  FatalError.throwExceptions();
  dictionary controlDict_(dictionary::null,IStringStream
  			  (
  			   cD
  			   )());
  Foam::Time *runTime = new Foam::Time(controlDict_, fileName(rootPath), fileName(caseName));

  runTimes[caseName] = runTime;
}


void foam_readMesh(std::string name)
{  
    
  Foam::Time *runTime = runTimes[name];
  Foam::word regionName = Foam::fvMesh::defaultRegion;

  new fvMesh(        
	     Foam::IOobject
	     (
	      regionName,
	      runTime->constant(),
	      *runTime,
	      Foam::IOobject::MUST_READ,
	      Foam::IOobject::NO_WRITE,
	      true                 // this to register object
	      ));
  
}

const fvMesh& getMeshFromDb(std::string name)
{
  Foam::Time *runTime = runTimes[name];
  Foam::word regionName = Foam::fvMesh::defaultRegion;
  return runTime->db().parent().lookupObject<fvMesh>(regionName); 
}

std::vector<double> foam_getPointCoordinates(std::string name, int label)
  {

    const fvMesh & mesh = getMeshFromDb(name);

    const pointField & pp = mesh.points();

    std::vector<double> coordinates(3); 
    coordinates[0]=pp[label].x();
    coordinates[1]=pp[label].y();
    coordinates[2]=pp[label].z();

    return coordinates;

  }

std::vector<double> foam_getPointCoordinates(std::string name)
  {

    const fvMesh & mesh = getMeshFromDb(name);
    const pointField & pp = mesh.points();

    std::vector<double> retValue(3*pp.size()); 
    int k=0;
    for(int i=0;i<pp.size();i++) {
      retValue[k]=pp[i].x();
      retValue[k+1]=pp[i].y();
      retValue[k+2]=pp[i].z();
      k+=3;
    }

    return retValue;

  }


std::vector<int> foam_getFacePoints(std::string name, int label)
{
  const fvMesh & mesh = getMeshFromDb(name);
  const faceList & faces = mesh.faces();
  std::vector<int> points(faces[label].size());
  for (std::vector<int>::size_type i=0;i<points.size();i++) points[i] = faces[label][i];
  return points;
}


std::vector<int> get_cell_points_in_order(const fvMesh & mesh, int label)
{
  
  const cellShapeList& cellShapes = mesh.cellShapes();
  const labelList& cellPoints = cellShapes[label];

  std::vector<int> points(cellPoints.size()); 

  for (std::vector<int>::size_type i=0;i<points.size();i++) points[i] = cellPoints[i];

  return points;

 
}



std::vector<int> foam_getCellPoints(std::string name, int label)
{
  const fvMesh & mesh = getMeshFromDb(name);
  return get_cell_points_in_order(mesh, label);
}


std::vector<int> foam_getCellPoints(std::string name)
{
  
  const fvMesh & mesh = getMeshFromDb(name);
  const labelListList &cellPoints = mesh.cellPoints();

  int k=0;

  // find out vector size
  for(int i=0;i< cellPoints.size();i++)
    {
      k+=1+cellPoints[i].size();
    }
  

  std::vector<int> retValue(k);

  k=0;
  for(int i=0;i< cellPoints.size();i++)
    {
      int nnodes = cellPoints[i].size();
      retValue[k]=nnodes;
      k++;

      std::vector<int> points = get_cell_points_in_order(mesh, i);
      for (std::vector<int>::size_type ii=0;ii<points.size();ii++) {
	retValue[k] = points[ii];
	k++;

      }

    }
       
      return retValue;
      
}


std::vector<std::string> foam_getCellDataNames(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  wordList dataNames = mesh.names("volScalarField");
  
  std::vector<std::string> names;

  for (int i=0;i<dataNames.size();i++) {
    names.push_back(dataNames[i]);
  }

  return names;
}

std::vector<std::string> foam_getCellVectorDataNames(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  wordList vectorDataNames = mesh.names("volVectorField");
  
  std::vector<std::string> names;

  for (int i=0;i<vectorDataNames.size();i++) {
    names.push_back(vectorDataNames[i]);
  }

  return names;
}

std::vector<std::string> foam_getCellTensorDataNames(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  wordList tensorDataNames = mesh.names("volTensorField");
  
  std::vector<std::string> names;

  for (int i=0;i<tensorDataNames.size();i++) {
    names.push_back(tensorDataNames[i]);
  }

  return names;
}

volVectorField& find_vectorData(std::string name,std::string dataname)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    volVectorField& vF = find_Data<vector>(mesh,dataname);
    return const_cast<volVectorField&>(vF);  
}

volTensorField& find_tensorData(std::string name,std::string dataname)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    volTensorField& vF = find_Data<tensor>(mesh,dataname);
    return const_cast<volTensorField&>(vF);  
}
   
volScalarField& find_scalarData(std::string name,std::string dataname)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    volScalarField& vS = find_Data<scalar>(mesh,dataname);
    return const_cast<volScalarField&>(vS);  

}




double foam_getCellData(std::string name, int label,std::string dataname)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    hashedWordList names(mesh.names());

    if (names.contains(word(dataname)))
      {
	volScalarField& field = find_Data<scalar>(mesh,dataname);
	return field[label];
      }
    else 
      {
	new volScalarField
	  (
	   IOobject
	   (
	    word(dataname),
	    runTimes[name]->timeName(),
	    mesh,
	    IOobject::MUST_READ_IF_MODIFIED,
	    //        IOobject::MUST_READ,
	    IOobject::NO_WRITE
	    ),
	   mesh);
	volScalarField& field = find_Data<scalar>(mesh,dataname);
	return field[label];	
      }
  
}



std::vector<int> foam_getBoundaryCells(std::string name, std::string boundaryname)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));

    label patchID = mesh.boundaryMesh().findPatchID(word(boundaryname));
    const fvPatch& patch = mesh.boundary()[patchID];
    const labelUList& labs = patch.faceCells();
    std::vector<int> cellIDs(labs.size());
    int i=0;
    forAll(labs, lab)
      {
	cellIDs[i] = labs[lab];
	i += 1;
      }

    return cellIDs;

}


std::vector<double> foam_getCellVectorData(std::string name, int label,std::string dataname)
{
  fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
  hashedWordList names(mesh.names());

  if (names.contains(word(dataname)))
    {
      volVectorField& field = find_Data<vector>(mesh,dataname);
      std::vector<double> values(3);
      values[0] = field[label].x(); 
      values[1] = field[label].y(); 
      values[2] = field[label].z();
      return values;
    }else
    { 
      new volVectorField
	(       IOobject
		(
		 word(dataname),
		 runTimes[name]->timeName(),
		 mesh,
		 IOobject::MUST_READ_IF_MODIFIED,
		 //		 IOobject::MUST_READ,
		 IOobject::NO_WRITE
		 ),
		mesh);
      volVectorField& field = find_Data<vector>(mesh,dataname);

      std::vector<double> values(3);
      values[0] = field[label].x(); 
      values[1] = field[label].y(); 
      values[2] = field[label].z();
      return values;
    }
}




std::vector<double> foam_getCellTensorData(std::string name, int label,std::string dataname)
{
  fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
  hashedWordList names(mesh.names());
  if (names.contains(word(dataname)))
    {
      volTensorField& field = find_Data<tensor>(mesh,dataname);
      
      std::vector<double> values(9);
      values[0] = field[label].xx(); 
      values[1] = field[label].xy(); 
      values[2] = field[label].xz();
      values[3] = field[label].yx(); 
      values[4] = field[label].yy(); 
      values[5] = field[label].yz();
      values[6] = field[label].zx(); 
      values[7] = field[label].zy(); 
      values[8] = field[label].zz();
      return values;
    }else
    {
    
      new volTensorField
	(
	 IOobject
	 (
	  word(dataname),
	  runTimes[name]->timeName(),
	  mesh,
	  IOobject::MUST_READ_IF_MODIFIED,
	  //	  IOobject::MUST_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh);
      volTensorField& field = find_Data<tensor>(mesh,dataname);
      std::vector<double> values(9);
      values[0] = field[label].xx(); 
      values[1] = field[label].xy(); 
      values[2] = field[label].xz();
      values[3] = field[label].yx(); 
      values[4] = field[label].yy(); 
      values[5] = field[label].yz();
      values[6] = field[label].zx(); 
      values[7] = field[label].zy(); 
      values[8] = field[label].zz();
      return values;
    }
}



std::vector<std::string> foam_getBoundaryPatchNames(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
      
  wordList patchNames = mesh.boundaryMesh().names();
      
  std::vector<std::string> names(patchNames.size());
      
  int i=0;
  forAll(patchNames, patchI)
    {
      names[i] = patchNames[patchI];
      i+=1;
    }
    
  return names;

}

std::vector<int> foam_getBoundaryPatchFaces(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  int vecsize = 0;
  forAll(mesh.boundary(), patchI)
    {
      const faceUList& patchFaces =  mesh.boundaryMesh()[patchI];
      vecsize+=patchFaces.size()+1;
    }
  
  std::vector<int> patchfacelabels(vecsize);
  
  int ind = 0;
  forAll(mesh.boundaryMesh().names(), patchName) {
    
    polyPatch cPatch = mesh.boundaryMesh()[patchName];
    patchfacelabels[ind] = cPatch.size();
    forAll(cPatch, faceI) {
      label faceId = cPatch.start() + faceI;
      ind+=1;
      patchfacelabels[ind] =faceId;
      
    }
    ind+=1;
  }
  
  return patchfacelabels;
  
}




std::vector<int> foam_getFacePoints(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  const faceList &faces = mesh.faces();
  
  int k=0;
  
  for(int i=0;i< faces.size();i++)
    {
      k+=1+faces[i].size();
    }
  
  std::vector<int> retValue(k);

  k=0;
  for(int i=0;i< faces.size();i++)
    {
      retValue[k]=faces[i].size();
      k++;
      for (int j=0;j<faces[i].size();j++) {
	retValue[k]=faces[i][j];
	k++;
      }
    }  
  
  return retValue;
  
}
  

std::vector<double> foam_getCellTensorData(std::string name, std::string dataname)
{
  fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
  hashedWordList names(mesh.names());
  if (names.contains(word(dataname)))
    {
      volTensorField& field = find_Data<tensor>(mesh,dataname);
      std::vector<double> retValue(9*field.size()); 
      int k=0;
      for(int i=0;i<field.size();i++) 
	{
	  retValue[k]=field[i].xx();
	  retValue[k+1]=field[i].xy();
	  retValue[k+2]=field[i].xz();
	  retValue[k+3]=field[i].yx();
	  retValue[k+4]=field[i].yy();
	  retValue[k+5]=field[i].yz();
	  retValue[k+6]=field[i].zx();
	  retValue[k+7]=field[i].zy();
	  retValue[k+8]=field[i].zz();
	  k+=9;
	}
	    
      return retValue;      

    }else
    {
    
      new volTensorField
	(
	 IOobject
	 (
	  word(dataname),
	  runTimes[name]->timeName(),
	  mesh,
	  IOobject::MUST_READ_IF_MODIFIED,
	  //	  IOobject::MUST_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh);
      volTensorField& field = find_Data<tensor>(mesh,dataname);
      std::vector<double> retValue(9*field.size()); 
      int k=0;
      for(int i=0;i<field.size();i++) 
	{
	  retValue[k]=field[i].xx();
	  retValue[k+1]=field[i].xy();
	  retValue[k+2]=field[i].xz();
	  retValue[k+3]=field[i].yx();
	  retValue[k+4]=field[i].yy();
	  retValue[k+5]=field[i].yz();
	  retValue[k+6]=field[i].zx();
	  retValue[k+7]=field[i].zy();
	  retValue[k+8]=field[i].zz();
	  k+=9;
	}
	    
      return retValue;

    }


}

std::vector<double> foam_getCellVectorData(std::string name, std::string dataname)
{
  fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
  hashedWordList names(mesh.names());

  if (names.contains(word(dataname)))
    {
      volVectorField& field = find_Data<vector>(mesh,dataname);
      std::vector<double> retValue(3*field.size()); 
      int k=0;
      for(int i=0;i<field.size();i++) 
	{
	  retValue[k]=field[i].x();
	  retValue[k+1]=field[i].y();
	  retValue[k+2]=field[i].z();
	  k+=3;
	}
	    
      return retValue;
    }else
    { 
      new volVectorField
	(       IOobject
		(
		 word(dataname),
		 runTimes[name]->timeName(),
		 mesh,
		 IOobject::MUST_READ_IF_MODIFIED,
		 //		 IOobject::MUST_READ,
		 IOobject::NO_WRITE
		 ),
		mesh);
      volVectorField& field = find_Data<vector>(mesh,dataname);

      std::vector<double> retValue(3*field.size()); 
      int k=0;
      for(int i=0;i<field.size();i++) 
	{
	  retValue[k]=field[i].x();
	  retValue[k+1]=field[i].y();
	  retValue[k+2]=field[i].z();
	  k+=3;
	}
	    
      return retValue;

    }

}




std::vector<double> foam_getCellData(std::string name, std::string dataname)
{

   fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    hashedWordList names(mesh.names());

    if (names.contains(word(dataname)))
      {
	volScalarField& field = find_Data<scalar>(mesh,dataname);
	std::vector<double> retValue(field.size()); 
	for(int i=0;i<field.size();i++) 
	  {
	    retValue[i]=field[i];
	  }
	return retValue;	
      }
    else 
      {
	new volScalarField
	  (
	   IOobject
	   (
	    word(dataname),
	    runTimes[name]->timeName(),
	    mesh,
	    IOobject::MUST_READ_IF_MODIFIED,
	    //        IOobject::MUST_READ,
	    IOobject::NO_WRITE
	    ),
	   mesh);
	volScalarField& field = find_Data<scalar>(mesh,dataname);
	std::vector<double> retValue(field.size()); 
	for(int i=0;i<field.size();i++) 
	  {
	    retValue[i]=field[i];
	  }
	return retValue;	
      }
 
}


void foam_updateCellVectorData(std::string name, std::string dataname)
  {

    Foam::Time *runTime = runTimes[name];
    const fvMesh & mesh = getMeshFromDb(name);
    word vectorName(dataname);
    IOobject vectorHeader
      (
       vectorName,
       runTime->timeName(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE,
       true                 // this to register object
       );
    new volVectorField(vectorHeader, mesh);

  }


void foam_updateCellData(std::string name, std::string dataname)
  {

    Foam::Time *runTime = runTimes[name];
    const fvMesh & mesh = getMeshFromDb(name);
    word scalarName(dataname);
    
    IOobject scalarHeader
      (
       scalarName,
       runTime->timeName(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE,
       true                 // this to register object
       );
    
    new volScalarField(scalarHeader, mesh);


  }
  

bool same_points(const face &f1, const face &f2, const pointField & points)
{
  const pointField &points1 = f1.points(points);
  const pointField &points2 = f2.points(points);

  int i=0;
  forAll(points1, pi1)
    {
      forAll(points2,pi2)
	{
	  if (pi1 == pi2) {
	    i++;
	    break;
	  }
	}
    }
  
  return (i == points1.size());

}

void foam_addMesh(std::string name,std::vector<double> points,  std::vector<int> cellpoints, std::vector<int> facepoints, std::vector<std::string> patchnames, std::vector<int> patchfaces, std::vector<std::string> patchtypes)
{  
  
  // add first points
  pointField  foamPoints(points.size()/3);

  int pind = 0;
  for (std::vector<double>::size_type i=0;i<points.size();i+=3) {    
    foamPoints[pind].x() = points[i];
    foamPoints[pind].y() = points[i+1];
    foamPoints[pind].z() = points[i+2];
    pind+=1;
  }

 
  // find out number of cells
  
  int nOfCells = 0;
  std::vector<int>::size_type cind = 0;
  
  while (cind < cellpoints.size()-1) {
    nOfCells+= 1;
    cind+=cellpoints[cind]+1;
  }

 
  // create cellshape list
  cellShapeList cellShapes(nOfCells);
    
  const cellModel& hex = *(cellModeller::lookup("hex"));
  const cellModel& prism = *(cellModeller::lookup("prism"));
  const cellModel& pyr = *(cellModeller::lookup("pyr"));
  const cellModel& tet = *(cellModeller::lookup("tet"));
  
  cind = 0;
  int icell = 0;
  while (cind < cellpoints.size()) {
    labelList labels(cellpoints[cind]);
    int np = cellpoints[cind];
    for (int i=0;i<np;i++) 
      {
	cind+=1;
	labels[i] = cellpoints[cind];	
      }
    
    if (np == 4)
      cellShapes[icell] = cellShape(tet,labels);	
    else if (np == 5)
      cellShapes[icell] = cellShape(pyr,labels);
    else if (np == 6)
      cellShapes[icell] = cellShape(prism,labels);
    else if (np == 8) {
      cellShapes[icell] = cellShape(hex,labels);
    }
    else {
      FatalErrorIn("simphonyInterface:foam_addMesh")
	<< "Cell type with "<< np << " points not supported " << exit(FatalError);
      throw;
    }
    cind+=1;		 
    icell+=1;
    
  }
    

  Foam::Time *runTime = runTimes[name];
 
  Foam::word regionName = Foam::fvMesh::defaultRegion;
 
  const word defaultFacesName = "defaultFaces";
  word defaultFacesType = emptyPolyPatch::typeName;

  wordList patchNames(patchnames.size());
  wordList patchTypes(patchtypes.size());

  PtrList<dictionary> patchDicts;
  
  preservePatchTypes
    (
        *runTime,
        runTime->constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

  forAll(patchNames, patchI)
    {
      patchNames[patchI] = word(patchnames[patchI]);
    }
  forAll(patchTypes, patchI)
    {
        patchTypes[patchI] = word(patchtypes[patchI]);
    }

  // Add patch information to dictionary
  forAll(patchNames, patchI)
    {
      if (!patchDicts.set(patchI))
        {
	  patchDicts.set(patchI, new dictionary());
        }
      
      patchDicts[patchI].add("type", patchTypes[patchI], false);
    }
 
 // find out number of faces
    
  int nOfFaces = 0;
  std::vector<int>::size_type find = 0;
  
  while (find < facepoints.size()-1) {
    nOfFaces+= 1;
    find+=facepoints[find]+1;
  }

  // create faces list
  faceList faces(nOfFaces);
  find = 0;
  nOfFaces = 0;
  
  while (find < facepoints.size()-1) {
    int nOfPoints=facepoints[find];
    faces[nOfFaces].setSize(nOfPoints);
    for (int i=0;i<nOfPoints;i++)
      {
	find+=1;
	label pl(facepoints[find]);
	faces[nOfFaces][i] = pl;
      }
    nOfFaces+= 1;
    find+=1;
  }
 

  // make patch face lists

  faceListList patchFaces(patchnames.size());

  std::vector<int>::size_type pi = 0;
  
  std::vector<std::string>::size_type paind =0;
  
  while(paind < patchnames.size())
    {
      int nf = patchfaces[pi];
      
      faceList flist(nf);
      for (int i=0;i<nf;i++)
	{
	  pi+=1;
	  flist[i] = faces[patchfaces[pi]];
	}
      patchFaces[paind]=flist;
      pi+=1;
      paind+=1;
      
    }


 const Foam::IOobject  meshregIO(
			regionName,
			runTime->constant(),
			*runTime,
			Foam::IOobject::NO_READ,
			Foam::IOobject::NO_WRITE,
			true                // this to register object
			    );


  
 new fvMesh
     (
      meshregIO,
      xferMove(foamPoints),
      cellShapes,
      patchFaces,
      patchNames,
      patchDicts,
      defaultFacesName,
      defaultFacesType
      );
}


void foam_setCellData(std::string name, int label,std::string dataname,double value)
{    

  volScalarField& field = find_scalarData(name,dataname);

  field[label] = value;
}



void foam_setCellVectorData(std::string name, int label,std::string dataname,std::vector<double> values)
{
  volVectorField& field = find_vectorData(name, dataname);
  
  vector newValues(values[0],values[1],values[2]);

  field[label] = newValues;

}

void foam_setCellTensorData(std::string name, int label,std::string dataname,std::vector<double> values)
{
  volTensorField& field = find_tensorData(name, dataname);
  
  tensor newValues(values[0],values[1],values[2],
		   values[3],values[4],values[5],
		   values[6],values[7],values[8]);

  field[label] = newValues;

}

void foam_setCellData(std::string name, std::string dataname, std::vector<double> values) 
  {

    volScalarField& field = find_scalarData(name,dataname);
    field.internalField() = Field<scalar>(UList<scalar>(&(values[0]),values.size()));
    //    field.correctBoundaryConditions();

  }

void foam_setAndWriteCellData(std::string name,std::string dataname, std::vector<double> values,std::vector<int> dimension) 
  {

    const fvMesh & mesh = getMeshFromDb(name);
    hashedWordList names(mesh.names());
    if (names.contains(word(dataname)))
      {
	volScalarField& field = find_scalarData(name,dataname);
	field.internalField() = Field<scalar>(UList<scalar>(&(values[0]),values.size()));
	//	field.correctBoundaryConditions();
	field.write();
      }
    else
      {
	new volScalarField
	(
	 IOobject
	 (
	  word(dataname),
	  runTimes[name]->timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
	 dimensionedScalar(word(dataname), dimensionSet(dimension[0], dimension[1], dimension[2], dimension[3], dimension[4], dimension[5], dimension[6]), 0)
	 );
	volScalarField& field = find_scalarData(name,dataname);
	field.internalField() = Field<scalar>(UList<scalar>(&(values[0]),values.size()));
	//	field.correctBoundaryConditions();
	field.write();
 
      }
  }

void foam_setCellVectorData(std::string name, std::string dataname, std::vector<double> values) 
  {
    volVectorField& field = find_vectorData(name,dataname);
    field.internalField() = Field<vector>(UList<vector>((vector*)&(values[0]),values.size()/3));
    //    field.correctBoundaryConditions();

  }


void foam_setAndWriteCellVectorData(std::string name, std::string dataname, std::vector<double> values,std::vector<int> dimension) 
  {
  const fvMesh & mesh = getMeshFromDb(name);
  hashedWordList names(mesh.names());
  if (names.contains(word(dataname)))
    {
      volVectorField& field = find_vectorData(name,dataname);
      field.internalField() = Field<vector>(UList<vector>((vector*)&(values[0]),values.size()/3));
      //      field.correctBoundaryConditions();
      field.write();
    }
  else
    {
      new volVectorField
	(
	 IOobject
	 (
	  word(dataname),
	  runTimes[name]->timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
	 dimensionedVector(word(dataname), dimensionSet(dimension[0], dimension[1], dimension[2], dimension[3], dimension[4], dimension[5], dimension[6]), vector::zero)
);
      volVectorField& field = find_vectorData(name,dataname);

      field.internalField() = Field<vector>(UList<vector>((vector*)&(values[0]),values.size()/3));
      //      field.correctBoundaryConditions();
      field.write();
      
    }

  }

void foam_setCellTensorData(std::string name, std::string dataname, std::vector<double> values) 
  {
    volTensorField& field = find_tensorData(name,dataname);
    field.internalField() = Field<tensor>(UList<tensor>((tensor*)&(values[0]),values.size()/9));
    //    field.correctBoundaryConditions();
    foam_extendToBoundaries(name, dataname);
  }



void foam_setAndWriteCellTensorData(std::string name,  std::string dataname, std::vector<double> values,std::vector<int> dimension) 
  {
  const fvMesh & mesh = getMeshFromDb(name);
  hashedWordList names(mesh.names());

  if (names.contains(word(dataname)))
    {
      volTensorField& field = find_tensorData(name,dataname);
      field.internalField() = Field<tensor>(UList<tensor>((tensor*)&(values[0]),values.size()/9));
      //      field.correctBoundaryConditions();
      foam_extendToBoundaries(name, dataname);
      field.write();
    }
  else
    {
      new volTensorField
	(
	 IOobject
	 (
	  word(dataname),
	  runTimes[name]->timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
	 dimensionedTensor(word(dataname), dimensionSet(dimension[0], dimension[1], dimension[2], dimension[3], dimension[4], dimension[5], dimension[6]), tensor::zero)
);
      volTensorField& field = find_tensorData(name,dataname);
      field.internalField() = Field<tensor>(UList<tensor>((tensor*)&(values[0]),values.size()/9));
      //      field.correctBoundaryConditions();
      foam_extendToBoundaries(name, dataname);
      field.write();

    }
  }


void foam_writeMesh(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  mesh.write();
}

void foam_deleteMesh(std::string name)
{
  Foam::Time *runTime = runTimes[name];
  runTime->clear();
  std::map<std::string,Foam::Time *>::iterator iter = runTimes.find(name);
  if( iter != runTimes.end() )
    runTimes.erase( iter );
}

void foam_writeCellData(std::string name,std::string dataname)
{
  volScalarField& field = find_scalarData(name,dataname);
  field.write();
}

void foam_writeCellVectorData(std::string name,std::string dataname)
{
  volVectorField& field = find_vectorData(name,dataname);
  field.write();
}

int foam_getFaceCount(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  return mesh.faces().size();
}

int foam_getPointCount(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  return mesh.points().size();
}

int foam_getCellCount(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  return mesh.cells().size();
}

void foam_createDefaultFields(std::string name, std::string solver, bool io)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    if (io) {
      if(solver=="pimpleFoam" || solver=="simpleFoam"){
        #include "createFieldsPimpleFoamIO.H"
      }
      else if(solver=="interFoam"){
        #include "createFieldsInterFoamIO.H"
      }
      else if(solver=="driftFluxFoam"){
        #include "createFieldsDriftFluxSimphonyFoamIO.H"
      }
      else 
	FatalErrorIn("simphonyInterface:run") << "Solver "<<solver<<" not supported for IO wrapper" <<exit(FatalError);        

    }else {
      if(solver=="pimpleFoam" || solver=="simpleFoam"){
        #include "createFieldsPimpleFoam.H"
      }
      else if(solver=="driftFluxFoam"){
        #include "createFieldsDriftFluxSimphonyFoam.H"
      }
      else
	FatalErrorIn("simphonyInterface:run") << "Solver "<<solver<<" not supported for Internal wrapper" <<exit(FatalError);        

    }
}

void foam_writeFields(std::string name)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));

    std::vector<std::string> names = foam_getCellDataNames(name);

    for (std::vector<std::string>::size_type i=0;i<names.size();i++) 
      {
	volScalarField& vS = find_Data<scalar>(mesh,names[i]);
	vS.write();
      }
    
    std::vector<std::string> vnames = foam_getCellVectorDataNames(name);
    for (std::vector<std::string>::size_type i=0;i<vnames.size();i++) 
      {
	volVectorField& vV = find_Data<vector>(mesh,vnames[i]);
	vV.write();
      }

    
    std::vector<std::string> tnames = foam_getCellTensorDataNames(name);
    for (std::vector<std::string>::size_type i=0;i<tnames.size();i++) 
      {
	volTensorField& vT = find_Data<tensor>(mesh,tnames[i]);
	vT.write();
      }

}

void foam_readFields(std::string name)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    std::vector<std::string> names = foam_getCellDataNames(name);

    
    for (std::vector<std::string>::size_type i=0;i<names.size();i++) 
      {

	volScalarField& vS = find_Data<scalar>(mesh,names[i]);
	volScalarField* vSF = new volScalarField
	  (
	   IOobject
	   (
	    word(names[i]),
	    runTimes[name]->timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE,
	    true
	    ),
	   mesh
	   );
	vS.internalField() = (*vSF).internalField();
	vS.boundaryField() = (*vSF).boundaryField();

      }
    
    std::vector<std::string> vnames = foam_getCellVectorDataNames(name);
    for (std::vector<std::string>::size_type i=0;i<vnames.size();i++) 
      {
	volVectorField& vV = find_Data<vector>(mesh,vnames[i]);
	volVectorField* vVF = new volVectorField
	  (
	   IOobject
	   (
	    word(vnames[i]),
	    runTimes[name]->timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE,
	    true
	    ),
	   mesh
	   );
	vV.internalField() = (*vVF).internalField();
	vV.boundaryField() = (*vVF).boundaryField();
      }

    
    std::vector<std::string> tnames = foam_getCellTensorDataNames(name);
    for (std::vector<std::string>::size_type i=0;i<tnames.size();i++) 
      {
	volTensorField& vT = find_Data<tensor>(mesh,tnames[i]);
	volTensorField* vVT = new volTensorField
	  (
	   IOobject
	   (
	    word(tnames[i]),
	    runTimes[name]->timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE,
	    true
	    ),
	   mesh
	   );
	vT.internalField() = (*vVT).internalField();
	vT.boundaryField() = (*vVT).boundaryField();

      }


}

void foam_writeNumerics(std::string name)
{
    const fvMesh & mesh = getMeshFromDb(name);

    dictionary& fvSchemesDict = const_cast<dictionary&>(mesh.schemesDict());
    dictionary& fvSolutionDict = const_cast<dictionary&>(mesh.solutionDict());
    dictionary& controlDict = const_cast<dictionary&>(runTimes[name]->controlDict());

    IOdictionary fvSc(
       IOobject
       (
        "fvSchemes",
        runTimes[name]->system(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
	),
       fvSchemesDict
       );
    fvSc.regIOobject::write();
 
    IOdictionary fvSo(
       IOobject
       (
        "fvSolution",
        runTimes[name]->system(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
	),
       fvSolutionDict
       );
    fvSo.regIOobject::write(); 
     
    IOdictionary cD(
       IOobject
       (
        "controlDict",
        runTimes[name]->system(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
	),
       controlDict
       );
    cD.regIOobject::write();     
     
    dictionary& TPDict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word("transportProperties"))); 

    IOdictionary tP(
       IOobject
       (
        "transportProperties",
        runTimes[name]->constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
	),
       TPDict
       );
    tP.regIOobject::write();     

}

void foam_writeDictionary(std::string name, std::string dictionaryName, bool constant)
{
  const fvMesh & mesh = getMeshFromDb(name);
     
  dictionary& dict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word(dictionaryName))); 

    if (constant) {
      IOdictionary d(
		     IOobject
		     (
		      dictionaryName,
		      runTimes[name]->constant(),
		      mesh,
		      IOobject::NO_READ,
		      IOobject::NO_WRITE
		      ),
		     dict
		     );
      d.regIOobject::write();     
    }
    else {
      IOdictionary d(
		     IOobject
		     (
		      dictionaryName,
		      runTimes[name]->system(),
		      mesh,
		      IOobject::NO_READ,
		      IOobject::NO_WRITE
		      ),
		     dict
		     );
    d.regIOobject::write();     
    }
}

void foam_writePathDictionary(std::string path,std::string dictionaryName,std::string head, std::string dictionaryContent)
{
    dictionary dict_(dictionary::null,IStringStream
    (
        dictionaryContent.c_str()
    )());

    fileName outputFile(fileName(path)/fileName(dictionaryName));
    OFstream os(outputFile);
    
    os <<head.c_str()<<endl;
    os <<dict_<<endl;

}



void foam_modifyNumerics(std::string name, std::string fvSch, std::string fvSol, std::string cD, std::string TP, bool io)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    
    IStringStream fvSchIS(fvSch.c_str());
    IStringStream fvSolIS(fvSol.c_str());
    IStringStream cDIS(cD.c_str());
    IStringStream TPIS(TP.c_str());

    dictionary& fvSchemesDict = const_cast<dictionary&>(mesh.schemesDict());
    dictionary& fvSolutionDict = const_cast<dictionary&>(mesh.solutionDict());
    dictionary& controlDict = const_cast<dictionary&>(runTimes[name]->controlDict());

    fvSchemesDict.read
    (
        fvSchIS()
    );
    fvSchemes* fvSc = &mesh;
    fvSc->read(true);

    fvSolutionDict.read
    (
        fvSolIS()
    );
    fvSolution* fvSo = &mesh;
    fvSo->read(true);

    controlDict.read
    (
        cDIS()
    );

    if (!io)
      runTimes[name]->read();

    dictionary& TPDict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word("transportProperties"))); 
    TPDict.read
    (
        TPIS()        
    );  

    dictionary& RPDict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word("RASProperties"))); 
    RPDict.read
    (
     IStringStream
     (
      "RASModel            laminar;"
      "turbulence            off;"      
      )()
    );    
}

void foam_modifyDictionary(std::string name, std::string dictionaryName, std::string dictionaryContent)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    
    IStringStream dictIS(dictionaryContent.c_str());
 
    dictionary& dict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word(dictionaryName))); 
    dict.read
    (
        dictIS()
    );    
}

void foam_modifyUniformDimensionedVectorField(std::string name, std::string fieldName, std::vector<double> value)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    
    uniformDimensionedVectorField& field = const_cast<uniformDimensionedVectorField&>(mesh.lookupObject<uniformDimensionedVectorField>(word(fieldName))); 
 
    vector newValue(value[0],value[1],value[2]);
    field.value() = newValue;
 
}


void foam_setBC(std::string name, std::string fieldname, std::string dict)
{

      fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));


	IStringStream dictIS(dict.c_str());

	Time& runTime = *(runTimes[name]);

	IOdictionary IOdict(runTime, dictIS);

	if(fieldname=="U"){
		volVectorField& vF = find_Data<vector>(mesh,fieldname);
		vF.boundaryField().readField(vF,IOdict);
	}
	
	if(fieldname=="p" || fieldname=="p_rgh" || fieldname=="alpha.phase1"){
		volScalarField& vS = find_Data<scalar>(mesh,fieldname);
		vS.boundaryField().readField(vS,IOdict);
	}

}


void foam_extendToBoundaries(std::string name, std::string fieldname) {

  fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
  volTensorField& vT = find_Data<tensor>(mesh,fieldname);
  forAll(vT.boundaryField(), patchi) {
    forAll(vT.boundaryField()[patchi], facei) {
      vT.boundaryField()[patchi][facei] = vT.internalField()[mesh.boundaryMesh()[patchi].faceCells()[facei]];
    }
  }
}

double foam_run(std::string name, int nproc, std::string solver){
    if(nproc<2){
        Time& runTime = *(runTimes[name]);
        fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
        if(solver=="pimpleFoam"){
            #include "pimpleFoam.H"
        }
	if(solver=="simpleFoam") {
            #include "simpleFoam.H"
	}
        if(solver=="driftFluxFoam"){
            #include "driftFluxSimphonyFoam.H"
        }
        return runTime.timeOutputValue();
    }else{
        FatalErrorIn("simphonyInterface:run") << "Parallel OpenFOAM Internal Interface not implemented yet" <<exit(FatalError);
	return 0.0;
    }

}

void foam_updateData(std::string name, double time){
  runTimes[name]->setTime(time,0);
  foam_readFields(name); 
}

void foam_updateTime(std::string name, double time){
  runTimes[name]->setTime(time,0);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

