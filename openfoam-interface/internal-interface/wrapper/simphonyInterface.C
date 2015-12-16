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
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "error.H"
#include "face.H"
#include "cell.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "preservePatchTypes.H"

//CIMEC adds 
#include "RASModel.H"
#include "kEpsilon.H"
#include "laminar.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.H"
#include "relativeVelocityModel.H"
#include "dictionaryEntry.H"
#include "UList.H"
#include "mpi.h"
#include "fvSchemes/fvSchemes.H"
#include "CMULES.H"
#include "subCycle.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::map<std::string,Foam::Time *> runTimes;
//std::map<std::string,Foam::fvMesh *> meshes;


void foam_init(std::string caseName,std::string cD)
{


    dictionary controlDict_(dictionary::null,IStringStream
    (
        cD
    )());

    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time *runTime = new Foam::Time(controlDict_, "..", ".");
    runTimes[caseName] = runTime;

}

void foam_init_IO(std::string caseName, std::string rootPath)
{
  FatalError.throwExceptions();
  word name = Foam::Time::controlDictName;
  Foam::Time *runTime = new Foam::Time(name, fileName(rootPath), fileName(caseName));
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

    //    const fvMesh & mesh = getMeshFromDb(name);
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



std::vector<int> order_cellPoints(int i, const fvMesh & mesh, const labelListList &cellPoints, const cellList &cellFaces, const faceList &facePoints, const labelListList &pointFaces,const labelListList &pointPoints)
{
  int nfaces = cellFaces[i].size();
  int nnodes = cellPoints[i].size();

  int k=0;
  std::vector<int> retValue(nnodes);
  std::string type="";
  if (nfaces==6)
    {
      if (nnodes==8)
	type="hex";
      else
	type="wedge";
    }
  else if (nfaces==5)
    {
      if (nnodes==5)
	type="pyr";
      else
	type="prism";
    }
  else if (nfaces==4)
    {
      if (nnodes==4)
	type="tet";
      else
	type="tetWedge";
    }
  
  
  if (type == "tet") {
    
    for (int j=0;j<nnodes;j++) {
      
      retValue[k]=cellPoints[i][j];
    }
    
  }
  else if (type == "hex" ) {
    
    // add first face points
    label fl = cellFaces[i][0];
    labelList fp = facePoints[fl];
    
    // face point ordering is clockwise, must be reversed    
    for (int j=fp.size()-1;j>-1;j--) {
      retValue[k]=fp[j];
      k++;
    }
    
    // find opposite points
    
    // opposite side has same number of points
    
    for (int l=fp.size()-1;l>-1;l--) {
      
      labelList conPoints = pointPoints[fp[l]];
      
      for (int j=0;j<conPoints.size();j++) {
	// point must be in current cell
	if (findIndex(cellPoints[i],conPoints[j]) > -1) {
	  labelList ff = pointFaces[conPoints[j]];
	  
	  // if first face not on list the point is on opposite side
	  // 
	  label ffl = findIndex(ff,fl);
	  
	  if (ffl == -1) {
	    retValue[k]=conPoints[j];
	    k++;
	    break;
	  }
	}
      }
    }
  }else if (type == "prism") {
    
    // find face which has three nodes
    label fl=-1;
    for (int j=0;j<nfaces;j++) {
      if (facePoints[cellFaces[i][j]].size() == 3) {
	fl = cellFaces[i][j];
	break;
      }
    }
    
    labelList fp = facePoints[fl];
    
    // add first face points (reversed order)
    for (int j=fp.size()-1;j>-1;j--) {
      retValue[k]=fp[j];
      k++;
    }
    
    // find opposite points
    
    // opposite side has same number of points
    
    for (int l=fp.size()-1;l>-1;l--) {
      
      labelList conPoints = pointPoints[fp[l]];
      
      for (int j=0;j<conPoints.size();j++) {
	// point must be in current cell
	if (findIndex(cellPoints[i],conPoints[j]) > -1) {
	  labelList ff = pointFaces[conPoints[j]];
	  
	  // if first face not on list the point is on opposite side
	  // 
	  label ffl = findIndex(ff,fl);
	  
	  if (ffl == -1) {
	    retValue[k]=conPoints[j];
	    k++;
	    break;
	  }
	}
      }
    }
    
  }else {
    
    std::string except = "Cell type "+type+" not supported yet";
    FatalErrorIn("simphonyInterface:foam_getCellPoints()") << except <<exit(FatalError);
    throw;
  }
  return retValue;
  
}




std::vector<int> get_ordered_points(const fvMesh & mesh, const labelListList &cellPoints, const cellList &cellFaces, int label)
{
  int nfaces = cellFaces[label].size();
  int nnodes = cellPoints[label].size();

  std::string type="";

  if (nfaces==6)
    {
      if (nnodes==8)
	type="hex";
      else
	type="wedge";
    }
  else if (nfaces==5)
    {
      if (nnodes==5)
	type="pyr";
      else
	type="prism";
    }
  else if (nfaces==4)
    {
      if (nnodes==4)
	type="tet";
      else
	type="tetWedge";
    }
  

  const cellModel& cellModel = *(cellModeller::lookup(type));

  cellShape cShape = cellShape(cellModel,cellPoints[label]);

  std::vector<int> points(cShape.nPoints()); 

  for (std::vector<int>::size_type i=0;i<points.size();i++) points[i] = cShape[i];

  return points;

 
}


std::vector<int> foam_getCellPoints(std::string name, int label)
{
  const fvMesh & mesh = getMeshFromDb(name);
  const labelListList &cellPoints = mesh.cellPoints();
  const cellList &cellFaces = mesh.cells();
  // gives vertices connected to vertex
  const labelListList &pointPoints = mesh.pointPoints();
  const labelListList &pointFaces = mesh.pointFaces();
  const faceList &facePoints = mesh.faces();

  return order_cellPoints(label, mesh, cellPoints, cellFaces, facePoints, pointFaces, pointPoints);

}


std::vector<int> foam_getCellPoints(std::string name)
{
  
  const fvMesh & mesh = getMeshFromDb(name);
  const labelListList &cellPoints = mesh.cellPoints();
  const cellList &cellFaces = mesh.cells();
  // gives vertices connected to vertex
  const labelListList &pointPoints = mesh.pointPoints();
  const labelListList &pointFaces = mesh.pointFaces();
  const faceList &facePoints = mesh.faces();

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

      std::vector<int> points =  order_cellPoints(i, mesh, cellPoints, cellFaces, facePoints, pointFaces, pointPoints);
  //      std::vector<int> points = get_ordered_points(mesh, cellPoints, cellFaces, i);

      for (std::vector<int>::size_type i=0;i<points.size();i++) {

	retValue[k] = points[i];
	k++;

      }

    }
       
      return retValue;
      
}


std::vector<std::string> foam_getCellDataNames(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  wordList dataNames = mesh.names("volScalarField");
  
  std::vector<std::string> names(dataNames.size());

  for (int i=0;i<dataNames.size();i++) {
    names.push_back(dataNames[i]);
  }

  return names;
}

std::vector<std::string> foam_getCellVectorDataNames(std::string name)
{
  const fvMesh & mesh = getMeshFromDb(name);
  wordList vectorDataNames = mesh.names("volVectorField");
  
  std::vector<std::string> names(vectorDataNames.size());

  for (int i=0;i<vectorDataNames.size();i++) {
    names.push_back(vectorDataNames[i]);
  }

  return names;
}

volVectorField& find_vectorData(std::string name,std::string dataname)
{
  const fvMesh & mesh = getMeshFromDb(name);
  return const_cast<volVectorField&>(mesh.lookupObject<volVectorField>(word(dataname))); 
}
   
volScalarField& find_scalarData(std::string name,std::string dataname)
{
  const fvMesh & mesh = getMeshFromDb(name);
  return const_cast<volScalarField&>(mesh.lookupObject<volScalarField>(word(dataname)));  
}



double foam_getCellData(std::string name, int label,std::string dataname)
{
  volScalarField& field = find_scalarData(name,dataname);
  return field[label];
}

std::vector<double> foam_getCellVectorData(std::string name, int label,std::string dataname)
{
  volVectorField& field = find_vectorData(name,dataname);

  std::vector<double> values(3);
  values[0] = field[label].x(); 
  values[1] = field[label].y(); 
  values[2] = field[label].z();
  return values;

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
  

std::vector<double> foam_getCellVectorData(std::string name, std::string dataname)
{

  volVectorField vect = find_vectorData(name,dataname);
  std::vector<double> retValue(3*vect.size()); 
  int k=0;
  for(int i=0;i<vect.size();i++) 
    {
      retValue[k]=vect[i].x();
      retValue[k+1]=vect[i].y();
      retValue[k+2]=vect[i].z();
      k+=3;
    }
	    
  return retValue;

}

std::vector<double> foam_getCellData(std::string name, std::string dataname)
{
  volScalarField scal = find_scalarData(name,dataname);
  std::vector<double> retValue(scal.size()); 
  for(int i=0;i<scal.size();i++) 
    {
      retValue[i]=scal[i];
    }
  return retValue;
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
       IOobject::MUST_READ,
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
       IOobject::MUST_READ,
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

void foam_setCellData(std::string name, std::string dataname, std::vector<double> values) 
  {

    volScalarField& field = find_scalarData(name,dataname);
    field.internalField() = Field<scalar>(UList<scalar>(&(values[0]),values.size()));

  }




void foam_setCellVectorData(std::string name, std::string dataname, std::vector<double> values) 
  {

    volVectorField& field = find_vectorData(name,dataname);

    field.internalField() = Field<vector>(UList<vector>((vector*)&(values[0]),values.size()/3));

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

void foam_createDefaultFields(std::string name, std::string solver)
{
    const fvMesh & mesh = getMeshFromDb(name);
    if(solver=="pimpleFoam"){
        #include "createFieldsPimpleFoam.H"
    }
    if(solver=="driftFluxSimphonyFoam"){
        #include "createFieldsDriftFluxSimphonyFoam.H"
    }
    return;
}

void foam_modifyNumerics(std::string name, std::string fvSch, std::string fvSol, std::string cD)
{
    fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
    runTimes[name]->setTime(0.0,0); // restarting times
    
    IStringStream fvSchIS(fvSch.c_str());
    IStringStream fvSolIS(fvSol.c_str());
    IStringStream cDIS(cD.c_str());

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

    TimeMod* t = static_cast<TimeMod*>(runTimes[name]);
    t->update();

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
	
	if(fieldname=="p" || fieldname=="p_rgh" || fieldname=="alpha1"){
		volScalarField& vS = find_Data<scalar>(mesh,fieldname);
		vS.boundaryField().readField(vS,IOdict);
	}

}


double foam_run(std::string name, int nproc, std::string solver){
    if(nproc<2){
        Time& runTime = *(runTimes[name]);
        fvMesh & mesh = const_cast<fvMesh&>(getMeshFromDb(name));
        if(solver=="pimpleFoam"){
            #include "pimpleFoam.H"
        }
        if(solver=="driftFluxSimphonyFoam"){
            #include "driftFluxSimphonyFoam.H"
        }
        return runTime.timeOutputValue();
    }else{
    
        FatalErrorIn("simphonyInterface:run") << "Parallel OpenFOAM Internal Interface not implemented yet" <<exit(FatalError);        

        return 0.0;         
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

