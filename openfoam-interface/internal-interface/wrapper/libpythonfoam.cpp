#define WM_DP 1
#define NoRepository 1
#include <sstream>
#include "Python.h"
#include "simphonyInterface.H"
//#include "mpi.h"


extern "C" {


  static PyObject* init(PyObject *self, PyObject *args)
  {
     char *name;
     char *cD;

     if (!PyArg_ParseTuple(args,"ss",&name,&cD)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_init(std::string(name), std::string(cD));
      return Py_BuildValue("");
    }
    
    
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
    
  } 


  static PyObject* init_IO(PyObject *self, PyObject *args)
  {   
     char *name;
     char *path;
     char *cD;

     //     if (!PyArg_ParseTuple(args,"ss",&name,&path)) {
     if (!PyArg_ParseTuple(args,"sss",&name,&path,&cD)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      //      foam_init_IO(std::string(name),std::string(path));
      foam_init_IO(std::string(name),std::string(path),std::string(cD));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 



  static PyObject* readMesh(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_readMesh(std::string(name));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getPointCoordinates(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    if (!PyArg_ParseTuple(args,"si",&name,&label)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getPointCoordinates(std::string(name), label);
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	 
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getAllPointCoordinates(PyObject *self, PyObject *args)
  {
    char *name;

    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getPointCoordinates(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	 
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }



  static PyObject* getFacePoints(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    
    if (!PyArg_ParseTuple(args,"si",&name,&label)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<int> values = foam_getFacePoints(std::string(name), label);
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<int>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("i",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getAllFacePoints(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<int> values = foam_getFacePoints(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<int>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("i",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }



  static PyObject* getCellPoints(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    
    if (!PyArg_ParseTuple(args,"si",&name,&label)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<int> values = foam_getCellPoints(std::string(name), label);
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<int>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("i",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getAllCellPoints(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<int> values = foam_getCellPoints(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<int>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("i",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }



  static PyObject* getCellDataNames(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<std::string> values = foam_getCellDataNames(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<std::string>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("s",values[i].c_str());
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getCellVectorDataNames(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<std::string> values = foam_getCellVectorDataNames(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<std::string>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("s",values[i].c_str());
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getCellTensorDataNames(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<std::string> values = foam_getCellTensorDataNames(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<std::string>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("s",values[i].c_str());
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }


  static PyObject* getCellData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    if (!PyArg_ParseTuple(args,"sis",&name,&label,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      double value = foam_getCellData(std::string(name), label,std::string(dataname));
      return Py_BuildValue("d",value);
      
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }  



  static PyObject* getBoundaryCells(PyObject *self, PyObject *args)
  {
    char *name;
    char *boundaryname;
    if (!PyArg_ParseTuple(args,"ss",&name,&boundaryname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<int> values = foam_getBoundaryCells(std::string(name),std::string(boundaryname));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<int>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("i",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }  


    
  static PyObject* getAllCellData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;

    if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getCellData(std::string(name), std::string(dataname));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }



  static PyObject* getCellVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;

    if (!PyArg_ParseTuple(args,"sis",&name,&label,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getCellVectorData(std::string(name), label, std::string(dataname));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
      
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getCellTensorData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;

    if (!PyArg_ParseTuple(args,"sis",&name,&label,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getCellTensorData(std::string(name), label, std::string(dataname));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
      
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }


  static PyObject* getAllCellVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;

    if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getCellVectorData(std::string(name), std::string(dataname));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }



  static PyObject* getAllCellTensorData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;

    if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getCellTensorData(std::string(name), std::string(dataname));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<double>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("d",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }




  static PyObject* getBoundaryPatchNames(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<std::string> values = foam_getBoundaryPatchNames(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<std::string>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("s",values[i].c_str());
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getBoundaryPatchFaces(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<int> values = foam_getBoundaryPatchFaces(std::string(name));
      PyObject *pylist = PyList_New(values.size());
      if (pylist != NULL) {
	for (std::vector<int>::size_type i=0;i<values.size();i++) {
	  PyObject *item = Py_BuildValue("i",values[i]);
	  PyList_SetItem(pylist, i, item);
	}
	return pylist;
      }else
	return NULL;	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* updateCellData(PyObject *self, PyObject *args)
  {
     char *name;
     char *dataname;
     if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_updateCellData(std::string(name),std::string(dataname));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 

  static PyObject* updateCellVectorData(PyObject *self, PyObject *args)
  {
     char *name;
     char *dataname;
     if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_updateCellVectorData(std::string(name),std::string(dataname));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 

  static PyObject* addMesh(PyObject *self, PyObject *args)
  {
    char *name;
    PyObject *cells;
    PyObject *faces;
    PyObject *points;
    PyObject *patchnames;
    PyObject *patchfaces;
    PyObject *patchtypes;
 
    if (!PyArg_ParseTuple(args,"sO!O!O!O!O!O!",&name,&PyList_Type,&points,&PyList_Type,&cells,&PyList_Type,&faces,&PyList_Type,&patchnames,&PyList_Type,&patchfaces,&PyList_Type,&patchtypes)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }    

   try{
     int lenpoints = PyList_Size(points);
     int lencells = PyList_Size(cells);
     int lenfaces = PyList_Size(faces);
     int lenpatchnames = PyList_Size(patchnames);
     int lenpatchfaces = PyList_Size(patchfaces);
     int lenpatchtypes = PyList_Size(patchtypes);
     
     std::vector<double> foampoints(lenpoints);
     std::vector<int> foamcells(lencells);
     std::vector<int> foamfaces(lenfaces);
     std::vector<std::string> foampatchnames(lenpatchnames);
     std::vector<std::string> foampatchtypes(lenpatchtypes);
     std::vector<int> foampatchfaces(lenpatchfaces);


     PyObject *strObj;
     for (int j=0;j<lenpoints;j++)
       {
	 strObj = PyList_GetItem(points,j);
	 foampoints[j] = PyFloat_AsDouble(strObj);
       }
     
     for (int j=0;j<lencells;j++)
       {
	 strObj = PyList_GetItem(cells,j);
	 foamcells[j] = PyInt_AsLong(strObj);
       }

     for (int j=0;j<lenfaces;j++)
       {
	 strObj = PyList_GetItem(faces,j);
	 foamfaces[j] = PyInt_AsLong(strObj);
       }

     for (int j=0;j<lenpatchfaces;j++)
       {
	 strObj = PyList_GetItem(patchfaces,j);
	 foampatchfaces[j] = PyInt_AsLong(strObj);
       }
  

     char *pname;
     for (int j=0;j<lenpatchnames;j++)
       {
	 strObj = PyList_GetItem(patchnames,j);
	 pname = PyString_AsString(strObj);
	 foampatchnames[j] = std::string(pname);
       }

     char *ptype;
     for (int j=0;j<lenpatchtypes;j++)
       {
	 strObj = PyList_GetItem(patchtypes,j);
	 ptype = PyString_AsString(strObj);
	 foampatchtypes[j] = std::string(ptype);
       }
     
         
     foam_addMesh(std::string(name),foampoints, foamcells, foamfaces, foampatchnames, foampatchfaces, foampatchtypes);
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* writeMesh(PyObject *self, PyObject *args)
  {
     char *name;
     if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_writeMesh(std::string(name));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 

  static PyObject* deleteMesh(PyObject *self, PyObject *args)
  {
     char *name;
     if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_deleteMesh(std::string(name));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 



 static PyObject* setCellData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    double value;
    if (!PyArg_ParseTuple(args,"sisd",&name,&label,&dataname,&value)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {

      foam_setCellData(std::string(name), label, std::string(dataname), value);
      return Py_BuildValue("");
	    	 
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

 static PyObject* setAllCellData(PyObject *self, PyObject *args)
  {
    char *name;
    int write;
    char *dataname;

    PyObject *values;
    PyObject *dime;
    if (!PyArg_ParseTuple(args,"ssiO!O!",&name,&dataname,&write,&PyList_Type,&values,&PyList_Type,&dime)) {
	PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
	return NULL;
    }
    if (!write)
      {

	try {
	  PyObject * strObj;
	  int valuessize = PyList_Size(values);
	  
	  std::vector<double> vals(valuessize);
	  for (int i=0;i<valuessize;i++) {
	    strObj = PyList_GetItem(values, i);
	    vals[i] = PyFloat_AsDouble(strObj);
	  }
	    
	  foam_setCellData(std::string(name),  std::string(dataname), vals);
	  return Py_BuildValue("");
	    
	}
	catch (Foam::error& fErr)
	  {
	    PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	    return NULL;
	  }
	catch (std::exception& e) {
	  PyErr_SetString(PyExc_RuntimeError,e.what());
	  return NULL;
	}
	catch (...) {
	  PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
	  return NULL;
	}
      }else
      {
    
	try {
	  PyObject * strObj;
	  int valuessize = PyList_Size(values);
	  
	  std::vector<double> vals(valuessize);
	  for (int i=0;i<valuessize;i++) {
	    strObj = PyList_GetItem(values, i);
	    vals[i] = PyFloat_AsDouble(strObj);
	  }


	  int dimensionsize = PyList_Size(dime);
	  std::vector<int> dimension(dimensionsize);
	  for (int i=0;i<dimensionsize;i++) {
	    strObj = PyList_GetItem(dime, i);
	    dimension[i] = PyInt_AsLong(strObj);
	  }
	  
	  foam_setAndWriteCellData(std::string(name), std::string(dataname), vals,dimension);
	  return Py_BuildValue("");
	  
	}
	catch (Foam::error& fErr)
	  {
	    PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	    return NULL;
	  }
	catch (std::exception& e) {
	  PyErr_SetString(PyExc_RuntimeError,e.what());
	  return NULL;
	}
	catch (...) {
	  PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
	  return NULL;
	}
      }
  }


static PyObject* setAllCellVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    int write;
    char *dataname;
    PyObject *values;
    PyObject *dime;

    if (!PyArg_ParseTuple(args,"ssiO!O!",&name,&dataname,&write,&PyList_Type,&values,&PyList_Type,&dime)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    if (!write) {
      try
	{
	  PyObject * strObj;
	  int valuessize = PyList_Size(values);
	  
	  std::vector<double> vals(valuessize);
	  for (int i=0;i<valuessize;i++) {
	    strObj = PyList_GetItem(values, i);
	    vals[i] = PyFloat_AsDouble(strObj);
	  }

	  foam_setCellVectorData(std::string(name), std::string(dataname),vals);
	  return Py_BuildValue("");
	  
	}
      catch (Foam::error& fErr)
	{
	  PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	  return NULL;
	}
      catch (std::exception& e) {
	PyErr_SetString(PyExc_RuntimeError,e.what());
	return NULL;
      }
      catch (...) {
	PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
	return NULL;
      }
    }
    else
      {
	try {
	  PyObject * strObj;
	  int valuessize = PyList_Size(values);
	  
	  std::vector<double> vals(valuessize);
	  for (int i=0;i<valuessize;i++) {
	    strObj = PyList_GetItem(values, i);
	    vals[i] = PyFloat_AsDouble(strObj);
	  }
	  int dimensionsize = PyList_Size(dime);
	  std::vector<int> dimension(dimensionsize);
	  for (int i=0;i<dimensionsize;i++) {
	    strObj = PyList_GetItem(dime, i);
	    dimension[i] = PyInt_AsLong(strObj);
	  }

	  foam_setAndWriteCellVectorData(std::string(name), std::string(dataname),vals,dimension);
	  return Py_BuildValue("");
	  
	}
	catch (Foam::error& fErr)
	  {
	    PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	    return NULL;
	  }
	catch (std::exception& e) {
	  PyErr_SetString(PyExc_RuntimeError,e.what());
	  return NULL;
	}
	catch (...) {
	  PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
	  return NULL;
	}
      }
  }

  static PyObject* setAllCellTensorData(PyObject *self, PyObject *args)
  {
    char *name;
    int write;
    char *dataname;
    PyObject *values;
    PyObject *dime;

    if (!PyArg_ParseTuple(args,"ssiO!O!",&name,&dataname,&write,&PyList_Type,&values,&PyList_Type,&dime)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    if (!write)
      {

	try
	  {
	  PyObject * strObj;
	  int valuessize = PyList_Size(values);
	  
	  std::vector<double> vals(valuessize);
	  for (int i=0;i<valuessize;i++) {
	    strObj = PyList_GetItem(values, i);
	    vals[i] = PyFloat_AsDouble(strObj);
	  }
	  foam_setCellTensorData(std::string(name), std::string(dataname),vals);
	    return Py_BuildValue("");
	    
	  }
	catch (Foam::error& fErr)
	  {
	    PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	    return NULL;
	  }
	catch (std::exception& e) {
	  PyErr_SetString(PyExc_RuntimeError,e.what());
	  return NULL;
	}
	catch (...) {
	  PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
	  return NULL;
	}
      }
    else
      {
	try {
	  PyObject * strObj;
	  int valuessize = PyList_Size(values);
	  
	  std::vector<double> vals(valuessize);
	  for (int i=0;i<valuessize;i++) {
	    strObj = PyList_GetItem(values, i);
	    vals[i] = PyFloat_AsDouble(strObj);
	  }
	  int dimensionsize = PyList_Size(dime);
	  std::vector<int> dimension(dimensionsize);
	  for (int i=0;i<dimensionsize;i++) {
	    strObj = PyList_GetItem(dime, i);
	    dimension[i] = PyInt_AsLong(strObj);
	  }

	  foam_setAndWriteCellTensorData(std::string(name), std::string(dataname),vals,dimension);
	  return Py_BuildValue("");
	  
	}
	catch (Foam::error& fErr)
	  {
	    PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	    return NULL;
	  }
	catch (std::exception& e) {
	  PyErr_SetString(PyExc_RuntimeError,e.what());
	  return NULL;
	}
	catch (...) {
	  PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
	  return NULL;
	}
      }
  }


 static PyObject* setCellVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"sisO!",&name,&label,&dataname,&PyList_Type,&values)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      PyObject * strObj;
      int valuessize = PyList_Size(values);

      std::vector<double> vals(valuessize);
      for (int i=0;i<valuessize;i++) {
	strObj = PyList_GetItem(values, i);
	vals[i] = PyFloat_AsDouble(strObj);
      }

 
      foam_setCellVectorData(std::string(name), label, std::string(dataname), vals);
      return Py_BuildValue("");
	    	 
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  
  }

 static PyObject* setCellTensorData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"sisO!",&name,&label,&dataname,&PyList_Type,&values)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      PyObject * strObj;
      int valuessize = PyList_Size(values);

      std::vector<double> vals(valuessize);
      for (int i=0;i<valuessize;i++) {
	strObj = PyList_GetItem(values, i);
	vals[i] = PyFloat_AsDouble(strObj);
      }

 
      foam_setCellTensorData(std::string(name), label, std::string(dataname), vals);
      return Py_BuildValue("");
	    	 
    }
    catch (Foam::error& fErr)
      {
	PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
	return NULL;
      }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  
  }


   static PyObject* writeCellData(PyObject *self, PyObject *args)
  {
     char *name;
     char *dataname;

     if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_writeCellData(std::string(name),std::string(dataname));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 
 
   static PyObject* writeCellVectorData(PyObject *self, PyObject *args)
  {
     char *name;
     char *dataname;

     if (!PyArg_ParseTuple(args,"ss",&name,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_writeCellVectorData(std::string(name),std::string(dataname));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 


  static PyObject* getFaceCount(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      int n = foam_getFaceCount(name);
      return Py_BuildValue("i",n);
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* getPointCount(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      int n = foam_getPointCount(name);
      return Py_BuildValue("i",n);
	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }
 

  static PyObject* getCellCount(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      int n = foam_getCellCount(name);
      return Py_BuildValue("i",n);
	
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* createDefaultFields(PyObject *self, PyObject *args)
  {

    char *name;
    char *solver;
    int io;

    if (!PyArg_ParseTuple(args,"ssi",&name,&solver,&io)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_createDefaultFields(std::string(name),std::string(solver),(io==1));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* modifyNumerics(PyObject *self, PyObject *args)
  {
    char *name;
    char *fvSch;
    char *fvSol;
    char *cD;
    char *TP;
    int io;

    if (!PyArg_ParseTuple(args,"sssssi",&name,&fvSch,&fvSol,&cD,&TP,&io)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_modifyNumerics(std::string(name),std::string(fvSch),std::string(fvSol),std::string(cD),std::string(TP),(io==1));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* modifyDictionary(PyObject *self, PyObject *args)
  {
    char *name;
    char *dictName;
    char *dictContent; 

    if (!PyArg_ParseTuple(args,"sss",&name,&dictName,&dictContent)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_modifyDictionary(std::string(name),std::string(dictName),std::string(dictContent));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  static PyObject* modifyUniformVectorField(PyObject *self, PyObject *args)
  {
    char *name;
    char *fieldName;
    PyObject *value;

    if (!PyArg_ParseTuple(args,"ssO!",&name,&fieldName,&PyList_Type,&value)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      PyObject * strObj;
      int valuesize = PyList_Size(value);

      std::vector<double> val(valuesize);
      for (int i=0;i<valuesize;i++) {
	strObj = PyList_GetItem(value, i);
	val[i] = PyFloat_AsDouble(strObj);
      }


      foam_modifyUniformDimensionedVectorField(std::string(name),std::string(fieldName),val);
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }


  
  static PyObject* extendTensorToBoundaries(PyObject *self, PyObject *args)
  {
    char *name;
    char *fieldname;
 

    if (!PyArg_ParseTuple(args,"ss",&name,&fieldname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_extendToBoundaries(std::string(name),std::string(fieldname));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }

  

  static PyObject* setBC(PyObject *self, PyObject *args)
  {
    char *name;
    char *fieldname;
    char *dict;
 

    if (!PyArg_ParseTuple(args,"sss",&name,&fieldname,&dict)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_setBC(std::string(name),std::string(fieldname),std::string(dict));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  }


  static PyObject* run(PyObject *self, PyObject *args)
  {   
     char *name;
     int nproc;
     double time;
     char *solver;

     if (!PyArg_ParseTuple(args,"sis",&name,&nproc,&solver)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      time = foam_run(std::string(name),nproc,std::string(solver));
      return Py_BuildValue("d",time);
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 

  static PyObject* writeNumerics(PyObject *self, PyObject *args)
  {   
     char *name;

     if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_writeNumerics(std::string(name));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 

  static PyObject* writeDictionary(PyObject *self, PyObject *args)
  {   
     char *name;
     char *dictionaryName;
     int constant;
     if (!PyArg_ParseTuple(args,"ssi",&name,&dictionaryName,&constant)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_writeDictionary(std::string(name),std::string(dictionaryName),(constant==1));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 


  static PyObject* writePathDictionary(PyObject *self, PyObject *args)
  {   
     char *path;
     char *dictionaryName;
     char *head;
     char *dictionaryContent;
     if (!PyArg_ParseTuple(args,"ssss",&path,&dictionaryName,&head,&dictionaryContent)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      
      foam_writePathDictionary(std::string(path),std::string(dictionaryName),std::string(head),std::string(dictionaryContent));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 


  static PyObject* writeFields(PyObject *self, PyObject *args)
  {   
     char *name;

     if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_writeFields(std::string(name));
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 


  static PyObject* updateData(PyObject *self, PyObject *args)
  {   
     char *name;
     double time;

     if (!PyArg_ParseTuple(args,"sd",&name,&time)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_updateData(std::string(name),time);
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }
  } 

  static PyObject* updateTime(PyObject *self, PyObject *args)
  {   
     char *name;
     double time;

     if (!PyArg_ParseTuple(args,"sd",&name,&time)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_updateTime(std::string(name),time);
      return Py_BuildValue("");
    }
    catch (Foam::error& fErr)
    {
      PyErr_SetString(PyExc_RuntimeError,fErr.message().c_str());
      return NULL;
    }
    catch (std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError,e.what());
      return NULL;
    }
    catch (...) {
      PyErr_SetString(PyExc_RuntimeError,"Unknown exception");
      return NULL;
    }

  } 



}
  

  static PyMethodDef FoamMethods[] = {
    {"init",init,METH_VARARGS,"Init wrapper"},
    {"init_IO",init_IO,METH_VARARGS,"Init IO wrapper"},
    {"readMesh",readMesh,METH_VARARGS,"Read foam mesh"},
    {"getPointCoordinates",getPointCoordinates,METH_VARARGS,"Get points coordinates"},
    {"getAllPointCoordinates",getAllPointCoordinates,METH_VARARGS,"Get every points coordinates"},
    {"getFacePoints",getFacePoints,METH_VARARGS,"Get face points"},
    {"getAllFacePoints",getAllFacePoints,METH_VARARGS,"Get every face points"},
    {"getCellPoints",getCellPoints,METH_VARARGS,"Get cell points"},
    {"getAllCellPoints",getAllCellPoints,METH_VARARGS,"Get every cell points"},
    {"getCellDataNames",getCellDataNames,METH_VARARGS,"Get names of data associated to cell"},
    {"getCellVectorDataNames",getCellVectorDataNames,METH_VARARGS,"Get names of vector data associated to cell"},
    {"getCellTensorDataNames",getCellTensorDataNames,METH_VARARGS,"Get names of tensor data associated to cell"},
    {"getCellData",getCellData,METH_VARARGS,"Get data associated to cell"},
    {"getAllCellData",getAllCellData,METH_VARARGS,"Get data associated to cells"},
    {"getCellVectorData",getCellVectorData,METH_VARARGS,"Get vector data associated to cell"},
    {"getCellTensorData",getCellTensorData,METH_VARARGS,"Get tensor data associated to cell"},
    {"getAllCellVectorData",getAllCellVectorData,METH_VARARGS,"Get vector data associated to cells"},
    {"getAllCellTensorData",getAllCellTensorData,METH_VARARGS,"Get tensor data associated to cells"},
    {"getBoundaryCells",getBoundaryCells,METH_VARARGS,"Get specified boundary's cell next to boundary face"},
     {"getBoundaryPatchNames",getBoundaryPatchNames,METH_VARARGS,"Get names of the mesh boundary patches"},
    {"getBoundaryPatchFaces",getBoundaryPatchFaces,METH_VARARGS,"Get mesh boundary patches faces"},
    {"updateCellData",updateCellData,METH_VARARGS,"Update (read) cell data"},
    {"updateCellVectorData",updateCellVectorData,METH_VARARGS,"Update (read) cell vector data"},
    {"addMesh",addMesh,METH_VARARGS,"Ads mesh to wrapper"},
    {"deleteMesh",deleteMesh,METH_VARARGS,"Deletes mesh"},
    {"writeMesh",writeMesh,METH_VARARGS,"Writes mesh to disk"},

    {"setCellData",setCellData,METH_VARARGS,"Sets data associated to cell"},
    {"setAllCellData",setAllCellData,METH_VARARGS,"Sets data associated to all cells"},
    {"setCellVectorData",setCellVectorData,METH_VARARGS,"Sets vector data associated to cell"},
    {"setCellTensorData",setCellTensorData,METH_VARARGS,"Sets tensor data associated to cell"},
    {"setAllCellVectorData",setAllCellVectorData,METH_VARARGS,"Sets vector data associated to cell"},
    {"setAllCellTensorData",setAllCellTensorData,METH_VARARGS,"Sets tensor data associated to cell"},
    {"writeCellData",writeCellData,METH_VARARGS,"Writes cell data to disk"},
    {"writeCellVectorData",writeCellVectorData,METH_VARARGS,"Writes cell vector data to disk"},
    {"getFaceCount",getFaceCount,METH_VARARGS,"Get mesh faces count"},
    {"getPointCount",getPointCount,METH_VARARGS,"Get mesh points count"},
    {"getCellCount",getCellCount,METH_VARARGS,"Get mesh cells count"},
    {"createDefaultFields",createDefaultFields,METH_VARARGS,"Create default fields depending on the solver"},
    {"modifyNumerics",modifyNumerics,METH_VARARGS,"Modify numerical dictionaries through memory"},
    {"modifyDictionary",modifyDictionary,METH_VARARGS,"Modify dictionary through memory"},
    {"modifyUniformVectorField",modifyUniformVectorField,METH_VARARGS,"Modify uniform vector field value through memory"},
    {"writeNumerics",writeNumerics,METH_VARARGS,"Write numerical dictionaries to disk"},
    {"writeDictionary",writeDictionary,METH_VARARGS,"Write dictionary to disk to system or constant directory"},
    {"writePathDictionary",writePathDictionary,METH_VARARGS,"Write given dictionary content to disk to path directory"},
    {"writeFields",writeFields,METH_VARARGS,"Write data fields to disk"},
    {"setBC",setBC,METH_VARARGS,"Modify boundary condition of the specified field through memory"},
    {"extendTensorToBoundaries",extendTensorToBoundaries,METH_VARARGS,"Extend tensor values to boundary face values"},
    {"run",run,METH_VARARGS,"Run the required time steps"},
    {"updateData",updateData,METH_VARARGS,"Update mesh data for given time from disk"},
    {"updateTime",updateTime,METH_VARARGS,"Update mesh time to given time without data reading"},
    
    {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  
  
  PyMODINIT_FUNC initsimphonyfoaminterface(void) {
    Py_InitModule("simphonyfoaminterface",FoamMethods);
  }

