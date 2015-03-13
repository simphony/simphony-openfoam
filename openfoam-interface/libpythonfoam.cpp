#define WM_DP 1
#define NoRepository 1
#include <sstream>
#include "Python.h"
#include "simphonyInterface.H"


extern "C" {

  static PyObject* init(PyObject *self, PyObject *args)
  {   
     char *name;
     char *path;
     if (!PyArg_ParseTuple(args,"ss",&name,&path)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      foam_init(std::string(name),std::string(path));
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


  static PyObject* getPointDataNames(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<std::string> values = foam_getPointDataNames(std::string(name));
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

  static PyObject* getPointData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;

    if (!PyArg_ParseTuple(args,"sis",&name,&label,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getPointData(std::string(name), label, std::string(dataname));
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


  static PyObject* getFaceDataNames(PyObject *self, PyObject *args)
  {
    char *name;
    
    if (!PyArg_ParseTuple(args,"s",&name)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<std::string> values = foam_getFaceDataNames(std::string(name));
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


  static PyObject* getFaceData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    if (!PyArg_ParseTuple(args,"sis",&name,&label,&dataname)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      std::vector<double> values = foam_getFaceData(std::string(name), label, std::string(dataname));
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
 
    if (!PyArg_ParseTuple(args,"sO!O!O!O!O!",&name,&PyList_Type,&points,&PyList_Type,&cells,&PyList_Type,&faces,&PyList_Type,&patchnames,&PyList_Type,&patchfaces)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }    

   try{
     int lenpoints = PyList_Size(points);
     int lencells = PyList_Size(cells);
     int lenfaces = PyList_Size(faces);
     int lenpatchnames = PyList_Size(patchnames);
     int lenpatchfaces = PyList_Size(patchfaces);
     
     std::vector<double> foampoints(lenpoints);
     std::vector<int> foamcells(lencells);
     std::vector<int> foamfaces(lenfaces);
     std::vector<std::string> foampatchnames(lenpatchnames);
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
     
         
     foam_addMesh(std::string(name),foampoints, foamcells, foamfaces, foampatchnames, foampatchfaces);
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


 static PyObject* setPointCoordinates(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    PyObject *coordinates;
    if (!PyArg_ParseTuple(args,"siO!",&name,&label,&PyList_Type,&coordinates)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      PyObject * strObj;
      int valuessize = PyList_Size(coordinates);

      std::vector<double> vals(valuessize);
      for (int i=0;i<valuessize;i++) {
	strObj = PyList_GetItem(coordinates, i);
	vals[i] = PyFloat_AsDouble(strObj);
      }

 
      foam_setPointCoordinates(std::string(name), label, vals);
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


 static PyObject* setFacePoints(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    PyObject *points;

    if (!PyArg_ParseTuple(args,"siO!",&name,&label,&PyList_Type,&points)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      PyObject * strObj;
      int pointssize = PyList_Size(points);
      std::vector<int> pts(pointssize);

      for (int i=0;i<pointssize;i++) {
	strObj = PyList_GetItem(points, i);
	pts[i] = PyInt_AsLong(strObj);
      }
      foam_setFacePoints(std::string(name), label, pts);
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



 static PyObject* setCellPoints(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    PyObject *points;

    if (!PyArg_ParseTuple(args,"siO!",&name,&label,&PyList_Type,&points)) {
      PyErr_SetString(PyExc_RuntimeError,"Invalid arguments");
      return NULL;
    }
    try {
      PyObject * strObj;
      int pointssize = PyList_Size(points);
      std::vector<int> pts(pointssize);

      for (int i=0;i<pointssize;i++) {
	strObj = PyList_GetItem(points, i);
	pts[i] = PyInt_AsLong(strObj);
      }
      foam_setCellPoints(std::string(name), label, pts);
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
    char *dataname;
    PyObject *dimensionset;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"ssO!O!",&name,&dataname,&PyList_Type,&dimensionset,&PyList_Type,&values)) {
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
      
      int setsize = PyList_Size(dimensionset);

      std::vector<int> dims(setsize);
      for (int i=0;i<setsize;i++) {
	strObj = PyList_GetItem(dimensionset, i);
	dims[i] = PyInt_AsLong(strObj);
      }
  
      foam_setCellData(std::string(name), std::string(dataname), dims, vals);
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


 static PyObject* setAllCellVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;
    PyObject *dimensionset;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"ssO!O!",&name,&dataname,&PyList_Type,&dimensionset,&PyList_Type,&values)) {
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
      
      int setsize = PyList_Size(dimensionset);

      std::vector<int> dims(setsize);
      for (int i=0;i<setsize;i++) {
	strObj = PyList_GetItem(dimensionset, i);
	dims[i] = PyInt_AsLong(strObj);
      }


      foam_setCellVectorData(std::string(name), std::string(dataname),dims,vals);
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


 static PyObject* setPointData(PyObject *self, PyObject *args)
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

      foam_setPointData(std::string(name), label, std::string(dataname), value);
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

 static PyObject* setAllPointData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;
    PyObject *dimensionset;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"ssO!O!",&name,&dataname,&PyList_Type,&dimensionset,&PyList_Type,&values)) {
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
      
      int setsize = PyList_Size(dimensionset);

      std::vector<int> dims(setsize);
      for (int i=0;i<setsize;i++) {
	strObj = PyList_GetItem(dimensionset, i);
	dims[i] = PyInt_AsLong(strObj);
      }
  
      foam_setPointData(std::string(name), std::string(dataname), dims, vals);
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


 static PyObject* setAllPointVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;
    PyObject *dimensionset;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"ssO!O!",&name,&dataname,&PyList_Type,&dimensionset,&PyList_Type,&values)) {
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
      
      int setsize = PyList_Size(dimensionset);

      std::vector<int> dims(setsize);
      for (int i=0;i<setsize;i++) {
	strObj = PyList_GetItem(dimensionset, i);
	dims[i] = PyInt_AsLong(strObj);
      }


      foam_setPointVectorData(std::string(name), std::string(dataname),dims,vals);
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

 static PyObject* setPointVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"sisO!",&name,&label,&dataname,&values,&PyList_Type)) {
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

 
      foam_setPointVectorData(std::string(name), label, std::string(dataname), vals);
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


 static PyObject* setFaceData(PyObject *self, PyObject *args)
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

      foam_setFaceData(std::string(name), label, std::string(dataname), value);
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

 static PyObject* setAllFaceData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;
    PyObject *dimensionset;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"ssO!O!",&name,&dataname,&PyList_Type,&dimensionset,&PyList_Type,&values)) {
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
      
      int setsize = PyList_Size(dimensionset);

      std::vector<int> dims(setsize);
      for (int i=0;i<setsize;i++) {
	strObj = PyList_GetItem(dimensionset, i);
	dims[i] = PyInt_AsLong(strObj);
      }
  
      foam_setFaceData(std::string(name), std::string(dataname), dims, vals);
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


 static PyObject* setAllFaceVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    char *dataname;
    PyObject *dimensionset;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"ssO!O!",&name,&dataname,&PyList_Type,&dimensionset,&PyList_Type,&values)) {
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
      
      int setsize = PyList_Size(dimensionset);

      std::vector<int> dims(setsize);
      for (int i=0;i<setsize;i++) {
	strObj = PyList_GetItem(dimensionset, i);
	dims[i] = PyInt_AsLong(strObj);
      }


      foam_setFaceVectorData(std::string(name), std::string(dataname),dims,vals);
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

 static PyObject* setFaceVectorData(PyObject *self, PyObject *args)
  {
    char *name;
    int label;
    char *dataname;
    PyObject *values;
    if (!PyArg_ParseTuple(args,"sisO!",&name,&label,&dataname,&values,&PyList_Type)) {
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

 
      foam_setFaceVectorData(std::string(name), label, std::string(dataname), vals);
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
 

  
  static PyMethodDef FoamMethods[] = {
    {"init",init,METH_VARARGS,"Init wrapper"},
    {"readMesh",readMesh,METH_VARARGS,"Read foam mesh"},
    {"getPointCoordinates",getPointCoordinates,METH_VARARGS,"Get points coordinates"},
    {"getAllPointCoordinates",getAllPointCoordinates,METH_VARARGS,"Get every points coordinates"},
    {"getPointDataNames",getPointDataNames,METH_VARARGS,"Get names of data associated to point"},
    {"getPointData",getPointData,METH_VARARGS,"Get data associated to point"},
    {"getFacePoints",getFacePoints,METH_VARARGS,"Get face points"},
    {"getAllFacePoints",getAllFacePoints,METH_VARARGS,"Get every face points"},
    {"getFaceDataNames",getFaceDataNames,METH_VARARGS,"Get names of data associated to face"},
    {"getFaceData",getFaceData,METH_VARARGS,"Get data associated to face"},
    {"getCellPoints",getCellPoints,METH_VARARGS,"Get cell points"},
    {"getAllCellPoints",getAllCellPoints,METH_VARARGS,"Get every cell points"},
    {"getCellDataNames",getCellDataNames,METH_VARARGS,"Get names of data associated to cell"},
    {"getCellVectorDataNames",getCellVectorDataNames,METH_VARARGS,"Get names of vector data associated to cell"},
    {"getCellData",getCellData,METH_VARARGS,"Get data associated to cell"},
    {"getAllCellData",getAllCellData,METH_VARARGS,"Get data associated to cells"},
    {"getCellVectorData",getCellVectorData,METH_VARARGS,"Get vector data associated to cell"},
    {"getAllCellVectorData",getAllCellVectorData,METH_VARARGS,"Get vector data associated to cells"},
    {"getBoundaryPatchNames",getBoundaryPatchNames,METH_VARARGS,"Get names of the mesh boundary patches"},
    {"getBoundaryPatchFaces",getBoundaryPatchFaces,METH_VARARGS,"Get mesh boundary patches faces"},
    {"updateCellData",updateCellData,METH_VARARGS,"Update (read) cell data"},
    {"updateCellVectorData",updateCellVectorData,METH_VARARGS,"Update (read) cell vector data"},
    {"addMesh",addMesh,METH_VARARGS,"Ads mesh to wrapper"},
    {"deleteMesh",deleteMesh,METH_VARARGS,"Deletes mesh"},
    {"writeMesh",writeMesh,METH_VARARGS,"Writes mesh to disk"},
    {"setPointCoordinates",setPointCoordinates,METH_VARARGS,"Sets points coordinates"},
    {"setFacePoints",setFacePoints,METH_VARARGS,"Set face points"},
    {"setCellPoints",setCellPoints,METH_VARARGS,"Sets cell points"},
    {"setCellData",setCellData,METH_VARARGS,"Sets data associated to cell"},
    {"setAllCellData",setAllCellData,METH_VARARGS,"Sets data associated to all cells"},
    {"setCellVectorData",setCellVectorData,METH_VARARGS,"Sets vector data associated to cell"},
    {"setAllCellVectorData",setAllCellVectorData,METH_VARARGS,"Sets vector data associated to cell"},
    {"setPointData",setPointData,METH_VARARGS,"Sets data associated to point"},
    {"setAllPointData",setAllPointData,METH_VARARGS,"Sets data associated to all points"},
    {"setPointVectorData",setPointVectorData,METH_VARARGS,"Sets vector data associated to point"},
    {"setAllPointVectorData",setAllPointVectorData,METH_VARARGS,"Sets vector data associated to point"},
    {"setFaceData",setFaceData,METH_VARARGS,"Sets data associated to face"},
    {"setAllFaceData",setAllFaceData,METH_VARARGS,"Sets data associated to all faces"},
    {"setFaceVectorData",setFaceVectorData,METH_VARARGS,"Sets vector data associated to face"},
    {"setAllFaceVectorData",setAllFaceVectorData,METH_VARARGS,"Sets vector data associated to face"},
    {"writeCellData",writeCellData,METH_VARARGS,"Writes cell data to disk"},
    {"writeCellVectorData",writeCellVectorData,METH_VARARGS,"Writes cell vector data to disk"},
    {"getFaceCount",getFaceCount,METH_VARARGS,"Get mesh faces count"},
    {"getPointCount",getPointCount,METH_VARARGS,"Get mesh points count"},
    {"getCellCount",getCellCount,METH_VARARGS,"Get mesh cells count"},

    {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  
  
  PyMODINIT_FUNC initsimphonyfoaminterface(void) {
    Py_InitModule("simphonyfoaminterface",FoamMethods);
  }
  
}
