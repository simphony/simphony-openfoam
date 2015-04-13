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
    {"getFacePoints",getFacePoints,METH_VARARGS,"Get face points"},
    {"getAllFacePoints",getAllFacePoints,METH_VARARGS,"Get every face points"},
    {"getCellPoints",getCellPoints,METH_VARARGS,"Get cell points"},
    {"getAllCellPoints",getAllCellPoints,METH_VARARGS,"Get every cell points"},
    {"getBoundaryPatchNames",getBoundaryPatchNames,METH_VARARGS,"Get names of the mesh boundary patches"},
    {"getBoundaryPatchFaces",getBoundaryPatchFaces,METH_VARARGS,"Get mesh boundary patches faces"},
    {"addMesh",addMesh,METH_VARARGS,"Ads mesh to wrapper"},
    {"deleteMesh",deleteMesh,METH_VARARGS,"Deletes mesh"},
    {"writeMesh",writeMesh,METH_VARARGS,"Writes mesh to disk"},
    {"getFaceCount",getFaceCount,METH_VARARGS,"Get mesh faces count"},
    {"getPointCount",getPointCount,METH_VARARGS,"Get mesh points count"},
    {"getCellCount",getCellCount,METH_VARARGS,"Get mesh cells count"},

    {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  
  
  PyMODINIT_FUNC initsimphonyfoaminterface(void) {
    Py_InitModule("simphonyfoaminterface",FoamMethods);
  }
  
}
