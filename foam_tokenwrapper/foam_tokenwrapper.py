""" foam_tokenwrapper module

Wrapper module for OpenFOAM tokens using native c++ interface
  
"""
from simphony.cuds.mesh import Mesh, Point, Face, Edge, Cell
import foamtokeninterface

class FoamTokenWrapper(object):

    def __init__(self, model):
        self.model = model

    def initArgs(self, argv):
        """ Sets OpenFoam options for reading
        """
        foamtokeninterface.init(argv)

    def readMesh(self):
        """ Reads OpenFoam mesh based on given options in initArgs
        """
        foamtokeninterface.readMesh()

    def getInternalScalarPointValues(self, name):
        """ Reads internal point values of scalar variable "name"
        """
        return foamtokeninterface.getInternalScalarPointValues(name)

    def getInternalScalarCellValues(self, name):
        """ Reads internal cell values of scalar variable "name"
        """
        return foamtokeninterface.getInternalScalarCellValues(name)

    def getInternalVectorPointValues(self, name):
        """ Reads internal point values of vector variable "name".  
        """
        return foamtokeninterface.getInternalVectorPointValues(name)

    def getInternalVectorCellValues(self, name):
        """ Reads internal cell values of vector variable "name".  
        """
        return foamtokeninterface.getInternalVectorCellValues(name)

    
