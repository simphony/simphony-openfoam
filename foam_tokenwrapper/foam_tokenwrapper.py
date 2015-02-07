""" foam_tokenwrapper module

Wrapper module for OpenFOAM tokens using native c++ interface

"""
from simphony.core.cuba import CUBA
from simphony.cuds.mesh import Mesh, Point, Face, Cell
from simphony.core.data_container import DataContainer
import foamtokeninterface


class FoamTokenWrapper(object):

    def __init__(self):
        self.CM = DataContainer()
        self.BC = DataContainer()
        self.SP = DataContainer()

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
        """ Reads internal point values of vector variable "name"
        """
        return foamtokeninterface.getInternalVectorPointValues(name)

    def getInternalVectorCellValues(self, name):
        """ Reads internal cell values of vector variable "name"
        """
        return foamtokeninterface.getInternalVectorCellValues(name)

    def get_point_data(self, simphonyMesh, foamMeshMap, cubaName):
        """
            Get nodal data from OpenFoam
        """

        if cubaName == CUBA.PRESSURE:
            data = self.getInternalScalarPointValues("p")
            for point in simphonyMesh.iter_points():
                point.data.update({CUBA.PRESSURE: data[foamMeshMap[point.uuid]-1]})
                simphonyMesh.update_point(point)
        elif cubaName == CUBA.VELOCITY:
            data = self.getInternalVectorPointValues("U")
            for point in simphonyMesh.iter_points():
                point.data.update({CUBA.VELOCITY: data[foamMeshMap[point.uuid]-1]})
                simphonyMesh.update_point(point)
        else:
            raise KeyError("Variable "+cubaName+" transfer not supported yet")

    def get_mesh(self, caseName, timeStep):
        """ OpenFoam mesh to Simphony mesh
        """
        simphonyMesh = Mesh(caseName)
        argv = ["test", "-case", caseName, "-time", timeStep]

        self.initArgs(argv)
        self.readMesh()

        foamPointCoordinates = foamtokeninterface.getPointCoordinates()

        foamLabelToUuid = {}
        pointUuids = []
        k = 1
        for i in range(0, len(foamPointCoordinates), 3):
            spoint = Point((foamPointCoordinates[i], foamPointCoordinates[i+1], foamPointCoordinates[i+2]))
            simphonyMesh.add_point(spoint)
            foamLabelToUuid[spoint.uuid] = k
            pointUuids.append(spoint.uuid)
            k += 1

        foamFacePoints = foamtokeninterface.getFacePoints()

        i = 0
        while i < len(foamFacePoints):
            k = 1
            points = []
            start = i+1
            end = i+1+foamFacePoints[i]
            i += 1
            for j in range(start, end):
                points.append(pointUuids[foamFacePoints[j]])
                i += 1
            face = Face(points)
            simphonyMesh.add_face(face)
            foamLabelToUuid[face.uuid] = k
            k += 1

        foamCellPoints = foamtokeninterface.getCellPoints()

        i = 0
        while i < len(foamCellPoints):
            k = 0
            points = []
            start = i+1
            end = i+1+foamCellPoints[i]
            i += 1
            for j in range(start, end):
                points.append(pointUuids[foamCellPoints[j]])
                i += 1

            cell = Cell(points)
            simphonyMesh.add_cell(cell)
            foamLabelToUuid[cell.uuid] = k
            k += 1

        del pointUuids

        return (simphonyMesh, foamLabelToUuid)
