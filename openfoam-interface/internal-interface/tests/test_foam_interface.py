""" test foam native interface

This module contains the unitary tests for the
OpenFoam native interface functionalities

"""

import unittest
import os
import shutil
from collections import OrderedDict

from simphony.cuds.mesh import Mesh, Face, Point, Cell
from simphony.core.cuba import CUBA

import simphonyfoaminterface as foamface

from foam_internalwrapper import foam_dicts


class FoamInterfaceTestCase(unittest.TestCase):
    """Test case for OpenFOAM native interface"""

    def setUp(self):

        print "Running FOAM interface test"

        self.name = 'mesh1'
        self.test_dir = 'test'
        self.path = os.path.join(self.test_dir, self.name)
        self.time = '0'
        self.mesh = Mesh(self.name)

        self.points = [
            Point(
                (0.0, 0.0, 0.0)),
            Point(
                (1.0, 0.0, 0.0)),
            Point(
                (1.0, 1.0, 0.0)),
            Point(
                (0.0, 1.0, 0.0)),
            Point(
                (0.0, 0.0, 1.0)),
            Point(
                (1.0, 0.0, 1.0)),
            Point(
                (1.0, 1.0, 1.0)),
            Point(
                (0.0, 1.0, 1.0))
        ]

        puids = self.mesh.add(self.points)

        self.faces = [
            Face([puids[0], puids[3], puids[7], puids[4]]),
            Face([puids[1], puids[2], puids[6], puids[5]]),
            Face([puids[0], puids[1], puids[5], puids[4]]),
            Face([puids[3], puids[2], puids[6], puids[7]]),
            Face([puids[0], puids[1], puids[2], puids[3]]),
            Face([puids[4], puids[5], puids[6], puids[7]])
        ]

        self.mesh.add(self.faces)

        self.cells = [
            Cell(puids)
        ]

        self.puids = puids

        self.mesh.add(self.cells)

        self._uuidToFoamLabelAndType = {}
        self._foamCellLabelToUuid = {}
        self._foamFaceLabelToUuid = {}
        self._foamEdgeLabelToUuid = {}
        self._foamPointLabelToUuid = {}
        # generate uuid mapping
        label = 0
        for point in self.mesh.iter(item_type=CUBA.POINT):
            uid = point.uid
            self._uuidToFoamLabelAndType[uid] = (label, CUBA.POINT)
            self._foamPointLabelToUuid[label] = uid
            label += 1

        label = 0
        for edge in self.mesh.iter(item_type=CUBA.EDGE):
            uid = edge.uid
            self._uuidToFoamLabelAndType[uid] = (label, CUBA.EDGE)
            self._foamEdgeLabelToUuid[label] = uid
            label += 1

        label = 0
        for face in self.mesh.iter(item_type=CUBA.FACE):
            uid = face.uid
            self._uuidToFoamLabelAndType[uid] = (label, CUBA.FACE)
            self._foamFaceLabelToUuid[label] = uid
            label += 1

        label = 0
        for cell in self.mesh.iter(item_type=CUBA.CELL):
            uid = cell.uid
            self._uuidToFoamLabelAndType[uid] = (label, CUBA.CELL)
            self._foamCellLabelToUuid[label] = uid
            label += 1

        # find out boundary patches
        patchNameFacesMap = OrderedDict()
        self.facePoints = []
        i = 0
        for face in self.mesh.iter(item_type=CUBA.FACE):
            boundary = 'boundary' + str(i)
            i += 1
            if boundary not in patchNameFacesMap:
                patchNameFacesMap[boundary] = []
            patchNameFacesMap[boundary].append(
                self._uuidToFoamLabelAndType[face.uid][0])
            # make compressed list of faces points
            self.facePoints.append(len(face.points))
            for puid in face.points:
                self.facePoints.append(self._uuidToFoamLabelAndType[puid][0])

        # make points coordinate list
        self.pointCoordinates = []
        for point in self.mesh.iter(item_type=CUBA.POINT):
            for coord in point.coordinates:
                self.pointCoordinates.append(coord)

        # make compressed list of cells points
        self.cellPoints = []
        for cell in self.mesh.iter(item_type=CUBA.CELL):
            self.cellPoints.append(len(cell.points))
            for puid in cell.points:
                self.cellPoints.append(self._uuidToFoamLabelAndType[puid][0])

        # make patch information
        self.patchNames = []
        self.patchFaces = []
        for patchName in patchNameFacesMap:
            self.patchNames.append(patchName)
            self.patchFaces.append(len(patchNameFacesMap[patchName]))
            for face in patchNameFacesMap[patchName]:
                self.patchFaces.append(face)

        if not self.patchNames:
            error_str = 'Could not initialize with mesh  {}. '
            error_str += 'Mesh has not boundary face definitions.'
            raise ValueError(error_str.format(self.mesh.name))

        self.patchTypes = []

        for i in range(len(self.patchNames)):
            self.patchTypes.append("patch")

        # this to have controlDict file for mesh definition
        mapContent = foam_dicts.dictionaryMaps['pimpleFoam']
        controlDict = foam_dicts.parse_map(mapContent['controlDict'])

        # init objectRegistry and map to mesh name
        foamface.init(self.name, controlDict)

        # add mesh to objectRegisty
        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        foamface.createDefaultFields(self.name, 'pimpleFoam', True)

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_add_mesh(self):
        """Test addMesh method

        """

        pointCoordinates = foamface.getAllPointCoordinates(self.name)
        self.assertEqual(self.pointCoordinates, pointCoordinates)
        cellPoints = foamface.getAllCellPoints(self.name)
        self.assertEqual(set(self.cellPoints), set(cellPoints))
        facePoints = foamface.getAllFacePoints(self.name)
        self.assertEqual(self.facePoints, facePoints)
        patchNames = foamface.getBoundaryPatchNames(self.name)
        self.assertEqual(self.patchNames, patchNames)
        patchFaces = foamface.getBoundaryPatchFaces(self.name)
        self.assertEqual(self.patchFaces, patchFaces)

    def test_get_point_coordinates(self):
        """Test getPointCoordinates method

        """

        for point in self.points:
            coords = foamface.getPointCoordinates(
                self.name,
                self._uuidToFoamLabelAndType[point.uid][0])
            self.assertEqual(point.coordinates, tuple(coords))

    def test_get_face_points(self):
        """Test getFacePoints method

        """

        for face in self.faces:
            pointLabels = foamface.getFacePoints(
                self.name,
                self._uuidToFoamLabelAndType[face.uid][0])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            self.assertEqual(puids, face.points)

    def test_get_cell_points(self):
        """Test getCellPoints method

        """

        for cell in self.cells:
            pointLabels = foamface.getCellPoints(
                self.name,
                self._uuidToFoamLabelAndType[cell.uid][0])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            self.assertEqual(set(puids), set(cell.points))

    def test_get_boundary_patch_names(self):
        """Test getBoundaryPatchNames and getBoundaryPatchFaces method

        """

        patchNames = foamface.getBoundaryPatchNames(self.name)
        patchFaces = foamface.getBoundaryPatchFaces(self.name)

        self.assertEqual(set(patchNames), set(self.patchNames))
        self.assertEqual(set(patchFaces), set(self.patchFaces))

    def test_get_face_count(self):
        """Test getFaceCount method

        """

        nFaces = foamface.getFaceCount(self.name)
        self.assertEqual(nFaces, len(self.faces))

    def test_get_point_count(self):
        """Test getPointCount method

        """

        nPoints = foamface.getPointCount(self.name)
        self.assertEqual(nPoints, len(self.points))

    def test_get_cell_count(self):
        """Test getCellCount method

        """

        nCells = foamface.getCellCount(self.name)
        self.assertEqual(nCells, len(self.cells))


if __name__ == '__main__':
    unittest.main()
