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
from simphony.core.data_container import DataContainer
import simphonyfoaminterfaceII as foamface

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

        puids = [self.mesh.add_point(point) for point in self.points]

        self.faces = [
            Face([puids[0], puids[3], puids[7], puids[4]],
                 data=DataContainer({CUBA.LABEL: 0})),
            Face([puids[1], puids[2], puids[6], puids[5]],
                 data=DataContainer({CUBA.LABEL: 1})),
            Face([puids[0], puids[1], puids[5], puids[4]],
                 data=DataContainer({CUBA.LABEL: 2})),
            Face([puids[3], puids[2], puids[6], puids[7]],
                 data=DataContainer({CUBA.LABEL: 3})),
            Face([puids[0], puids[1], puids[2], puids[3]],
                 data=DataContainer({CUBA.LABEL: 4})),
            Face([puids[4], puids[5], puids[6], puids[7]],
                 data=DataContainer({CUBA.LABEL: 5}))

        ]

        [self.mesh.add_face(face) for face in self.faces]

        self.cells = [
            Cell(puids)
        ]

        self.puids = puids

        [self.mesh.add_cell(cell) for cell in self.cells]

        self._uuidToFoamLabel = {}
        self._foamCellLabelToUuid = {}
        self._foamFaceLabelToUuid = {}
        self._foamEdgeLabelToUuid = {}
        self._foamPointLabelToUuid = {}
        # generate uuid mapping
        label = 0
        for point in self.mesh.iter_points():
            uid = point.uid
            self._uuidToFoamLabel[uid] = label
            self._foamPointLabelToUuid[label] = uid
            label += 1

        label = 0
        for edge in self.mesh.iter_edges():
            uid = edge.uid
            self._uuidToFoamLabel[uid] = label
            self._foamEdgeLabelToUuid[label] = uid
            label += 1

        label = 0
        for face in self.mesh.iter_faces():
            uid = face.uid
            self._uuidToFoamLabel[uid] = label
            self._foamFaceLabelToUuid[label] = uid
            label += 1

        label = 0
        for cell in self.mesh.iter_cells():
            uid = cell.uid
            self._uuidToFoamLabel[uid] = label
            self._foamCellLabelToUuid[label] = uid
            label += 1

        # find out boundary patches
        patchNameFacesMap = OrderedDict()
        self.facePoints = []
        for face in self.mesh.iter_faces():
            if CUBA.LABEL in face.data:
                boundary = 'boundary' + str(face.data[CUBA.LABEL])
                if boundary not in patchNameFacesMap:
                    patchNameFacesMap[boundary] = []
                patchNameFacesMap[boundary].append(
                    self._uuidToFoamLabel[face.uid])

            # make compressed list of faces points
            self.facePoints.append(len(face.points))
            for puid in face.points:
                self.facePoints.append(self._uuidToFoamLabel[puid])

        # make points coordinate list
        self.pointCoordinates = []
        for point in self.mesh.iter_points():
            for coord in point.coordinates:
                self.pointCoordinates.append(coord)

        # make compressed list of cells points
        self.cellPoints = []
        for cell in self.mesh.iter_cells():
            self.cellPoints.append(len(cell.points))
            for puid in cell.points:
                self.cellPoints.append(self._uuidToFoamLabel[puid])

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
        foamface.init(self.name, os.path.abspath(os.path.join(self.path,
                      os.pardir)), controlDict)

        # add mesh to objectRegisty
        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        foamface.createDefaultFields(self.name, 'pimpleFoam')

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_add_mesh(self):
        """Test addMesh method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

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

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        for point in self.points:
            coords = foamface.getPointCoordinates(
                self.name,
                self._uuidToFoamLabel[point.uid])
            self.assertEqual(point.coordinates, tuple(coords))

    def test_get_face_points(self):
        """Test getFacePoints method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        for face in self.faces:
            pointLabels = foamface.getFacePoints(
                self.name,
                self._uuidToFoamLabel[face.uid])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            self.assertEqual(puids, face.points)

    def test_get_cell_points(self):
        """Test getCellPoints method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        for cell in self.cells:
            pointLabels = foamface.getCellPoints(
                self.name,
                self._uuidToFoamLabel[cell.uid])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            self.assertEqual(set(puids), set(cell.points))

    def test_get_boundary_patch_names(self):
        """Test getBoundaryPatchNames and getBoundaryPatchFaces method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        patchNames = foamface.getBoundaryPatchNames(self.name)
        patchFaces = foamface.getBoundaryPatchFaces(self.name)
        i = 0
        k = 0
        facePatchMap = {}
        while i < len(patchFaces):
            start = i+1
            end = start+patchFaces[i]
            i += 1
            for j in range(start, end):
                # here we assume that the boundaries are named
                # as boundary0,...
                # this to overcome limitation in tableextensio.pyx
                # at a moment
                facePatchMap[patchFaces[j]] = \
                    int(patchNames[k].replace('boundary', ''))
                i += 1
            k += 1

        for face in self.faces:
            self.assertEqual(
                face.data[CUBA.LABEL],
                facePatchMap[self._uuidToFoamLabel[face.uid]])

    def test_get_face_count(self):
        """Test getFaceCount method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        nFaces = foamface.getFaceCount(self.name)
        self.assertEqual(nFaces, len(self.faces))

    def test_get_point_count(self):
        """Test getPointCount method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        nPoints = foamface.getPointCount(self.name)
        self.assertEqual(nPoints, len(self.points))

    def test_get_cell_count(self):
        """Test getCellCount method

        """

        foamface.addMesh(self.name,
                         self.pointCoordinates,
                         self.cellPoints,
                         self.facePoints,
                         self.patchNames,
                         self.patchFaces,
                         self.patchTypes)

        nCells = foamface.getCellCount(self.name)
        self.assertEqual(nCells, len(self.cells))

if __name__ == '__main__':
    unittest.main()
