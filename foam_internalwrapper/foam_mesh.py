""" Mesh module

This module contains the implementation to store, access,
and modify a mesh and related data

"""
import uuid
import tempfile
import os


from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.mesh import Point, Face, Cell
from simphony.core.cuba import CUBA

import simphony.core.data_container as dc
import simphonyfoaminterface as foamface

from .foam_dicts import (get_dictionary_maps, parse_map, check_boundary_names)
from foam_controlwrapper.foam_variables import dataNameMap
from foam_controlwrapper.foam_variables import (dataKeyMap, dataTypeMap)
from .mesh_utils import (create_dummy_celldata, set_cells_data)


class FoamMesh(ABCMesh):
    """ Proxy class to communicate with OpenFoam engines mesh data.

    Parameters
    ----------
    name : str
       name of mesh

    BC : dictionary
       boundary conditions

    solver : str
       name of OpenFoam solver

    mesh : ABCMesh
       mesh to store

    Attributes
    ----------
    name : str
        name of mesh
    data : DataContainer
        DataContainer to store mesh global data (midpoint, ...).
    _uuidToFoamLabelAndType : dictionary
        Mapping from uuid to OpenFoam label number and type
    _foamCellLabelToUuid : dictionary
        Mapping from OpenFoam cell label number to uuid
    _foamFaceLabelToUuid : dictionary
        Mapping from  OpenFoam face label number to uuid
    _foamEdgeLabelToUuid : dictionary
        Mapping from  OpenFoam edge label number to uuid
    _foamPointLabelToUuid : dictionary
        Mapping from  OpenFoam point label number to uuid
    _boundaries : dictionary
        Mapping from boundary name to list of boundary faces

    """

    def __init__(self, name, BC, solver, mesh=None, path=None):
        super(FoamMesh, self).__init__()
        self.name = name
        self.data = dc.DataContainer()
        self._uuidToFoamLabelAndType = {}
        self._foamCellLabelToUuid = {}
        self._foamFaceLabelToUuid = {}
        self._foamEdgeLabelToUuid = {}
        self._foamPointLabelToUuid = {}
        self._boundaries = {}
        self._time = 0.0

        if path:
            self.path = path
        else:
            self.path = os.path.join(tempfile.mkdtemp(), name)
        if mesh:

            # generate uuid mapping
            label = 0
            pointCoordinates = []
            pointMap = {}
            for point in mesh.iter(item_type=CUBA.POINT):
                pointMap[point.uid] = label
                uid = self._generate_uuid()
                self._uuidToFoamLabelAndType[uid] = (label, CUBA.POINT)
                self._foamPointLabelToUuid[label] = uid
                for coord in point.coordinates:
                    pointCoordinates.append(coord)
                label += 1

            label = 0
            for edge in mesh.iter(item_type=CUBA.EDGE):
                uid = self._generate_uuid()
                self._uuidToFoamLabelAndType[uid] = (label, CUBA.EDGE)
                self._foamEdgeLabelToUuid[label] = uid
                label += 1

            label = 0
            facePoints = []
            faceMap = {}
            for face in mesh.iter(item_type=CUBA.FACE):
                faceMap[face.uid] = label
                uid = self._generate_uuid()
                self._uuidToFoamLabelAndType[uid] = (label, CUBA.FACE)
                self._foamFaceLabelToUuid[label] = uid
                # make compressed list of faces points
                facePoints.append(len(face.points))
                for puid in face.points:
                    facePoints.append(pointMap[puid])
                label += 1

            label = 0
            cellPoints = []
            for cell in mesh.iter(item_type=CUBA.CELL):
                uid = self._generate_uuid()
                self._uuidToFoamLabelAndType[uid] = (label, CUBA.CELL)
                self._foamCellLabelToUuid[label] = uid
                cellPoints.append(len(cell.points))
                for puid in cell.points:
                    cellPoints.append(pointMap[puid])
                label += 1

            pointMap.clear()

            # make patch information
            patchNames = []
            patchFaces = []
            if hasattr(mesh, '_boundaries'):
                for patchName in mesh._boundaries:
                    patchNames.append(patchName)
                    self._boundaries[patchName] = []
                    patchFaces.append(len(mesh._boundaries[patchName]))
                    for fuid in mesh._boundaries[patchName]:
                        flabel = faceMap[fuid]
                        new_fuid = self._foamFaceLabelToUuid[flabel]
                        self._boundaries[patchName].append(new_fuid)
                        patchFaces.append(flabel)

            faceMap.clear()

            patchTypes = []
            for patchName in patchNames:
                if BC:
                    first_key = BC.keys()[0]
                    check_boundary_names(BC[first_key].keys(), patchNames,
                                         first_key)
                    if BC[first_key][patchName] == "empty":
                        patchTypes.append("empty")
                    else:
                        patchTypes.append("patch")
                else:
                    patchTypes.append("patch")

            mapContent = get_dictionary_maps(solver, False)
            controlDict = parse_map(mapContent['controlDict'])

            # init objectRegistry and map to mesh name
            foamface.init(name, controlDict)

            foamface.updateTime(name, self._time)
            # add mesh to objectRegisty
            foamface.addMesh(name, pointCoordinates, cellPoints,
                             facePoints, patchNames, patchFaces, patchTypes)

            foamface.createDefaultFields(name, solver, False)

            # write possible cell data to time directory
            self.copy_cells(mesh.iter(item_type=CUBA.CELL))
            # correct boundary face labels
            patchNames = foamface.getBoundaryPatchNames(name)
            patchFaces = foamface.getBoundaryPatchFaces(name)
            boundaries = {}
            i = 0
            k = 0
            while i < len(patchFaces):
                boundaries[patchNames[k]] = []
                start = i+1
                end = start+patchFaces[i]
                i += 1
                for j in range(start, end):
                    boundaries[patchNames[k]].append(
                        self._foamFaceLabelToUuid[patchFaces[j]])
                    i += 1
                k += 1
            self._boundaries = boundaries

    def _get_point(self, uuid):
        """ Returns a point with a given uuid.

        Returns the point stored in the mesh
        identified by uuid. If such point do not
        exists an exception is raised.

        Parameters
        ----------
        uuid
            uuid of the desired point.

        Returns
        -------
        Point
            Mesh point identified by uuid

        Raises
        ------
        Exception
            If the point identified by uuid was not found

        """

        label, type = self._uuidToFoamLabelAndType[uuid]
        if type != CUBA.POINT:
            raise KeyError("No point with uuid {}".format(uuid))

        coords = foamface.getPointCoordinates(self.name, label)
        point = Point(coords, uuid)
        return point

    def _get_edge(self, uuid):
        """ Returns an edge with a given uuid.

        Returns the edge stored in the mesh
        identified by uuid. If such edge do not
        exists an exception is raised.

        Parameters
        ----------
        uuid
            uuid of the desired edge.

        Returns
        -------
        Edge
            Edge identified by uuid

        Raises
        ------
        Exception
            If the edge identified by uuid was not found

        """
        message = "Edges are not supported yet in OpenFoam engine"
        raise NotImplementedError(message)

    def _get_face(self, uuid):
        """ Returns a face with a given uuid.

        Returns the face stored in the mesh
        identified by uuid. If such face do not
        exists an exception is raised.

        Parameters
        ----------
        uuid
            uuid of the desired face.

        Returns
        -------
        Face
            Face identified by uuid

        Raises
        ------
        Exception
            If the face identified by uuid was not found

        """

        label, type = self._uuidToFoamLabelAndType[uuid]
        if type != CUBA.FACE:
            raise KeyError("No face with uuid {}".format(uuid))

        pointLabels = foamface.getFacePoints(self.name, label)
        puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]

        face = Face(puids, uuid)
        return face

    def get_boundary_cells(self, boundary):
        """Returns boundar cells for a given boundary.

        Returns the face stored in the mesh
        identified by uuid. If such face do not
        exists an exception is raised.

        Parameters
        ----------
        boundary : str
            boundary name

        Returns
        -------
        Iterator Cell
            Iterator to boundary cells

        Raises
        ------
        Exception
            If the face identified by uuid was not found

        """

        cells = foamface.getBoundaryCells(self.name, boundary)
        for label in cells:
            yield self._get_cell(self._foamCellLabelToUuid[label])

    def _get_cell(self, uuid):
        """ Returns a cell with a given uuid.

        Returns the cell stored in the mesh
        identified by uuid . If such cell do not
        exists an exception is raised.

        Parameters
        ----------
        uuid
            uuid of the desired cell.

        Returns
        -------
        Cell
            Cell with id identified by uuid

        Raises
        ------
        Exception
            If the cell identified by uuid was not found

        """

        label, type = self._uuidToFoamLabelAndType[uuid]
        if type != CUBA.CELL:
            raise KeyError("No Cell with uuid {}".format(uuid))

        pointLabels = foamface.getCellPoints(self.name, label)
        puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
        cell = Cell(puids, uuid)
        label, _ = self._uuidToFoamLabelAndType[uuid]

        dataNames = foamface.getCellDataNames(self.name)
        dataNames += foamface.getCellVectorDataNames(self.name)
        dataNames += foamface.getCellTensorDataNames(self.name)
        for dataName in set(dataKeyMap.keys()).intersection(dataNames):
            if dataTypeMap[dataKeyMap[dataName]] == "scalar":
                cell.data[dataKeyMap[dataName]] = \
                    foamface.getCellData(self.name,
                                         label,
                                         dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "vector":
                cell.data[dataKeyMap[dataName]] = \
                    foamface.getCellVectorData(self.name,
                                               label,
                                               dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "tensor":
                cell.data[dataKeyMap[dataName]] = \
                    foamface.getCellTensorData(self.name,
                                               label,
                                               dataName)

        return cell

    def _add_points(self, points):
        message = 'Points addition not supported yet'
        raise NotImplementedError(message)

    def _add_edges(self, edges):
        message = 'Edges addition not supported yet'
        raise NotImplementedError(message)

    def _add_faces(self, face):
        message = 'Faces addition not supported yet'
        raise NotImplementedError(message)

    def _add_cells(self, cells):
        message = 'Cells addition not supported yet'
        raise NotImplementedError(message)

    def _update_points(self, points):
        message = 'Point update not supported yet'
        raise NotImplementedError(message)

    def _update_edges(self, edges):
        message = 'Edges update not supported yet'
        raise NotImplementedError(message)

    def _update_faces(self, faces):
        message = 'Faces update not supported yet'
        raise NotImplementedError(message)

    def _update_cells(self, cells):
        """ Updates the information of a set of cells.

        Gets the mesh cell identified by the same
        uuid as the provided cell and updates its information
        with the one provided with the new cell.

        Parameters
        ----------
        cells : iterable of Cell
            Cell set to be updated

        Raises
        ------
        KeyError
            If the cell was not found in the mesh

        TypeError
            If the object provided is not a cell

        """

        dataNames = foamface.getCellDataNames(self.name)
        dataNames += foamface.getCellVectorDataNames(self.name)
        dataNames += foamface.getCellTensorDataNames(self.name)

        # if cell data does not exists in the mesh at all, initialize it
        newDataNames = []
        dataNameKeyMap = {}
        cellList = list(cells)
        for cell in cellList:
            for data in cell.data:
                if data not in dataNameMap:
                    error_str = "Data named "+data+" not supported"
                    raise NotImplementedError(error_str)
                dataName = dataNameMap[data]
                if dataName not in dataNameKeyMap:
                    dataNameKeyMap[dataName] = data
                if dataName not in dataNames and dataName not in newDataNames:
                    newDataNames.append(dataName)
        for dataName in newDataNames:
            create_dummy_celldata(self.name, dataName)

        for cell in cellList:
            if cell.uid not in self._uuidToFoamLabelAndType:
                error_str = "Trying to update a non-existing cell with uuid: "\
                    + str(cell.uid)
                raise KeyError(error_str)

            # if points are changed raise warning
            pointLabels = \
                foamface.getCellPoints(
                    self.name,
                    self._uuidToFoamLabelAndType[cell.uid][0])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            if set(puids) != set(cell.points):
                raise Warning("Cell points can't be updated")

        set_cells_data(self.name, cellList, dataNameKeyMap)

    def copy_cells(self, cells):
        """ Copy the information of a set of cells.

        Gets the mesh cell and copy its data
        with the one provided with the new cell. Cell points
        are not copied.

        Parameters
        ----------
        cells : iterable of Cell
            Cell set to be updated

        Raises
        ------
        KeyError
            If the cell was not found in the mesh

        TypeError
            If the object provided is not a cell

        """

        dataNames = foamface.getCellDataNames(self.name)
        dataNames += foamface.getCellVectorDataNames(self.name)
        dataNames += foamface.getCellTensorDataNames(self.name)

        # if cell data does not exists in the mesh at all, initialize it
        newDataNames = []
        dataNameKeyMap = {}
        cellList = list(cells)
        for cell in cellList:
            for data in cell.data:
                if data not in dataNameMap:
                    error_str = "Data named "+data+" not supported"
                    raise NotImplementedError(error_str)
                dataName = dataNameMap[data]
                if dataName not in dataNameKeyMap:
                    dataNameKeyMap[dataName] = data
                if dataName not in dataNames and dataName not in newDataNames:
                    newDataNames.append(dataName)
        for dataName in newDataNames:
            create_dummy_celldata(self.name, dataName)

        set_cells_data(self.name, cellList, dataNameKeyMap)

    def _iter_points(self, point_uuids=None):
        """ Returns an iterator over the selected points.

        Returns an iterator over the points with uuid in
        point_ids. If none of the ids in point_ids exists,
        an empty iterator is returned. If there is no ids
        inside point_ids, a iterator over all points of
        the mesh is returned instead.

        Parameters
        ----------
        point_uuids : list of uuids, optional
            uuids of the desired points, default empty

        Returns
        -------
        iter
            Iterator over the selected points

        """

        if point_uuids is None:
            pointCount = foamface.getPointCount(self.name)
            for label in range(pointCount):
                point = self._get_point(self._foamPointLabelToUuid[label])
                yield Point.from_point(point)
        else:
            for uid in point_uuids:
                point = self._get_point(uid)
                yield Point.from_point(point)

    def _iter_edges(self, edge_uuids=None):
        """ Return empty list while edges are not supported yet

        """

        return []

    def _iter_faces(self, face_uuids=None):
        """ Returns an iterator over the selected faces.

        Returns an iterator over the faces with uuid in
        face_uuids. If none of the uuids in face_uuids exists,
        an empty iterator is returned. If there is no uuids
        inside face_uuids, a iterator over all faces of
        the mesh is returned instead.

        Parameters
        ----------
        face_uuids : list of uuids, optional
            Uuids of the desired faces, default empty

        Returns
        -------
        iter
            Iterator over the selected faces

        """

        if face_uuids is None:
            faceCount = foamface.getFaceCount(self.name)
            for label in range(faceCount):
                face = self._get_face(self._foamFaceLabelToUuid[label])
                yield Face.from_face(face)
        else:
            for uid in face_uuids:
                face = self._get_face(uid)
                yield Face.from_face(face)

    def _iter_cells(self, cell_uuids=None):
        """ Returns an iterator over the selected cells.

        Returns an iterator over the cells with uuid in
        cell_uuids. If none of the uuids in cell_uuids exists,
        an empty iterator is returned. If there is no uuids
        inside cell_uuids, a iterator over all cells of
        the mesh is returned instead.

        Parameters
        ----------
        cell_uuids : list of uuids, optional
            Uuids of the desired cell, default empty

        Returns
        -------
        iter
            Iterator over the selected cells

        """

        if cell_uuids is None:
            cellCount = foamface.getCellCount(self.name)
            for label in range(cellCount):
                yield self._get_cell(self._foamCellLabelToUuid[label])
        else:
            for uid in cell_uuids:
                cell = self._get_cell(uid)
                yield cell

    def _has_points(self):
        """Returns False if there are no points.
        """
        return foamface.getPointCount(self.name) > 0

    def _has_edges(self):
        """ Return false while edges are not supported yet

        """

        return False

    def _has_faces(self):
        """ Check if the mesh has faces

        Returns
        -------
        bool
            True of there are faces inside the mesh,
            False otherwise

        """
        numberFaces = foamface.getFaceCount(self.name)
        return numberFaces > 0

    def _has_cells(self):
        """ Check if the mesh has cells

        Returns
        -------
        bool
            True of there are cells inside the mesh,
            False otherwise

        """
        numberCells = foamface.getCellCount(self.name)
        return numberCells > 0

    def count_of(self, item_type):
        """ Return the count of item_type in the container.

        Parameters
        ----------
        item_type : CUBA item
            The CUBA enum of the type of the items to return the count of.

        Returns
        -------
        count : int
            The number of items of item_type in the container.

        Raises
        ------
        ValueError :
            If the type of the item is not supported in the current
            container.

        """

        if item_type == CUBA.POINT:
            return foamface.getPointCount(self.name)
        elif item_type == CUBA.EDGE:
            return 0
        elif item_type == CUBA.FACE:
            return foamface.getFaceCount(self.name)
        elif item_type == CUBA.CELL:
            return foamface.getCellCount(self.name)
        else:
            error_str = 'Item type {} not supported'
            raise ValueError(error_str.format(item_type))

    def generate_uuidmapping(self, nPoints, nEdges, nFaces, nCells):
        '''generate uuid mapping assuming continuous ordering of mesh objects

        '''

        for point in range(nPoints):
            uuid = self._generate_uuid()
            self._uuidToFoamLabelAndType[uuid] = (point, CUBA.POINT)
            self._foamPointLabelToUuid[point] = uuid

        for edge in range(nEdges):
            uuid = self._generate_uuid()
            self._uuidToFoamLabelAndType[uuid] = (edge, CUBA.EDGE)
            self._foamEdgeLabelToUuid[edge] = uuid

        for face in range(nFaces):
            uuid = self._generate_uuid()
            self._uuidToFoamLabelAndType[uuid] = (face, CUBA.FACE)
            self._foamFaceLabelToUuid[face] = uuid

        for cell in range(nCells):
            uuid = self._generate_uuid()
            self._uuidToFoamLabelAndType[uuid] = (cell, CUBA.CELL)
            self._foamCellLabelToUuid[cell] = uuid

    def _generate_uuid(self):
        """ Provides an uuid for the object

        Provides san uuid as defined in the standard RFC 4122
        """

        return uuid.uuid4()

    def add_boundaries(self, boundaries):
        """adds boundaries to boundary mapping

        """

        for bname, blist in boundaries:
            if bname not in self._boundaries:
                self.boundaries[bname] = blist

    def get_boundaries(self):
        """ get boundaries

        """

        return self._boundaries
