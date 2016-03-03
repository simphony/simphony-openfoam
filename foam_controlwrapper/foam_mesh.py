""" Mesh module

This module contains the implementation to store, access,
and modify a mesh and related data

"""
import uuid
import tempfile
import os

from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.mesh import Point, Face, Cell
from simphony.core.cuds_item import CUDSItem

import simphony.core.data_container as dc
import simphonyfoaminterface as foamface

from foam_internalwrapper.mesh_utils import set_cells_data
from foam_internalwrapper.mesh_utils import create_dummy_celldata
from .foam_variables import (dataNameMap, dataKeyMap)
from .foam_variables import dataTypeMap
from foam_internalwrapper.foam_dicts import (dictionaryMaps, parse_map)


class FoamMesh(ABCMesh):
    """ Proxy class to communicate with OpenFoam engines mesh data.

    Parameters
    ----------
    name : str
        name of mesh

    BC : dictionary
       boundary conditions

    mesh : ABCMesh
       mesh to store

    Attributes
    ----------
    name : str
        name of mesh
    data : DataContainer
        DataContainer to store mesh global data (midpoint, ...).
    _time : int
        current time step
    _uuidToFoamLabel : dictionary
        Mapping from uuid to OpenFoam label number
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
        self._time = 0
        self._uuidToFoamLabel = {}
        self._foamCellLabelToUuid = {}
        self._foamFaceLabelToUuid = {}
        self._foamEdgeLabelToUuid = {}
        self._foamPointLabelToUuid = {}
        self._boundaries = {}

        if path:
            self.path = path
        else:
            self.path = os.path.join(tempfile.mkdtemp(), name)
        if mesh:
            # generate uuid mapping
            label = 0
            pointCoordinates = []
            pointMap = {}
            for point in mesh.iter_points():
                pointMap[point.uid] = label
                uid = self._generate_uuid()
                self._uuidToFoamLabel[uid] = label
                self._foamPointLabelToUuid[label] = uid
                for coord in point.coordinates:
                    pointCoordinates.append(coord)
                label += 1

            label = 0
            for edge in mesh.iter_edges():
                uid = self._generate_uuid()
                self._uuidToFoamLabel[uid] = label
                self._foamEdgeLabelToUuid[label] = uid
                label += 1

            label = 0
            facePoints = []
            faceMap = {}
            for face in mesh.iter_faces():
                faceMap[face.uid] = label
                uid = self._generate_uuid()
                self._uuidToFoamLabel[uid] = label
                self._foamFaceLabelToUuid[label] = uid
                # make compressed list of faces points
                facePoints.append(len(face.points))
                for puid in face.points:
                    facePoints.append(pointMap[puid])
                label += 1

            label = 0
            cellPoints = []
            for cell in mesh.iter_cells():
                uid = self._generate_uuid()
                self._uuidToFoamLabel[uid] = label
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
                    if BC[first_key][patchName] == "empty":
                        patchTypes.append("empty")
                    else:
                        patchTypes.append("patch")
                else:
                    patchTypes.append("patch")

            mapContent = dictionaryMaps[solver]
            controlDict = parse_map(mapContent['controlDict'])

            # init objectRegistry and map to mesh name
            foamface.init_IO(name, os.path.abspath(
                os.path.join(self.path, os.pardir)), controlDict)

            # add mesh to objectRegisty
            foamface.addMesh(name, pointCoordinates, cellPoints, facePoints,
                             patchNames, patchFaces, patchTypes)

            foamface.createDefaultFields(name, solver, True)
            # write mesh to disk
            foamface.writeMesh(name)

            # write possible cell data to time directory
            self.copy_cells(mesh.iter_cells())
            foamface.writeFields(name)

    def get_point(self, uuid):
        """Returns a point with a given uuid.

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

        try:
            coords = foamface.getPointCoordinates(self.name,
                                                  self._uuidToFoamLabel[uuid])
            point = Point(coords, uuid)
            return point
        except KeyError:
            error_str = "Trying to get an non-existing point with uuid: {}"
            raise ValueError(error_str.format(uuid))

    def get_edge(self, uuid):
        """Returns an edge with a given uuid.

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

    def get_face(self, uuid):
        """Returns a face with a given uuid.

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

        try:
            pointLabels = foamface.getFacePoints(self.name,
                                                 self._uuidToFoamLabel[uuid])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]

            face = Face(puids, uuid)
            return face
        except KeyError:
            error_str = "Trying to get an non-existing edge with uuid: {}"
            raise ValueError(error_str.format(uuid))

    def get_cell(self, uuid):
        """Returns a cell with a given uuid.

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

        try:

            pointLabels = foamface.getCellPoints(self.name,
                                                 self._uuidToFoamLabel[uuid])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            cell = Cell(puids, uuid)
            label = self._uuidToFoamLabel[uuid]

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

        except KeyError:
            error_str = "Trying to get an non-existing cell with uuid: {}"
            raise ValueError(error_str.format(uuid))

    def add_points(self, points):
        message = 'Points addition not supported yet'
        raise NotImplementedError(message)

    def add_edges(self, edges):
        message = 'Edges addition not supported yet'
        raise NotImplementedError(message)

    def add_faces(self, faces):
        message = 'Faces addition not supported yet'
        raise NotImplementedError(message)

    def add_cells(self, cells):
        message = 'Cells addition not supported yet'
        raise NotImplementedError(message)

    def update_points(self, points):
        message = 'Points update not supported yet'
        raise NotImplementedError(message)

    def update_edges(self, edges):
        message = 'Edges update not supported yet'
        raise NotImplementedError(message)

    def update_faces(self, faces):
        message = 'Faces update not supported yet'
        raise NotImplementedError(message)

    def update_cells(self, cells):
        """ Updates the information of a set of cells.

        Gets the mesh cell identified by the same
        uuid as the provided cell and updates its data
        with the one provided with the new cell. Cell points
        are not updated.

        Parameters
        ----------
        cells : iterable of Cell
            Cell set to be updated

        Raises
        ------
        KeyError
            If the cell was not found in the mesh

        TypeError
            If the object provided is not an iterable of Cell objects

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
                    error_str = "Data named "+data.name+" not supported"
                    raise NotImplementedError(error_str)
                dataName = dataNameMap[data]
                if dataName not in dataNameKeyMap:
                    dataNameKeyMap[dataName] = data
                if dataName not in dataNames and dataName not in newDataNames:
                    newDataNames.append(dataName)
        for dataName in newDataNames:
            create_dummy_celldata(self.name, dataName, True)

        for cell in cellList:
            if cell.uid not in self._uuidToFoamLabel:
                error_str = "Trying to update a non-existing cell with uuid: "\
                    + str(cell.uid)
                raise KeyError(error_str)

            # if points are changed raise warning
            pointLabels = \
                foamface.getCellPoints(self.name,
                                       self._uuidToFoamLabel[cell.uid])
            puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]
            if set(puids) != set(cell.points):
                raise Warning("Cell points can't be updated")

        set_cells_data(self.name, cellList, dataNameKeyMap, True)

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

        TypeError
            If the object provided is not an iterable of Cell objects

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
                    error_str = "Data named "+data.name+" not supported"
                    raise NotImplementedError(error_str)
                dataName = dataNameMap[data]
                if dataName not in dataNameKeyMap:
                    dataNameKeyMap[dataName] = data
                if dataName not in dataNames and dataName not in newDataNames:
                    newDataNames.append(dataName)
        for dataName in newDataNames:
            create_dummy_celldata(self.name, dataName, True)

        set_cells_data(self.name, cellList, dataNameKeyMap, True)

    def iter_points(self, point_uuids=None):
        """Returns an iterator over the selected points.

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
                yield self.get_point(self._foamPointLabelToUuid[label])
        else:
            for uid in point_uuids:
                point = self.get_point(uid)
                yield point

    def iter_edges(self, edge_uuids=None):
        """Return empty list while edges are not supported yet
        Needs to return empty list to get for example H5CUDS add_mesh
        method working with FoamMesh

        """

        return []

    def iter_faces(self, face_uuids=None):
        """Returns an iterator over the selected faces.

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
                face = self.get_face(
                    self._foamFaceLabelToUuid[label])
                yield face
        else:
            for uid in face_uuids:
                face = self.get_face(uid)
                yield face

    def iter_cells(self, cell_uuids=None):
        """Returns an iterator over the selected cells.

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
                yield self.get_cell(self._foamCellLabelToUuid[label])
        else:
            for uid in cell_uuids:
                cell = self.get_cell(uid)
                yield cell

    def has_edges(self):
        """Return false while edges are not supported yet

        """

        return False

    def has_faces(self):
        """Check if the mesh has faces

        Returns
        -------
        bool
            True of there are faces inside the mesh,
            False otherwise

        """
        numberFaces = foamface.getFaceCount(self.name)
        return numberFaces > 0

    def has_cells(self):
        """Check if the mesh has cells

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
        item_type : CUDSItem
            The CUDSItem enum of the type of the items to return the count of.

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

        if item_type == CUDSItem.POINT:
            return foamface.getPointCount(self.name)
        elif item_type == CUDSItem.EDGE:
            error_str = 'Item type {} not supported'
            raise ValueError(error_str.format(item_type))
        elif item_type == CUDSItem.FACE:
            return foamface.getFaceCount(self.name)
        elif item_type == CUDSItem.CELL:
            return foamface.getCellCount(self.name)
        else:
            error_str = 'Item type {} not supported'
            raise ValueError(error_str.format(item_type))

    def write(self):
        """Writes mesh from OpenFOAM's objectRegistry to disk

        """

        foamface.writeMesh(self.name)

    def generate_uuidmapping(self, nPoints, nEdges, nFaces, nCells):
        """Generate uuid mapping assuming continuous ordering of mesh objects

        """

        for point in range(nPoints):
            uuid = self._generate_uuid()
            self._uuidToFoamLabel[uuid] = point
            self._foamPointLabelToUuid[point] = uuid

        for edge in range(nEdges):
            uuid = self._generate_uuid()
            self._uuidToFoamLabel[uuid] = edge
            self._foamEdgeLabelToUuid[edge] = uuid

        for face in range(nFaces):
            uuid = self._generate_uuid()
            self._uuidToFoamLabel[uuid] = face
            self._foamFaceLabelToUuid[face] = uuid

        for cell in range(nCells):
            uuid = self._generate_uuid()
            self._uuidToFoamLabel[uuid] = cell
            self._foamCellLabelToUuid[cell] = uuid

    def _generate_uuid(self):
        """Provides an uuid for the object

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
