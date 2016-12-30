""" Mesh module

This module contains the implementation to store, access,
and modify a mesh and related data

"""
import uuid
import tempfile
import os
import time

from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.mesh import Point, Face, Cell
from simphony.cuds.meta import api
from simphony.core.cuba import CUBA

import simphony.core.data_container as dc
import simphonyfoaminterface as foamface

from foam_internalwrapper.mesh_utils import set_cells_data
from foam_internalwrapper.mesh_utils import create_dummy_celldata
from .foam_variables import (dataNameMap, dataKeyMap, dataTypeMap,
                             dataDimensionMap)
from foam_internalwrapper.foam_dicts import (get_dictionary_maps, parse_map,
                                             check_boundary_names,
                                             get_foam_boundary_condition)


class FoamMesh(ABCMesh):
    """ Proxy class to communicate with OpenFoam engines mesh data.

    Parameters
    ----------
    name : str
        name of mesh

    cuds : CUDS
       CUDS

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

    def __init__(self, name, cuds, solver, mesh=None, path=None):
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
        self._foamMaterialLabelToUuid = None
        if path:
            self.path = path
        else:
            self.path = os.path.join(tempfile.mkdtemp(), name)
        if mesh:
            # generate uuid mapping
            if hasattr(mesh, '_get_packed_coordinate_list'):
                pointCoordinates = mesh._get_packed_coordinate_list()
                i = 0
                label = 0
                while i < len(pointCoordinates):
                    uid = self._generate_uuid()
                    self._uuidToFoamLabel[uid] = label
                    self._foamPointLabelToUuid[label] = uid
                    label += 1
                    i += 3
            else:
                label = 0
                pointCoordinates = []
                pointMap = {}
                for point in mesh._iter_points():
                    pointMap[point.uid] = label
                    uid = self._generate_uuid()
                    self._uuidToFoamLabel[uid] = label
                    self._foamPointLabelToUuid[label] = uid
                    for coord in point.coordinates:
                        pointCoordinates.append(coord)
                    label += 1

            label = 0
            for edge in mesh._iter_edges():
                uid = self._generate_uuid()
                self._uuidToFoamLabel[uid] = label
                self._foamEdgeLabelToUuid[label] = uid
                label += 1

            faceMap = {}
            if hasattr(mesh, '_get_packed_face_list'):
                label = -1
                facePoints = mesh._get_packed_face_list()
                i = 0
                while i < len(facePoints):
                    label += 1
                    n_points = facePoints[i]
                    i += 1 + n_points
                    face_uid = mesh._foamFaceLabelToUuid[label]
                    faceMap[face_uid] = label
                    uid = self._generate_uuid()
                    self._uuidToFoamLabel[uid] = label
                    self._foamFaceLabelToUuid[label] = uid
            else:
                label = 0
                facePoints = []
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

            if hasattr(mesh, '_get_packed_cell_list'):
                cellPoints = mesh._get_packed_cell_list()
                cell_label = -1
                i = 0
                while i < len(cellPoints):
                    cell_label += 1
                    n_points = cellPoints[i]
                    i += 1
                    uid = self._generate_uuid()
                    self._uuidToFoamLabel[uid] = cell_label
                    self._foamCellLabelToUuid[cell_label] = uid
                    i += n_points
            else:
                label = 0
                cellPoints = []
                cellMap = {}
                for cell in mesh.iter_cells():
                    uid = self._generate_uuid()
                    cellMap[cell.uid] = uid
                    self._uuidToFoamLabel[uid] = label
                    self._foamCellLabelToUuid[label] = uid
                    cellPoints.append(len(cell.points))
                    for puid in cell.points:
                        cellPoints.append(pointMap[puid])
                    label += 1
                pointMap.clear()

            if hasattr(mesh, '_get_cell_data_map'):
                cell_data_map = mesh._get_cell_data_map()
                self._foamMaterialLabelToUuid = mesh._foamMaterialLabelToUuid
            else:
                cell_data_map = {}
                nCells = mesh.count_of(CUBA.CELL) 
                for cell in mesh.iter_cells():
                    label = self._uuidToFoamLabel[cellMap[cell.uid]]
                    for key in cell.data:
                        if key not in cell_data_map:
                            if  dataTypeMap[key] == "scalar":
                                cell_data_map[key] = [0] * nCells
                            elif dataTypeMap[key] == "vector":
                                cell_data_map[key] = [0] * (nCells * 3)
                            elif dataTypeMap[key] == "tensor":
                                cell_data_map[key] = [0] * (nCells * 9)
                        if key == CUBA.MATERIAL:
                            self._foamMaterialLabelToUuid = cell.data[key]
                            cell_data_map[key][label] = 0
                        else:
                            if  dataTypeMap[key] == "scalar":
                                cell_data_map[key][label] = cell.data[key]
                            elif dataTypeMap[key] == "vector":
                                for i in range(len(cell.data[key])):
                                    k = 3 * label + i
                                    cell_data_map[key][k] = cell.data[key][i]  
                            elif dataTypeMap[key] == "tensor":
                                for i in range(len(cell.data[key])):
                                    k = 9 * label + i
                                    cell_data_map[key][k] = cell.data[key][i]  

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
            if cuds:
                bcs = {}
                for boundary in cuds.iter(api.Boundary):
                    bcs[boundary.name] = \
                        get_foam_boundary_condition(boundary.condition[0])
                check_boundary_names(bcs.keys(), patchNames)

                for patchName in patchNames:
                    if patchName in bcs and bcs[patchName] == "empty":
                        patchTypes.append("empty")
                    else:
                        patchTypes.append("patch")
            else:
                for patchName in patchNames:
                    patchTypes.append("patch")

            mapContent = get_dictionary_maps(solver, False)
            controlDict = parse_map(mapContent['controlDict'])

            # init objectRegistry and map to mesh name
            foamface.init_IO(name, os.path.abspath(
                os.path.join(self.path, os.pardir)), controlDict)

            # update time
            foamface.updateTime(name, self._time)
            # add mesh to objectRegisty
            foamface.addMesh(name, pointCoordinates, cellPoints, facePoints,
                             patchNames, patchFaces, patchTypes)

            # create default fields
            foamface.createDefaultFields(name, solver, True)
            # write mesh to disk
            foamface.writeMesh(name)

            # copy possible cell data to time register
            self.copy_cells(cell_data_map)
            cell_data_map.clear()
            foamface.writeFields(name)

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

    def _get_edge(self, uuid):
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

    def _get_all_faces(self):
        """Returns all faces at once.

        Returns
        -------
        faces generator
            All faces in mesh as generator
        """

        pointLabels = foamface.getAllFacePoints(self.name)
        face_label = -1
        i = 0
        while i < len(pointLabels):
            face_label += 1
            n_points = pointLabels[i]
            i += 1
            puids = []
            for j in range(n_points):
                puids.append(self._foamPointLabelToUuid[pointLabels[i]])
                i += 1
            yield Face(puids, self._foamFaceLabelToUuid[face_label])

    def _get_face(self, uuid):
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
                    if dataKeyMap[dataName] == CUBA.MATERIAL:
                        cell.data[dataKeyMap[dataName]] = \
                            self._foamMaterialLabelToUuid
                    else:
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

    def _add_points(self, points):
        message = 'Points addition not supported yet'
        raise NotImplementedError(message)

    def _add_edges(self, edges):
        message = 'Edges addition not supported yet'
        raise NotImplementedError(message)

    def _add_faces(self, faces):
        message = 'Faces addition not supported yet'
        raise NotImplementedError(message)

    def _add_cells(self, cells):
        message = 'Cells addition not supported yet'
        raise NotImplementedError(message)

    def _update_points(self, points):
        message = 'Points update not supported yet'
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
        for data_name in newDataNames:
            create_dummy_celldata(self.name, data_name, True)


        if CUBA.MATERIAL in cellList[0].data:
                self._foamMaterialLabelToUuid = cellList[0].data[CUBA.MATERIAL]

        for cell in cellList:
            if CUBA.MATERIAL in cell.data:
                    cell.data[CUBA.MATERIAL] = 0
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

    def copy_cells(self, cell_data_map):
        """ Copy the information of a set of cells.

        Gets the mesh cell and copy its data
        with the one provided with the new cell. Cell points
        are not copied.

        Parameters
        ----------
        cell_data_map : dictionary
             map from data key to cell data values in
             mesh internal order

        """

        dataNames = foamface.getCellDataNames(self.name)
        dataNames += foamface.getCellVectorDataNames(self.name)
        dataNames += foamface.getCellTensorDataNames(self.name)

        # if cell data does not exists in the mesh at all, initialize it
        newDataNames = [] 
        for key in cell_data_map:
            if key not in dataNameMap:
                error_str = "Data named "+data.name+" not supported"
                raise NotImplementedError(error_str)
            dataName = dataNameMap[key]
            if dataName not in dataNames and dataName not in newDataNames:
                newDataNames.append(dataName)
        for dataName in newDataNames:
            create_dummy_celldata(self.name, dataName, True)

        for key, data in cell_data_map.iteritems():
            dimension = dataDimensionMap[key]
            dataName = dataNameMap[key]
            if dataTypeMap[key] == "scalar":
                foamface.setAllCellData(self.name, dataName, 1,
                                        data, dimension)
            elif dataTypeMap[key] == "vector":
                foamface.setAllCellVectorData(self.name, dataName, 1,
                                              data, dimension)
            elif dataTypeMap[key] == "tensor":
                foamface.setAllCellTensorData(self.name, dataName, 1,
                                              data, dimension)

    def _iter_points(self, point_uuids=None):
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
            point_coordinates = foamface.getAllPointCoordinates(self.name)
            i = 0
            label = 0
            while i < len(point_coordinates):
                yield Point(point_coordinates[i:i+3],
                            self._foamPointLabelToUuid[label])
                label += 1
                i += 3
        else:
            for uid in point_uuids:
                point = self._get_point(uid)
                yield point

    def _iter_edges(self, edge_uuids=None):
        """Return empty list while edges are not supported yet
        Needs to return empty list to get for example H5CUDS _add_mesh
        method working with FoamMesh

        """

        return []

    def _iter_faces(self, face_uuids=None):
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
            pointLabels = foamface.getAllFacePoints(self.name)
            face_label = -1
            i = 0
            while i < len(pointLabels):
                face_label += 1
                n_points = pointLabels[i]
                i += 1
                puids = []
                for j in range(n_points):
                    puids.append(self._foamPointLabelToUuid[pointLabels[i]])
                    i += 1
                yield Face(puids, self._foamFaceLabelToUuid[face_label])
        else:
            for uid in face_uuids:
                face = self._get_face(uid)
                yield face

    def _get_cells(self):
        """Returns all cells in the mesh.

        Returns
        -------
        Cells
            list of Cell

        """

        try:

            pointLabels = foamface.getAllCellPoints(self.name)

            cells = []

            i = 0
            cell_label = 0
            while i < len(pointLabels):
                puids = []
                for ip in range(pointLabels[i]):
                    i += 1
                    puids.append(self._foamPointLabelToUuid[pointLabels[i]])
                i += 1
                uuid = self._foamCellLabelToUuid[cell_label]
                cells.append(Cell(puids, uuid))
                cell_label += 1

            dataNames = foamface.getCellDataNames(self.name)
            dataNames += foamface.getCellVectorDataNames(self.name)
            dataNames += foamface.getCellTensorDataNames(self.name)
            for dataName in set(dataKeyMap.keys()).intersection(dataNames):
                if dataTypeMap[dataKeyMap[dataName]] == "scalar":
                    if dataKeyMap[dataName] == CUBA.MATERIAL:
                        for cell in cells:
                            cell.data[dataKeyMap[dataName]] =\
                                self._foamMaterialLabelToUuid
                    dataS = foamface.getAllCellData(self.name, dataName)
                    for cell in cells:
                        cell.data[dataKeyMap[dataName]] =\
                            dataS[self._uuidToFoamLabel[cell.uid]]

                elif dataTypeMap[dataKeyMap[dataName]] == "vector":
                    dataS = foamface.getAllCellVectorData(self.name, dataName)
                    for cell in cells:
                        cell_label = self._uuidToFoamLabel[cell.uid]
                        istart = 3*cell_label
                        cell.data[dataKeyMap[dataName]] =\
                            dataS[istart:istart+3]

                elif dataTypeMap[dataKeyMap[dataName]] == "tensor":
                    dataS = foamface.getAllCellTensorData(self.name, dataName)
                    for cell in cells:
                        cell_label = self._uuidToFoamLabel[cell.uid]
                        istart = 9*cell_label
                        cell.data[dataKeyMap[dataName]] =\
                            dataS[istart:istart+9]

            return cells

        except KeyError:
            error_str = "Trying to get an non-existing cell with uuid: {}"
            raise ValueError(error_str.format(uuid))

    def _iter_cells(self, cell_uuids=None):
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
            pointLabels = foamface.getAllCellPoints(self.name)
            start = time.time()
            data_map = self._get_cell_data_map()
            print "Time spend in data map: ",time.time()-start
            cell_label = -1
            i = 0
            while i < len(pointLabels):
                cell_label += 1
                n_points = pointLabels[i]
                i += 1
                puids = []
                for j in range(n_points):
                    puids.append(self._foamPointLabelToUuid[pointLabels[i]])
                    i += 1
                cell = Cell(puids, self._foamCellLabelToUuid[cell_label])
                for dataKey, data in data_map.iteritems():
                    if dataTypeMap[dataKey] == "scalar":
                        if dataKey == CUBA.MATERIAL:
                            cell.data[dataKey] = \
                                self._foamMaterialLabelToUuid
                        else:
                            cell.data[dataKey] = data[cell_label]
                    elif dataTypeMap[dataKey] == "vector":
                        cell.data[dataKey] = \
                            [data[cell_label * 3 + k] for k in range(3)]
                    elif dataTypeMap[dataKey] == "tensor":
                        cell.data[dataKey] = \
                            [data[cell_label * 9 + k] for k in range(9)]
                yield cell
            data_map.clear()
        else:
            for uid in cell_uuids:
                cell = self._get_cell(uid)
                yield cell

    def _has_points(self):
        """Check if the mesh has points

        Returns
        -------
        bool
            True of there are points in the mesh,
            False otherwise

        """
        numberPoints = foamface.getPointCount(self.name)
        return numberPoints > 0

    def _has_edges(self):
        """Return false while edges are not supported yet

        """

        return False

    def _has_faces(self):
        """Check if the mesh has faces

        Returns
        -------
        bool
            True of there are faces inside the mesh,
            False otherwise

        """
        numberFaces = foamface.getFaceCount(self.name)
        return numberFaces > 0

    def _has_cells(self):
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
        """ Return the count of points, edges, faces or cells in the container.

        Parameters
        ----------
        item_type : CUBA
            The CUBA keywordto give type of the items to return the count of.

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

    def _get_boundaries(self):
        """ get boundaries

        """

        return self._boundaries

    def _get_packed_coordinate_list(self):
        """ get packed list of points coordinate values
        """
        return foamface.getAllPointCoordinates(self.name)

    def _get_packed_face_list(self):
        """ get packed list of faces point labels
        """
        return foamface.getAllFacePoints(self.name)

    def _get_packed_cell_list(self):
        """ get packed list of cells point labels
        """
        return foamface.getAllCellPoints(self.name)

    def _get_cell_data_map(self):
        """ get map for mesh data 
        """

        dataNames = foamface.getCellDataNames(self.name)
        dataNames += foamface.getCellVectorDataNames(self.name)
        dataNames += foamface.getCellTensorDataNames(self.name)
        
        dataMap = {}
        for dataName in  set(dataKeyMap.keys()).intersection(dataNames):
            if dataTypeMap[dataKeyMap[dataName]] == "scalar":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllCellData(self.name, dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "vector":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllCellVectorData(self.name, dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "tensor":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllCellTensorData(self.name, dataName)
        return dataMap
