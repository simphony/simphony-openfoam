""" Mesh module

This module contains the implementation to store, access,
and modify a mesh and related data

"""
import uuid
import tempfile
import os
from multiprocessing import Pool, cpu_count

from simphony.cuds.mesh import Mesh
from simphony.cuds.mesh_items import Point, Face, Cell

from simphony.core.cuba import CUBA
from simphony.cuds.meta.api import PhaseVolumeFraction

import simphonyfoaminterface as foamface

from foam_internalwrapper.mesh_utils import (set_cells_data,
                                             create_dummy_celldata,
                                             get_cells_in_range)

from .foam_variables import (dataNameMap, dataKeyMap, dataTypeMap,
                             dataDimensionMap, phaseNames)
from foam_internalwrapper.foam_dicts import (get_dictionary_maps, parse_map,
                                             check_boundary_names,
                                             get_foam_boundary_condition,
                                             not_empty)


class FoamMesh(Mesh):
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
    _uuidToFoamLabelAndType : dictionary
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
    _foamPhaseNameToMaterial : dictionary
        Mapping from phase name to material

    """

    def __init__(self, name, cuds, solver, mesh=None, path=None):
        super(FoamMesh, self).__init__(name)

        self._time = 0
        self._uuidToFoamLabelAndType = {}
        self._foamCellLabelToUuid = {}
        self._foamFaceLabelToUuid = {}
        self._foamEdgeLabelToUuid = {}
        self._foamPointLabelToUuid = {}
        self._boundaries = {}
        self._foamPhaseNameToMaterial = {}

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
                    self._uuidToFoamLabelAndType[uid] = (label, CUBA.POINT)
                    self._foamPointLabelToUuid[label] = uid
                    label += 1
                    i += 3
            else:
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
                    self._uuidToFoamLabelAndType[uid] = (label, CUBA.FACE)
                    self._foamFaceLabelToUuid[label] = uid
            else:
                label = 0
                facePoints = []
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

            if hasattr(mesh, '_get_packed_cell_list'):
                cellPoints = mesh._get_packed_cell_list()
                cell_label = -1
                i = 0
                while i < len(cellPoints):
                    cell_label += 1
                    n_points = cellPoints[i]
                    i += 1
                    uid = self._generate_uuid()
                    self._uuidToFoamLabelAndType[uid] = (cell_label, CUBA.CELL)
                    self._foamCellLabelToUuid[cell_label] = uid
                    i += n_points
            else:
                label = 0
                cellPoints = []
                cellMap = {}
                for cell in mesh.iter(item_type=CUBA.CELL):
                    uid = self._generate_uuid()
                    cellMap[cell.uid] = uid
                    self._uuidToFoamLabelAndType[uid] = (label, CUBA.CELL)
                    self._foamCellLabelToUuid[label] = uid
                    cellPoints.append(len(cell.points))
                    for puid in cell.points:
                        cellPoints.append(pointMap[puid])
                    label += 1
                pointMap.clear()

            if hasattr(mesh, '_foamPhaseNameToMaterial') and \
                    mesh._foamPhaseNameToMaterial:
                self._foamPhaseNameToMaterial = mesh._foamPhaseNameToMaterial
            elif cuds and not_empty(cuds.iter(item_type=CUBA.MATERIAL)):
                materials = list(cuds.iter(item_type=CUBA.MATERIAL))
                im = 0
                for material in materials:
                    self._foamPhaseNameToMaterial[phaseNames[im]] = material
                    im += 1

            if hasattr(mesh, '_get_cell_data_map'):
                cell_data_map = mesh._get_cell_data_map()
            else:
                cell_data_map = {}
                nCells = mesh.count_of(CUBA.CELL)
                for cell in mesh.iter(item_type=CUBA.CELL):
                    label, _ = self._uuidToFoamLabelAndType[cellMap[cell.uid]]
                    for key in cell.data:
                        if key not in cell_data_map:
                            if dataTypeMap[key] == "scalar":
                                cell_data_map[key] = [0] * nCells
                            elif dataTypeMap[key] == "vector":
                                cell_data_map[key] = [0] * (nCells * 3)
                            elif dataTypeMap[key] == "tensor":
                                cell_data_map[key] = [0] * (nCells * 9)
                        if dataTypeMap[key] == "scalar":
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
                for boundary in cuds.iter(item_type=CUBA.BOUNDARY):
                    bcs[boundary.name] = \
                        get_foam_boundary_condition(
                        boundary.condition[0], self._foamPhaseNameToMaterial,
                        solver)
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
        label, type = self._uuidToFoamLabelAndType[uuid]
        if type != CUBA.POINT:
            raise KeyError("No point with uuid {}".format(uuid))

        coords = foamface.getPointCoordinates(self.name, label)
        point = Point(coords, uuid)
        return point

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

        # Edges are not supported yet in OpenFoam engine
        raise KeyError("Edge not found for uuid {}".format(uuid))

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
        label, type = self._uuidToFoamLabelAndType[uuid]
        if type != CUBA.FACE:
            raise KeyError("No face with uuid {}".format(uuid))

        pointLabels = foamface.getFacePoints(self.name, label)
        puids = [self._foamPointLabelToUuid[lbl] for lbl in pointLabels]

        face = Face(puids, uuid)

        return face

    def get_boundary_cells(self, boundary):
        """Returns boundary cells for a given boundary.

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
        if len(self._foamPhaseNameToMaterial) > 1:
            dataNames.append(dataNameMap[CUBA.VOLUME_FRACTION])
        for dataName in set(dataKeyMap.keys()).intersection(dataNames):
            if dataName == dataNameMap[CUBA.VOLUME_FRACTION]:
                # currently this is only for volume_fraction and for two phase
                dName = dataName + '.' + phaseNames[0]
                vol_frac1 = foamface.getCellData(self.name, label, dName)
                material1 = self._foamPhaseNameToMaterial[phaseNames[0]]
                material2 = self._foamPhaseNameToMaterial[phaseNames[1]]
                phase1_vol_frac = PhaseVolumeFraction(material1, vol_frac1)
                phase2_vol_frac = PhaseVolumeFraction(material2, 1 - vol_frac1)
                cell.data[dataKeyMap[dataName]] = [phase1_vol_frac,
                                                   phase2_vol_frac]

            elif dataTypeMap[dataKeyMap[dataName]] == "scalar":
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

        # check that if volume fraction is in data, phase names are mapped to
        # materials
        if CUBA.VOLUME_FRACTION in dataNameKeyMap.values() \
                and not self._foamPhaseNameToMaterial:
            cell = cellList[0]
            phase_vol_fracs = cell.data[CUBA.VOLUME_FRACTION]
            im = 0
            for phase_vol_frac in phase_vol_fracs:
                self._foamPhaseNameToMaterial[phaseNames[im]] = \
                    phase_vol_frac.material
                im += 1
        set_cells_data(self.name, cellList, dataNameKeyMap,
                       self._foamPhaseNameToMaterial, True)

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
                error_str = "Data named "+key+" not supported"
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
                if key == CUBA.VOLUME_FRACTION:
                    foamface.setAllCellData(self.name,
                                            dataName + '.' + phaseNames[0], 1,
                                            data, dimension)
                else:
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
            data_map = self._get_point_data_map()

            i = 0
            label = 0
            while i < len(point_coordinates):
                point = Point(point_coordinates[i:i+3],
                            self._foamPointLabelToUuid[label])
                for dataKey, data in data_map.iteritems():
                    if dataTypeMap[dataKey] == "scalar":
                        if dataKey == CUBA.VOLUME_FRACTION:
                            if self._foamPhaseNameToMaterial:
                                material1 = self._foamPhaseNameToMaterial[
                                    phaseNames[0]]
                                material2 = self._foamPhaseNameToMaterial[
                                    phaseNames[1]]
                                vol_frac1 = data[label]
                                phase1_vol_frac = PhaseVolumeFraction(
                                    material1, vol_frac1)
                                phase2_vol_frac = PhaseVolumeFraction(
                                    material2, 1 - vol_frac1)
                                point.data[dataKey] = [phase1_vol_frac,
                                                      phase2_vol_frac]
                        else:
                            point.data[dataKey] = data[label]
                    elif dataTypeMap[dataKey] == "vector":
                        point.data[dataKey] = \
                            [data[label * 3 + k] for k in range(3)]
                    elif dataTypeMap[dataKey] == "tensor":
                        point.data[dataKey] = \
                            [data[label * 9 + k] for k in range(9)]
                
                yield point
                label += 1
                i += 3
            data_map.clear()
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

    def _find_cells_puids(self, point_labels):
        clist = []
        i = 0
        while i < len(point_labels):
            n_points = point_labels[i]
            i += 1
            puids = []
            for j in range(n_points):
                puids.append(self._foamPointLabelToUuid[point_labels[i]])
                i += 1
            clist.append(puids)
        return clist

    def _iter_cells_parallel(self):
        """Returns an iterator over all cells.

        Returns an iterator over the cells. Cell instances are
        made parallell.


        Returns
        -------
        iter
            Iterator over cells

        """

        pointLabels = foamface.getAllCellPoints(self.name)
        data_map = self._get_cell_data_map()

        n_jobs = cpu_count()
        cells_puids = self._find_cells_puids(pointLabels)

        n_cells = len(cells_puids)
        group_size = n_cells / n_jobs
        last_group_size = group_size + n_cells % n_jobs
        cell_indeces = []
        for i in range(n_jobs - 1):
            cellis = []
            cellis.append(i * group_size)
            cellis.append((i+1) * group_size - 1)
            cell_indeces.append(cellis)
        cellis = []
        cellis.append((n_jobs - 1) * group_size)
        cellis.append((n_jobs - 1) * group_size + last_group_size - 1)
        cell_indeces.append(cellis)

        pool = Pool(n_jobs)
        args = [[cell_indeces[i][0], cell_indeces[i][1],
                 cells_puids, data_map, self._foamCellLabelToUuid,
                 self._foamPhaseNameToMaterial]
                for i in range(n_jobs)]
        results = pool.map(get_cells_in_range, args)

        for res in results:
            for item in res:
                yield item
        pool.close()
        data_map.clear()

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
            data_map = self._get_cell_data_map()
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
                        if dataKey == CUBA.VOLUME_FRACTION:
                            if self._foamPhaseNameToMaterial:
                                material1 = self._foamPhaseNameToMaterial[
                                    phaseNames[0]]
                                material2 = self._foamPhaseNameToMaterial[
                                    phaseNames[1]]
                                vol_frac1 = data[cell_label]
                                phase1_vol_frac = PhaseVolumeFraction(
                                    material1, vol_frac1)
                                phase2_vol_frac = PhaseVolumeFraction(
                                    material2, 1 - vol_frac1)
                                cell.data[dataKey] = [phase1_vol_frac,
                                                      phase2_vol_frac]
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
        """Returns False if there are no points
        """

        return (foamface.getPointCount(self.name) > 0)

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

    def write(self):
        """Writes mesh from OpenFOAM's objectRegistry to disk

        """

        foamface.writeMesh(self.name)

    def generate_uuidmapping(self, nPoints, nEdges, nFaces, nCells):
        """Generate uuid mapping assuming continuous ordering of mesh objects

        """

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
        # take phase dependance away from the name
        dataNames = [dname.partition('.')[0] for dname in dataNames]

        for dataName in set(dataKeyMap.keys()).intersection(dataNames):
            if dataTypeMap[dataKeyMap[dataName]] == "scalar":
                if dataKeyMap[dataName] == CUBA.VOLUME_FRACTION:
                    dName = dataName + '.' + phaseNames[0]
                    dataMap[dataKeyMap[dataName]] = \
                        foamface.getAllCellData(self.name, dName)
                else:
                    dataMap[dataKeyMap[dataName]] = \
                        foamface.getAllCellData(self.name, dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "vector":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllCellVectorData(self.name, dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "tensor":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllCellTensorData(self.name, dataName)
        return dataMap

    def _get_point_data_map(self):
        """ get map for mesh pointwise data
        """

        dataNames = foamface.getCellDataNames(self.name)
        dataNames += foamface.getCellVectorDataNames(self.name)
        dataNames += foamface.getCellTensorDataNames(self.name)

        dataMap = {}
        # take phase dependance away from the name
        dataNames = [dname.partition('.')[0] for dname in dataNames]

        for dataName in set(dataKeyMap.keys()).intersection(dataNames):
            if dataTypeMap[dataKeyMap[dataName]] == "scalar":
                if dataKeyMap[dataName] == CUBA.VOLUME_FRACTION:
                    dName = dataName + '.' + phaseNames[0]
                    dataMap[dataKeyMap[dataName]] = \
                        foamface.getAllPointData(self.name, dName)
                else:
                    dataMap[dataKeyMap[dataName]] = \
                        foamface.getAllPointData(self.name, dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "vector":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllPointVectorData(self.name, dataName)
            elif dataTypeMap[dataKeyMap[dataName]] == "tensor":
                dataMap[dataKeyMap[dataName]] = \
                    foamface.getAllPointTensorData(self.name, dataName)
        return dataMap
