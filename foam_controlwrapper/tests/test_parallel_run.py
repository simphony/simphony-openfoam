import unittest
import os
import re
import shutil
import tempfile

from simphony.api import CUDS, Simulation
from simphony.core.cuba import CUBA
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

from foam_controlwrapper.blockmesh_utils import create_quad_mesh


class WrapperRunTestCase(unittest.TestCase):
    def setUp(self):

        case_name = "simplemesh_parallel"
        mesh_name = "simplemesh_parallel_mesh"
        cuds = CUDS(name=case_name)
        # physics model
        cfd = api.Cfd(name='default model')
        cuds.add(cfd)

        self.sp = api.SolverParameter(name='solver_parameters')
        self.sp._data[CUBA.NUMBER_OF_CORES] = 4

        cuds.add(self.sp)

        sim_time = api.IntegrationTime(name='simulation_time',
                                       current=0.0,
                                       final=1.0,
                                       size=1.0)
        cuds.add(sim_time)

        mat = api.Material(name='a_material')
        mat._data[CUBA.DENSITY] = 1.0
        mat._data[CUBA.DYNAMIC_VISCOSITY] = 1.0
        cuds.add(mat)

        vel_inlet = api.Dirichlet(name='vel_inlet')
        vel_inlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
        vel_inlet._data[CUBA.VELOCITY] = (0.1, 0, 0)
        pres_inlet = api.Neumann(name='pres_inlet')
        pres_inlet._data[CUBA.VARIABLE] = CUBA.PRESSURE

        vel_outlet = api.Neumann(name='vel_outlet')
        vel_outlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
        pres_outlet = api.Dirichlet(name='pres_outlet')
        pres_outlet._data[CUBA.VARIABLE] = CUBA.PRESSURE
        pres_outlet._data[CUBA.PRESSURE] = 0.0

        vel_walls = api.Dirichlet(name='vel_walls')
        vel_walls._data[CUBA.VARIABLE] = CUBA.VELOCITY
        vel_walls._data[CUBA.VELOCITY] = (0, 0, 0)
        pres_walls = api.Neumann(name='pres_walls')
        pres_walls._data[CUBA.VARIABLE] = CUBA.PRESSURE

        vel_frontAndBack = api.Empty(name='vel_frontAndBack')
        vel_frontAndBack._data[CUBA.VARIABLE] = CUBA.VELOCITY
        pres_frontAndBack = api.Empty(name='pres_frontAndBack')
        pres_frontAndBack._data[CUBA.VARIABLE] = CUBA.PRESSURE

        inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet])
        walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls])
        outlet = api.Boundary(name='outlet', condition=[vel_outlet,
                                                        pres_outlet])
        frontAndBack = api.Boundary(name='frontAndBack',
                                    condition=[vel_frontAndBack,
                                               pres_frontAndBack])

        cuds.add(inlet)
        cuds.add(walls)
        cuds.add(outlet)
        cuds.add(frontAndBack)

        corner_points = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0),
                         (5.0, 5.0, 0.0), (0.0, 5.0, 0.0),
                         (0.0, 0.0, 1.0), (5.0, 0.0, 1.0),
                         (5.0, 5.0, 1.0), (0.0, 5.0, 1.0)]
        self.mesh_path = tempfile.mkdtemp()
        mesh = create_quad_mesh(self.mesh_path, mesh_name,
                                corner_points, 5, 5, 5)
        cuds.add(mesh)
        self.cuds = cuds
        self.sim = Simulation(cuds, 'OpenFOAM',
                              engine_interface=EngineInterface.FileIO)
        self.mesh_in_cuds = self.cuds.get(mesh_name)

    def tearDown(self):
        if os.path.exists(self.mesh_in_cuds.path):
            shutil.rmtree(self.mesh_in_cuds.path)
        if os.path.exists(self.mesh_path):
            shutil.rmtree(self.mesh_path)

    def test_parallel_run(self):
        """Test parallel running of OpenFoam

        """

        self.sim.run()

        self.assertEqual(self.sp._data[CUBA.NUMBER_OF_CORES],
                         len([d for d in
                              os.listdir(self.mesh_in_cuds.path)
                              if re.match(r'processor*', d)]))


if __name__ == '__main__':
    unittest.main()
