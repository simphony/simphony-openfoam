"""Example to solve 2D poiseuille flow

"""

#from __future__ import print_function

from simphony.engine import openfoam_file_io

from simphony.core.cuba import CUBA

#from mayavi.scripts import mayavi2

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
#from simphony.engine import openfoam_file_io
from simphony.engine import EngineInterface

import tempfile

case_name = 'poiseuille'
mesh_name = 'poiseuille_mesh'

cuds = CUDS(name=case_name)

mat = api.Material(name='a_material')
mat.data[CUBA.DENSITY] = 1.0
mat.data[CUBA.DYNAMIC_VISCOSITY] = 1.0
cuds.add(mat)

rheo = api.NewtonianFluidModel(name='newtonian')
cuds.add(rheo)

compr = api.IncompressibleFluidModel(name='incompressible')
cuds.add(compr)

turb = api.LaminarFlowModel(name='laminar')
cuds.add(turb)

sim_time = api.IntegrationTime(name='simulation_time',
                               current = 0.0,
                               final = 1.0,
                               size = 0.01)
cuds.add(sim_time)

#sp = api.SolverParameter(name = 'solver_parameters')


# create computational mesh

corner_points = [(0.0, 0.0, 0.0), (20.0e-3, 0.0, 0.0),
                 (20.0e-3, 1.0e-3, 0.0), (0.0, 1.0e-3, 0.0),
                 (0.0, 0.0, 0.1e-3), (20.0e-3, 0.0, 0.1e-3),
                 (20.0e-3, 1.0e-3, 0.1e-3), (0.0, 1.0e-3, 0.1e-3)]
# elements in x -direction
nex = 8
# elements in y -direction
ney = 4
# this routine creates one block quad mesh
#  - boundaries have predescribed names (inlet, walls, outlet, frontAndBack)

mesh_cuds=CUDS(name=mesh_name)
sim = Simulation(mesh_cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)

mesh = sim.create_quad_mesh(tempfile.mkdtemp(), mesh_name,
                                         corner_points, nex, ney, 1)
cuds.add(mesh)


vel_inlet = api.Dirichlet(name='vel_inlet')
vel_inlet.data[CUBA.VELOCITY] = (0.1, 0, 0)

vel_outlet = api.Neumann(name='vel_outlet')
vel_outlet.data[CUBA.VELOCITY] = (0, 0, 0)

vel_walls = api.Dirichlet(name='vel_walls')
vel_walls.data[CUBA.VELOCITY] = (0, 0, 0)

vel_frontAndBack = api.Empty(name='vel_frontAndBack')
vel_frontAndBack.data[CUBA.VELOCITY] = (0, 0, 0)

pres_inlet = api.Neumann(name='pres_inlet')
pres_inlet.data[CUBA.PRESSURE] = 0

pres_outlet = api.Dirichlet(name='pres_outlet')
pres_outlet.data[CUBA.PRESSURE] = 0

pres_walls = api.Neumann(name='pres_walls')
pres_walls.data[CUBA.PRESSURE] = 0

pres_frontAndBack = api.Empty(name='pres_frontAndBack')
pres_frontAndBack.data[CUBA.PRESSURE] = (0, 0, 0)

inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet])
frontAndBack = api.Boundary(name='frontAndBack', condition=[vel_frontAndBack, pres_frontAndBack])
cuds.add(inlet)
cuds.add(walls)
cuds.add(outlet)
cuds.add(frontAndBack)


sim = Simulation(cuds, 'OpenFOAM', engine_interface = EngineInterface.FileIO)


sim.run()

average_pressure = 0.0
mesh_inside_wrapper = cuds.get(mesh_name)
for cell in mesh_inside_wrapper.get_boundary_cells('inlet'):
    average_pressure += cell.data[CUBA.PRESSURE]

average_pressure /= len(mesh_inside_wrapper._boundaries['inlet'])

#print "Average pressure on inlet: ", average_pressure
