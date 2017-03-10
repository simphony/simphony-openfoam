"""Example to solve 2D vortex shedding flow

"""

from mayavi.scripts import mayavi2

import foam_controlwrapper as wrapper
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

import vortex_shedding_mesh

import tempfile
import math

case_name = 'vortex_shedding'
mesh_name = 'vortex_shedding_mesh'


# create cuds
cuds = CUDS(name=case_name)

# physics model
cfd = api.Cfd(name='default model')

# these are already by default set in CFD
cfd.rheology_model = api.NewtonianFluidModel(name='newtonian')
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')
# add to cuds
cuds.add([cfd])

# material
mat = api.Material(name='a_material')
mat._data[CUBA.DENSITY] = 1.0
mat._data[CUBA.DYNAMIC_VISCOSITY] = 2.0e-2
cuds.add([mat])

# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=150.0,
                               size=0.025)
cuds.add([sim_time])

# create mesh
mesh = wrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                 vortex_shedding_mesh.blockMeshDict)
cuds.add([mesh])

vel_b = [None]*6
pres_b = [None]*6
vel_b[0] = api.Dirichlet(name='vel_b0')
vel_b[0]._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_b[0]._data[CUBA.VELOCITY] = (1.0, 0, 0)
pres_b[0] = api.Neumann(name='pres_b0')
pres_b[0]._data[CUBA.VARIABLE] = CUBA.PRESSURE

vel_b[1] = api.Neumann(name='vel_b1')
vel_b[1]._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_b[1] = api.Dirichlet(name='pres_b1')
pres_b[1]._data[CUBA.VARIABLE] = CUBA.PRESSURE
pres_b[1]._data[CUBA.PRESSURE] = 0.0

vel_b[2] = api.SlipVelocity(name='vel_b2')
vel_b[2]._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_b[2] = api.Neumann(name='pres_b2')
pres_b[2]._data[CUBA.VARIABLE] = CUBA.PRESSURE

vel_b[3] = api.SlipVelocity(name='vel_b3')
vel_b[3]._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_b[3] = api.Neumann(name='pres_b3')
pres_b[3]._data[CUBA.VARIABLE] = CUBA.PRESSURE

vel_b[4] = api.Neumann(name='vel_b4')
vel_b[4]._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_b[4] = api.Dirichlet(name='pres_b4')
pres_b[4]._data[CUBA.VARIABLE] = CUBA.PRESSURE
pres_b[4]._data[CUBA.PRESSURE] = 0.0

vel_b[5] = api.Empty(name='vel_b5')
vel_b[5]._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_b[5] = api.Empty(name='pres_b5')
pres_b[5]._data[CUBA.VARIABLE] = CUBA.PRESSURE

boundaries = []
for i in range(6):
    boundaries.append(api.Boundary(name='boundary' + str(i),
                                   condition=[vel_b[i], pres_b[i]]))
cuds.add(boundaries)

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.Internal)

sim.run()

mesh_in_engine = cuds.get_by_name(mesh_name)


# compute velocity magnitude for visualization
updated_cells = []

for cell in mesh_in_engine.iter(item_type=CUBA.CELL):
    velo = cell.data[CUBA.VELOCITY]
    cell.data[CUBA.MAGNITUDE] = math.sqrt(sum(velo[i]*velo[i]
                                              for i in range(3)))
    updated_cells.append(cell)
mesh_in_engine._update_cells(updated_cells)


@mayavi2.standalone
def view(dataset):
    from mayavi.modules.surface import Surface
    from simphony_mayavi.sources.api import CUDSSource

    mayavi.new_scene()  # noqa
    src = CUDSSource(cuds=dataset)
    mayavi.add_source(src)  # noqa
    s = Surface()
    mayavi.add_module(s)  # noqa


if __name__ == '__main__':
    view(mesh_in_engine)
