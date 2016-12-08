"""Example to solve 2D poiseuille flow

"""

from simphony.core.cuba import CUBA

from mayavi.scripts import mayavi2

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

import tempfile
from foam_controlwrapper import create_quad_mesh

case_name = 'poiseuille'
mesh_name = 'poiseuille_mesh'

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
cuds.add(cfd)

# material
mat = api.Material(name='a_material')
mat._data[CUBA.DENSITY] = 1.0
mat._data[CUBA.DYNAMIC_VISCOSITY] = 1.0
cuds.add(mat)

# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=1.0,
                               size=0.01)
cuds.add(sim_time)

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

mesh = create_quad_mesh(tempfile.mkdtemp(), mesh_name,
                        corner_points, nex, ney, 1)
cuds.add(mesh)

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
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet])
frontAndBack = api.Boundary(name='frontAndBack',
                            condition=[vel_frontAndBack, pres_frontAndBack])

cuds.add(inlet)
cuds.add(walls)
cuds.add(outlet)
cuds.add(frontAndBack)

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)

sim.run()

average_pressure = 0.0
mesh_inside_wrapper = cuds.get(mesh_name)
print "working directory ", mesh_inside_wrapper.path
for cell in mesh_inside_wrapper.get_boundary_cells(inlet.name):
    average_pressure += cell.data[CUBA.PRESSURE]

average_pressure /= len(mesh_inside_wrapper._boundaries[inlet.name])

print "Average pressure on " + inlet.name + ": ", average_pressure


@mayavi2.standalone
def view():
    from mayavi.modules.surface import Surface
    from simphony_mayavi.sources.api import CUDSSource

    mayavi.new_scene()  # noqa
    src = CUDSSource(cuds=mesh_inside_wrapper)
    mayavi.add_source(src)  # noqa
    s = Surface()
    mayavi.add_module(s)  # noqa


if __name__ == '__main__':
    view()

