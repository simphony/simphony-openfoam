"""Example to solve 2D 2 -phase poiseuille flow

"""

import foam_controlwrapper
from simphony.core.cuba import CUBA

from mayavi.scripts import mayavi2

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

import tempfile

case_name = 'poiseuille_vof'
mesh_name = 'poiseuille_mesh'

# create cuds
cuds = CUDS(name=case_name)


# physics model
cfd = api.Cfd(name='default model')

# these are already bt default se in CFD
cfd.rheology_model = api.NewtonianFluidModel(name='newtonian')
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')
# add to cuds
cuds.add(cfd)


# materials
water = api.Material(name='water')
water._data[CUBA.DENSITY] = 1000.0
water._data[CUBA.DYNAMIC_VISCOSITY] = 0.001
cuds.add(water)

air = api.Material(name='air')
air._data[CUBA.DENSITY] = 1.0
air._data[CUBA.DYNAMIC_VISCOSITY] = 1.0e-5
cuds.add(air)


# surface tension

st = api.SurfaceTensionRelation(material=[water, air],
                                surface_tension=72.86e-3)
cuds.add(st)

# free surface model

fsm = api.FreeSurfaceModel(name='vof')
cuds.add(fsm)

# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=0.1,
                               size=0.001)
cuds.add(sim_time)

# solver parameters
sp = api.SolverParameter(name='solver_parameters')
sp._data[CUBA.MAXIMUM_COURANT_NUMBER] = 0.5

cuds.add(sp)

# computational mesh
len_x = 20.0e-3
len_y = 1.0e-3
len_z = 0.1e-3
corner_points = [(0.0, 0.0, 0.0), (len_x, 0.0, 0.0),
                 (len_x, len_y, 0.0), (0.0, len_y, 0.0),
                 (0.0, 0.0, len_z), (len_x, 0.0, len_z),
                 (len_x, len_y, len_z), (0.0, len_y, len_z)]
# elements in x -direction
nex = 100
# elements in y -direction
ney = 10
mesh = foam_controlwrapper.create_quad_mesh(tempfile.mkdtemp(), mesh_name,
                                            corner_points, nex, ney, 1)
cuds.add(mesh)

# boundary conditions
vel_inlet = api.Neumann(name='vel_inlet')
vel_inlet._data[CUBA.VARIABLE] = CUBA.VELOCITY

pres_inlet = api.Dirichlet(name='pres_inlet')
pres_inlet._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_inlet._data[CUBA.DYNAMIC_PRESSURE] = 1.0

vf_inlet = api.Dirichlet(name='vf_inlet')
vf_inlet._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_inlet._data[CUBA.VOLUME_FRACTION] = 1.0

vel_outlet = api.Neumann(name='vel_outlet')
vel_outlet._data[CUBA.VARIABLE] = CUBA.VELOCITY

pres_outlet = api.Dirichlet(name='pres_outlet')
pres_outlet._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_outlet._data[CUBA.DYNAMIC_PRESSURE] = 0.0

vf_outlet = api.Neumann(name='vf_outlet')
vf_outlet._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_walls = api.Dirichlet(name='vel_walls')
vel_walls._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_walls._data[CUBA.VELOCITY] = (0, 0, 0)

pres_walls = api.Neumann(name='pres_walls')
pres_walls._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE

vf_walls = api.Neumann(name='vf_walls')
vf_walls._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_frontAndBack = api.Empty(name='vel_frontAndBack')
vel_frontAndBack._data[CUBA.VARIABLE] = CUBA.VELOCITY

pres_frontAndBack = api.Empty(name='pres_frontAndBack')
pres_frontAndBack._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE

vf_frontAndBack = api.Empty(name='vf_frontAndBack')
vf_frontAndBack._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

# boundaries
inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet, vf_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls, vf_walls])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet,
                                                vf_outlet])
frontAndBack = api.Boundary(name='frontAndBack', condition=[vel_frontAndBack,
                                                            pres_frontAndBack,
                                                            vf_frontAndBack])

cuds.add(inlet)
cuds.add(walls)
cuds.add(outlet)
cuds.add(frontAndBack)

# initial state. In VOF only one velocity and pressure field
mesh_in_cuds = cuds.get(mesh_name)
updated_cells = []
for cell in mesh_in_cuds.iter(item_type=CUBA.CELL):
    xmid = sum(mesh_in_cuds._get_point(puid).coordinates[0]
               for puid in cell.points)
    xmid /= sum(1.0 for _ in cell.points)
    if xmid < len_x/3.:
        cell.data[CUBA.VOLUME_FRACTION] = 1.0
    else:
        cell.data[CUBA.VOLUME_FRACTION] = 0.0

    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]

    updated_cells.append(cell)

mesh_in_cuds.update(updated_cells)


sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)

sim.run()

mesh_inside_wrapper = cuds.get(mesh_name)

print "Case directory ", mesh_inside_wrapper.path


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

