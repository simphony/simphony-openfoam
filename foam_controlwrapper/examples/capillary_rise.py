"""Example to solve 2D 2 -phase poiseuille flow

"""


import foam_controlwrapper
from simphony.core.cuba import CUBA

# from mayavi.scripts import mayavi2

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface


import tube_mesh

import tempfile


case_name = 'capillary_rise'
mesh_name = 'capillary_rise_mesh'

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
                               final=0.5,
                               size=1.0e-5)
cuds.add(sim_time)

# solver parameters
sp = api.SolverParameter(name='solver_parameters')
sp._data[CUBA.MAXIMUM_COURANT_NUMBER] = 0.2
cuds.add(sp)

gm = api.GravityModel(name='gravitation')
gm.acceleration = (0, -9.81, 0)
cuds.add(gm)

# boundary conditions
vel_inlet = api.InletOutlet(name='vel_inlet')
vel_inlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_inlet._data[CUBA.VELOCITY] = (0, 0, 0)
pres_inlet = api.Dirichlet(name='pres_inlet')
pres_inlet._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_inlet._data[CUBA.DYNAMIC_PRESSURE] = 0.0
vf_inlet = api.InletOutlet(name='vf_inlet')
vf_inlet._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_inlet._data[CUBA.VOLUME_FRACTION] = 1.0

vel_atm = api.InletOutlet(name='vel_atm')
vel_atm._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_atm._data[CUBA.VELOCITY] = (0, 0, 0)
pres_atm = api.Dirichlet(name='pres_atm')
pres_atm._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_atm._data[CUBA.DYNAMIC_PRESSURE] = 0.0
vf_atm = api.Neumann(name='vf_atm')
vf_atm._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_walls = api.Dirichlet(name='vel_walls')
vel_walls._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_walls._data[CUBA.VELOCITY] = (0, 0, 0)
pres_walls = api.Neumann(name='pres_walls')
pres_walls._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_walls = api.WettingAngle(name='vf_walls')
vf_walls.contact_angle = 45.0
vf_walls._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_frontAndBack = api.Empty(name='vel_frontAndBack')
vel_frontAndBack._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_frontAndBack = api.Empty(name='pres_frontAndBack')
pres_frontAndBack._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_frontAndBack = api.Empty(name='vf_frontAndBack')
vf_frontAndBack._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet, vf_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls, vf_walls])
atm = api.Boundary(name='atmosphere', condition=[vel_atm, pres_atm,
                                                 vf_atm])
frontAndBack = api.Boundary(name='frontAndBack', condition=[vel_frontAndBack,
                                                            pres_frontAndBack,
                                                            vf_frontAndBack])

cuds.add(inlet)
cuds.add(walls)
cuds.add(atm)
cuds.add(frontAndBack)


# create mesh
mesh = foam_controlwrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                             tube_mesh.blockMeshDict)

cuds.add(mesh)
# initial state. In VOF only one velocity and pressure field

mesh_in_cuds = cuds.get(mesh_name)

updated_cells = []
for cell in mesh_in_cuds._iter_cells():
    ymid = sum(mesh_in_cuds.get_point(puid).coordinates[1]
               for puid in cell.points)
    ymid /= sum(1.0 for _ in cell.points)
    if ymid < 8e-3:
        cell.data[CUBA.VOLUME_FRACTION] = 1.0
    else:
        cell.data[CUBA.VOLUME_FRACTION] = 0.0

    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]

    updated_cells.append(cell)

mesh_in_cuds._update_cells(updated_cells)

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)

sim.run()
mesh_inside_wrapper = cuds.get(mesh_name)

print "Case directory ", mesh_inside_wrapper.path

"""
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
"""
