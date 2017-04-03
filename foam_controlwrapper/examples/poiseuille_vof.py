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

# these are already bt default set in CFD
cfd.rheology_model = api.NewtonianFluidModel(name='newtonian')
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')
# add to cuds
cuds.add([cfd])


# materials
water = api.Material(name='water')
water.data[CUBA.DENSITY] = 1000.0
water.data[CUBA.DYNAMIC_VISCOSITY] = 0.001
cuds.add([water])

air = api.Material(name='air')
air.data[CUBA.DENSITY] = 1.0
air.data[CUBA.DYNAMIC_VISCOSITY] = 1.0e-5
cuds.add([air])


# surface tension

st = api.SurfaceTensionRelation(material=[water, air],
                                surface_tension=72.86e-3)
cuds.add([st])

# free surface model

fsm = api.FreeSurfaceModel(name='vof')
cuds.add([fsm])

# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=0.1,
                               size=0.001)
cuds.add([sim_time])

# solver parameters
sp = api.SolverParameter(name='solver_parameters')
sp.data[CUBA.MAXIMUM_COURANT_NUMBER] = 0.5

cuds.add([sp])

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
cuds.add([mesh])

# boundary conditions
vel_inlet = api.Neumann(water, name='vel_inlet')
vel_inlet.data[CUBA.VARIABLE] = CUBA.VELOCITY

#pres_inlet = api.ConstantPressureCondition(10.0, water, name='pres_inlet')
pres_inlet = api.Dirichlet(water, name='pres_inlet')
pres_inlet.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_inlet.data[CUBA.DYNAMIC_PRESSURE] = 10.0

vf_inlet = api.Dirichlet(water, name='vf_inlet')
vf_inlet.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_inlet.data[CUBA.VOLUME_FRACTION] = 1.0

vel_outlet = api.Neumann(water, name='vel_outlet')
vel_outlet.data[CUBA.VARIABLE] = CUBA.VELOCITY

pres_outlet = api.Dirichlet(water, name='pres_outlet')
pres_outlet.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_outlet.data[CUBA.DYNAMIC_PRESSURE] = 0.0

vf_outlet = api.Neumann(water, name='vf_outlet')
vf_outlet.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_walls = api.Dirichlet(water, name='vel_walls')
vel_walls.data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_walls.data[CUBA.VELOCITY] = (0, 0, 0)

pres_walls = api.Neumann(water, name='pres_walls')
pres_walls.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE

vf_walls = api.Neumann(water, name='vf_walls')
vf_walls.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_frontAndBack = api.EmptyCondition(name='vel_frontAndBack')
vel_frontAndBack.data[CUBA.VARIABLE] = CUBA.VELOCITY

pres_frontAndBack = api.EmptyCondition(name='pres_frontAndBack')
pres_frontAndBack.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE

vf_frontAndBack = api.EmptyCondition(name='vf_frontAndBack')
vf_frontAndBack.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

# boundaries
inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet, vf_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls, vf_walls])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet,
                                                vf_outlet])
frontAndBack = api.Boundary(name='frontAndBack', condition=[vel_frontAndBack,
                                                            pres_frontAndBack,
                                                            vf_frontAndBack])

cuds.add([inlet, walls, outlet, frontAndBack])

# initial state. In VOF only one velocity and pressure field
mesh_in_cuds = cuds.get_by_name(mesh_name)
updated_cells = []
zero_water = api.PhaseVolumeFraction(water, 0)
zero_air= api.PhaseVolumeFraction(air, 0)
one_water = api.PhaseVolumeFraction(water, 1)
one_air = api.PhaseVolumeFraction(air, 1)
for cell in mesh_in_cuds.iter(item_type=CUBA.CELL):
#for cell in mesh_in_cuds._iter_cells_parallel():
    xmid = sum(mesh_in_cuds._get_point(puid).coordinates[0]
               for puid in cell.points)
    xmid /= sum(1.0 for _ in cell.points)
    if xmid < len_x/3.:
        cell.data[CUBA.VOLUME_FRACTION] = [one_water, zero_air]
    else:
        cell.data[CUBA.VOLUME_FRACTION] = [zero_water, one_air]

    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]

    updated_cells.append(cell)

mesh_in_cuds._update_cells(updated_cells)

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)

sim.run()

mesh_in_engine = cuds.get_by_name(mesh_name)

print "Case directory ", mesh_in_engine.path


@mayavi2.standalone
def view():
    from mayavi.modules.surface import Surface
    from simphony_mayavi.sources.api import CUDSSource

    mayavi.new_scene()  # noqa
    src = CUDSSource(cuds=mesh_in_engine)
    mayavi.add_source(src)  # noqa
    s = Surface()
    mayavi.add_module(s)  # noqa


if __name__ == '__main__':
    view()

