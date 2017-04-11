"""Example to solve 2D 2 -phase poiseuille flow

"""


import foam_controlwrapper
from simphony.core.cuba import CUBA

from mayavi.scripts import mayavi2

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

import tube_mesh

import tempfile

import time

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
                               final=0.5,
                               size=1.0e-5)
cuds.add([sim_time])

# solver parameters
sp = api.SolverParameter(name='solver_parameters')
sp.data[CUBA.MAXIMUM_COURANT_NUMBER] = 0.2
cuds.add([sp])

gm = api.GravityModel(name='gravitation')
gm.acceleration = (0, -9.81, 0)
cuds.add([gm])

# boundary conditions
vel_inlet = api.InletOutletVelocity((0, 0, 0), water, name='vel_inlet')
pres_inlet = api.ConstantPressureCondition(0.0, water, name='pres_inlet')
vf_inlet = api.InletOutletVolumeFraction(1.0, water, name='vf_inlet')

vel_atm = api.InletOutletVelocity((0, 0, 0), water, name='vel_atm')
pres_atm = api.ConstantPressureCondition(0.0, water, name='pres_atm')
vf_atm = api.ZeroGradientVolumeFractionCondition(0.0, water, name='vf_atm')

vel_walls = api.ConstantVelocityCondition((0, 0, 0), water, name='vel_walls')
pres_walls = api.ZeroGradientPressureCondition(0.0, water, name='pres_walls')

vf_walls = api.WettingAngle([water, air], contact_angle=45.0, name='vf_walls')

vel_frontAndBack = api.EmptyCondition(name='vel_frontAndBack')
vel_frontAndBack.data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_frontAndBack = api.EmptyCondition(name='pres_frontAndBack')
pres_frontAndBack.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_frontAndBack = api.EmptyCondition(name='vf_frontAndBack')
vf_frontAndBack.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet, vf_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls, vf_walls])
atm = api.Boundary(name='atmosphere', condition=[vel_atm, pres_atm,
                                                 vf_atm])
frontAndBack = api.Boundary(name='frontAndBack', condition=[vel_frontAndBack,
                                                            pres_frontAndBack,
                                                            vf_frontAndBack])

cuds.add([inlet, walls, atm, frontAndBack])

# create mesh
mesh = foam_controlwrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                             tube_mesh.blockMeshDict)
cuds.add([mesh])

mesh_in_cuds = cuds.get_by_name(mesh_name)

start = time.time()

updated_cells = []
zero_water = api.PhaseVolumeFraction(water, 0)
zero_air = api.PhaseVolumeFraction(air, 0)
one_water = api.PhaseVolumeFraction(water, 1)
one_air = api.PhaseVolumeFraction(air, 1)

# initial state. In VOF only one velocity and pressure field
# for cell in mesh_in_cuds._iter_cells_parallell():
for cell in mesh_in_cuds.iter(item_type=CUBA.CELL):
    ymid = sum(mesh_in_cuds._get_point(puid).coordinates[1]
               for puid in cell.points)
    ymid /= sum(1.0 for _ in cell.points)
    if ymid < 8e-3:
        cell.data[CUBA.VOLUME_FRACTION] = [one_water, zero_air]
    else:
        cell.data[CUBA.VOLUME_FRACTION] = [zero_water, one_air]

    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0, 0.0, 0.0]

    updated_cells.append(cell)

print "Time spend in initialization 1 : ", time.time()-start
start = time.time()

mesh_in_cuds._update_cells(updated_cells)

print "Time spend in initialization 2 : ", time.time()-start
start = time.time()

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)

print "Time spend in wrapper init: ", time.time()-start
start = time.time()
mesh_in_engine = cuds.get_by_name(mesh_name)

print "Case directory ", mesh_in_engine.path

sim.run()
print "Time spend in run: ", time.time()-start


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
