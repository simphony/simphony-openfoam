"""Example to solve mixture model
"""

import foam_controlwrapper
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

from mayavi.scripts import mayavi2

import dahl_mesh
import tempfile

case_name = 'dahl'
mesh_name = 'dahl_mesh'

cuds = CUDS(name=case_name)

# physics model
cfd = api.Cfd(name='default model')

# these are already by default set in CFD
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')

# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=0.1,
                               size=0.1)
cuds.add([sim_time])

# create computational mesh
mesh = foam_controlwrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                             dahl_mesh.blockMeshDict)
cuds.add([mesh])

# materials
sludge = api.Material(name='sludge')
sludge._data[CUBA.DENSITY] = 1900.0
sludge._data[CUBA.DYNAMIC_VISCOSITY] = 0.01
cuds.add([sludge])

water = api.Material(name='water')
water._data[CUBA.DENSITY] = 1000.0
water._data[CUBA.DYNAMIC_VISCOSITY] = 1.0e-3
cuds.add([water])

# mixture model
mm = api.MixtureModel(name='mixture')
mm.disperse = cuds.get_by_name('sludge').uid
cuds.add([mm])

# use Bingham plastic for sludge rheology model
bp = api.BinghamPlasticModel(name='sludge_rheology')
bp.linear_constant = (0.00023143, 0.0005966)
bp.power_law_index = (179.26, 1050.8)
bp.maximum_viscosity = 10.0
bp.material = cuds.get_by_name('sludge').uid

# water is Newtonian so no need to add that to rheology models
cfd.rheology_model = bp

cuds.add([cfd])

sm = api.StandardStressModel(name='standard_stress_model')
cuds.add([sm])

rvm = api.MesoScaleRelativeVelocityModel(name='meso_rel_vel_model')
cuds.add([rvm])

gm = api.GravityModel(name='gravitation')
gm.acceleration = (0, -9.81, 0)
cuds.add([gm])

vel_inlet = api.Dirichlet(name='vel_inlet')
vel_inlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_inlet._data[CUBA.VELOCITY] = (0.0191, 0, 0)
pres_inlet = api.Neumann(name='pres_inlet')
pres_inlet._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_inlet = api.Dirichlet(name='vf_inlet')
vf_inlet._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_inlet._data[CUBA.VOLUME_FRACTION] = 0.001

vel_outlet = api.InletOutlet(name='vel_outlet')
vel_outlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_outlet._data[CUBA.VELOCITY] = (0, 0, 0)
pres_outlet = api.Dirichlet(name='pres_outlet')
pres_outlet._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_outlet._data[CUBA.DYNAMIC_PRESSURE] = 0.0
vf_outlet = api.InletOutlet(name='vf_outlet')
vf_outlet._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_outlet._data[CUBA.VOLUME_FRACTION] = 0.001

vel_walls = api.Dirichlet(name='vel_walls')
vel_walls._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_walls._data[CUBA.VELOCITY] = (0, 0, 0)
pres_walls = api.Neumann(name='pres_walls')
pres_walls._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_walls = api.Neumann(name='vf_walls')
vf_walls._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_top = api.SlipVelocity(name='vel_top')
vel_top._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_top = api.Neumann(name='pres_top')
pres_top._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_top = api.Neumann(name='vf_top')
vf_top._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_frontAndBack = api.Empty(name='vel_frontAndBack')
vel_frontAndBack._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_frontAndBack = api.Empty(name='pres_frontAndBack')
pres_frontAndBack._data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_frontAndBack = api.Empty(name='vf_frontAndBack')
vf_frontAndBack._data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

# boundaries
inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet, vf_inlet])
bottom_wall = api.Boundary(name='bottomWall', condition=[vel_walls, pres_walls,
                                                         vf_walls])
end_wall = api.Boundary(name='endWall', condition=[vel_walls, pres_walls,
                                                   vf_walls])
top = api.Boundary(name='top', condition=[vel_top, pres_top, vf_top])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet,
                                                vf_outlet])
frontAndBack = api.Boundary(name='frontAndBack', condition=[vel_frontAndBack,
                                                            pres_frontAndBack,
                                                            vf_frontAndBack])

cuds.add([inlet, bottom_wall, end_wall, top, outlet, frontAndBack])

mesh_in_cuds = cuds.get_by_name(mesh_name)

updated_cells = []
for cell in mesh_in_cuds.iter(item_type=CUBA.CELL):
    cell.data[CUBA.VOLUME_FRACTION] = 0.001
    cell.data[CUBA.MATERIAL] = water.uid
    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0191, 0.0, 0.0]
    updated_cells.append(cell)

mesh_in_cuds._update_cells(updated_cells)

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.Internal)

# must get mesh again, while mesh is changed on cuds on addition to engine
mesh_in_engine = cuds.get_by_name(mesh_name)

V0 = [0.0, -0.002, 0.0]
a = 285.0

number_of_outer_timesteps = 6400

for time_i in range(number_of_outer_timesteps):
    # solve macroscopic scale
    print "Solve cfd"
    sim.run()
    print "Time: ", mesh_in_engine._time
    print "Mesoscale as analytic coupling"
    print " Update relative velocity"

    updated_cells = []
    for cell in mesh_in_engine.iter(item_type=CUBA.CELL):
        alphad = cell.data[CUBA.VOLUME_FRACTION]
        vr = [V*pow(10.0, -a*max(alphad, 0.0))/(1-alphad) for V in V0]
        cell.data[CUBA.RELATIVE_VELOCITY] = vr
        updated_cells.append(cell)

    mesh_in_engine._update_cells(updated_cells)


average_vf = 0.0
for cell in mesh_in_engine.get_boundary_cells(outlet.name):
    average_vf += cell.data[CUBA.VOLUME_FRACTION]

average_vf /= len(mesh_in_engine._boundaries[outlet.name])

print "Average volume fraction on outlet: ", average_vf


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
