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
                               final=640.0,
                               size=0.1)
cuds.add([sim_time])


# create computational mesh
mesh = foam_controlwrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                             dahl_mesh.blockMeshDict)
cuds.add([mesh])


# materials
sludge = api.Material(name='sludge')
sludge.data[CUBA.DENSITY] = 1900.0
sludge.data[CUBA.DYNAMIC_VISCOSITY] = 0.01
cuds.add([sludge])

water = api.Material(name='water')
water.data[CUBA.DENSITY] = 1000.0
water.data[CUBA.DYNAMIC_VISCOSITY] = 1.0e-3
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

rvm = api.SimpleRelativeVelocityModel(name='simple_rel_vel_model')
rvm.diffusion_velocity = (0, -0.002, 0)
rvm.linear_constant = 285.0
cuds.add([rvm])

gm = api.GravityModel(name='gravitation')
gm.acceleration = (0, -9.81, 0)
cuds.add([gm])

# boundary condtitions
vel_inlet = api.Dirichlet(water, name='vel_inlet')
vel_inlet.data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_inlet.data[CUBA.VELOCITY] = (0.0191, 0, 0)
pres_inlet = api.Neumann(water, name='pres_inlet')
pres_inlet.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_inlet = api.Dirichlet(water, name='vf_inlet')
vf_inlet.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_inlet.data[CUBA.VOLUME_FRACTION] = 0.001

vel_outlet = api.InletOutlet(name='vel_outlet')
vel_outlet.data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_outlet.data[CUBA.VELOCITY] = (0, 0, 0)

pres_outlet = api.Dirichlet(water, name='pres_outlet')
pres_outlet.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
pres_outlet.data[CUBA.DYNAMIC_PRESSURE] = 0.0
vf_outlet = api.InletOutlet(name='vf_outlet')
vf_outlet.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION
vf_outlet.data[CUBA.VOLUME_FRACTION] = 0.001

vel_walls = api.Dirichlet(water, name='vel_walls')
vel_walls.data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_walls.data[CUBA.VELOCITY] = (0, 0, 0)
pres_walls = api.Neumann(water, name='pres_walls')
pres_walls.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_walls = api.Neumann(water, name='vf_walls')
vf_walls.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_top = api.SlipVelocity([0, 0, 0], water, name='vel_top')
vel_top.data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_top = api.Neumann(water, name='pres_top')
pres_top.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_top = api.Neumann(water, name='vf_top')
vf_top.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

vel_frontAndBack = api.EmptyCondition(name='vel_frontAndBack')
vel_frontAndBack.data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_frontAndBack = api.EmptyCondition(name='pres_frontAndBack')
pres_frontAndBack.data[CUBA.VARIABLE] = CUBA.DYNAMIC_PRESSURE
vf_frontAndBack = api.EmptyCondition(name='vf_frontAndBack')
vf_frontAndBack.data[CUBA.VARIABLE] = CUBA.VOLUME_FRACTION

# boundaries
inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet, vf_inlet])
bottom_wall = api.Boundary(name='bottomWall', condition=[vel_walls,
                                                         pres_walls,
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
vf_sludge = api.PhaseVolumeFraction(sludge, 0.001)
vf_water = api.PhaseVolumeFraction(water, 1 - 0.001)

for cell in mesh_in_cuds.iter(item_type=CUBA.CELL):
    cell.data[CUBA.VOLUME_FRACTION] = [vf_sludge, vf_water]
    cell.data[CUBA.DYNAMIC_PRESSURE] = 0.0
    cell.data[CUBA.VELOCITY] = [0.0191, 0.0, 0.0]
    updated_cells.append(cell)

mesh_in_cuds._update_cells(updated_cells)

sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.Internal)

sim.run()

mesh_in_engine = cuds.get_by_name(mesh_name)

print "Working directory: ", mesh_in_engine.path

average_vf = 0.0
for cell in mesh_in_engine.get_boundary_cells(outlet.name):
    pvf = cell.data[CUBA.VOLUME_FRACTION]
    i = 0
    if pvf[1].material == sludge:
        i = 1
    average_vf += pvf[i].volume_fraction

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
