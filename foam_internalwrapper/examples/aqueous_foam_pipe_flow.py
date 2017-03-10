"""Example to solve 3D aqueous foam pipe flow using rheological
Herschel-Bulkley power law for bulk and wall shear stress dependent
slip velocity law for wall layer

"""

import foam_controlwrapper
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

from mayavi.scripts import mayavi2

import pipe_mesh
import tempfile
import time

start = time.time()

case_name = 'aqueous_foam'
mesh_name = 'aqueous_foam_mesh'

cuds = CUDS(name=case_name)

# physics model
cfd = api.Cfd(name='default model')

# these are already bt default set in CFD
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')

# material
foam = api.Material(name='foam')
foam._data[CUBA.DENSITY] = 250.0
foam._data[CUBA.DYNAMIC_VISCOSITY] = 4.37
cuds.add([foam])

# use Herschel Bulkley viscosity model for aqueous foam
hb = api.HerschelBulkleyModel(name='foam_rheology')
hb.initial_viscosity = 0.01748
hb.relaxation_time = 0.0148
hb.linear_constant = 0.00268
hb.power_law_index = 0.5
hb.material = cuds.get_by_name('foam').uid
cfd.rheology_model = hb

cuds.add([cfd])

# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=0.3,
                               size=0.0001)
cuds.add([sim_time])

end = time.time()
print "Time spend in initialization: ", end-start

start = time.time()
# create computational mesh
mesh = foam_controlwrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                             pipe_mesh.blockMeshDict)
end = time.time()
print "Time spend in blockmesh: ", end-start

start = time.time()
cuds.add([mesh])
end = time.time()
print "Time spend in add mesh to cuds: ", end-start

start = time.time()
# boundary conditions
vel_inlet = api.Dirichlet(name='vel_inlet')
vel_inlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_inlet._data[CUBA.VELOCITY] = (0, 0, 0.53)
pres_inlet = api.Neumann(name='pres_inlet')
pres_inlet._data[CUBA.VARIABLE] = CUBA.PRESSURE

vel_outlet = api.Neumann(name='vel_outlet')
vel_outlet._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_outlet = api.Dirichlet(name='pres_outlet')
pres_outlet._data[CUBA.VARIABLE] = CUBA.PRESSURE
pres_outlet._data[CUBA.PRESSURE] = 0.0

vel_walls = api.ShearStressPowerLawSlipVelocity(name='vel_walls')
vel_walls.density = 250.0
vel_walls.linear_constant = 3.1e-3
vel_walls.power_law_index = 1.16
vel_walls._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_walls = api.Neumann(name='pres_walls')
pres_walls._data[CUBA.VARIABLE] = CUBA.PRESSURE

inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet])

cuds.add([inlet, walls, outlet])

end = time.time()
print "Time spend in boundary settings: ", end-start

start = time.time()
sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.Internal)
end = time.time()
print "Time spend in Simulation initialization: ", end-start

start = time.time()

sim.run()
end = time.time()
print "Time spend in run: ", end-start

start = time.time()
mesh_in_engine = cuds.get_by_name(mesh_name)

average_pressure = 0.0
for cell in mesh_in_engine.get_boundary_cells(inlet.name):
    average_pressure += cell.data[CUBA.PRESSURE]

average_pressure /= len(mesh_in_engine._boundaries[inlet.name])
end = time.time()
print "Time spend in post processing: ", end-start

print "Average pressure on inlet: ", average_pressure


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
