"""Example to solve 3D aqueous foam pipe flow using rheological
Herschel-Bulkley power law for bulk and wall shear stress dependent
slip velocity law for wall layer

"""
import numpy as np

import foam_controlwrapper
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api
from simphony.engine import EngineInterface

from mayavi.scripts import mayavi2

import slit_mesh
import tempfile
import time

start = time.time()

case_name = 'aqueous_foam'
mesh_name = 'aqueous_foam_mesh'

cuds = CUDS(name=case_name)

# physics model
cfd = api.Cfd(name='cfd model')

# these are already bt default set in CFD
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')

# material
foam = api.Material(name='foam')
foam._data[CUBA.DENSITY] = 250.0
foam._data[CUBA.DYNAMIC_VISCOSITY] = 4.37  # initial_viscosity of HB model
cuds.add([foam])

# use Herschel Bulkley viscosity model for aqueous foam
hb = api.HerschelBulkleyModel(name='foam_rheology')
hb.initial_viscosity = 0.01748 * foam._data[CUBA.DENSITY]
hb.relaxation_time = 0.0148 * foam._data[CUBA.DENSITY]
hb.linear_constant = 0.00268 * foam._data[CUBA.DENSITY]
hb.power_law_index = 0.5
hb.material = cuds.get_by_name('foam').uid

cfd.rheology_model = hb

cuds.add([cfd])


sol_par = api.SolverParameter(name='steady_state')
sol_par._data[CUBA.STEADY_STATE] = True
cuds.add([sol_par])

end = time.time()
print "Time spend in initialization: ", end-start

start = time.time()
# create computational mesh
mesh = foam_controlwrapper.create_block_mesh(tempfile.mkdtemp(), mesh_name,
                                             slit_mesh.blockMeshDict)
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
vel_inlet._data[CUBA.VELOCITY] = (0.53, 0, 0)
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

vel_frontAndBack = api.Empty(name='vel_frontAndBack')
vel_frontAndBack._data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_frontAndBack = api.Empty(name='pres_frontAndBack')
pres_frontAndBack._data[CUBA.VARIABLE] = CUBA.PRESSURE


inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet])
frontAndBack = api.Boundary(name='frontAndBack',
                            condition=[vel_frontAndBack, pres_frontAndBack])

cuds.add([inlet, walls, outlet, frontAndBack])

end = time.time()
print "Time spend in boundary settings: ", end-start

start = time.time()
sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.Internal)
end = time.time()
print "Time spend in Simulation initialization: ", end-start


# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=67,
                               size=1)
cuds.add([sim_time])
sim.run()


mesh_in_engine = cuds.get_by_name(mesh_name)


sm = api.MesoscopicStressModel(name='meso_stress_model')
cuds.add([sm])

cuds.remove([sim_time.uid])
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=1,
                               size=1)
cuds.add([sim_time])

start = time.time()

number_of_outer_timesteps = 20

for time_i in range(number_of_outer_timesteps):
    updated_cells = []
    for cell in mesh_in_engine._iter_cells():
        strain = cell.data[CUBA.STRAIN_TENSOR]
        strain_rate = np.linalg.norm(strain)
        nu = min(hb.initial_viscosity,
                 (hb.relaxation_time +
                  hb.linear_constant*pow(strain_rate, hb.power_law_index))
                 / (max(strain_rate, 1.0e-6))
                 )
        hb_stress = [nu*sri for sri in strain]
        cell.data[CUBA.HOMOGENIZED_STRESS_TENSOR] = hb_stress
        updated_cells.append(cell)
    mesh_in_engine._update_cells(updated_cells)

    # solve macroscopic scale
    print "Solve cfd"
    sim.run()
    print "Time: ", mesh_in_engine._time
    print "Mesoscale as analytic coupling"
    print " Update stress"

end = time.time()
print "Time spend in run: ", end-start

start = time.time()

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
