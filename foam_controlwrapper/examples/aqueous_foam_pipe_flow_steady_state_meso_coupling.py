"""Example to solve 3D aqueous foam pipe.
This example tests meso model coupling.
Meso model is analytical Herschel-Bulkley power law.
Wall shear stress dependent slip velocity law for wall layer.
"""

import numpy as np

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
foam.data[CUBA.DENSITY] = 250.0
foam.data[CUBA.DYNAMIC_VISCOSITY] = 4.37  # initial_viscosity of HB model
cuds.add([foam])

# use Herschel Bulkley viscosity model for aqueous foam
hb = api.HerschelBulkleyModel(name='foam_rheology')
hb.initial_viscosity = 0.01748 * foam.data[CUBA.DENSITY]
hb.relaxation_time = 0.0148 * foam.data[CUBA.DENSITY]
hb.linear_constant = 0.00268 * foam.data[CUBA.DENSITY]
hb.power_law_index = 0.5
hb.material = cuds.get_by_name('foam').uid
cfd.rheology_model = hb

cuds.add([cfd])

sol_par = api.SolverParameter(name='steady_state')
sol_par.data[CUBA.STEADY_STATE] = True
cuds.add([sol_par])

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
vel_inlet = api.Dirichlet(foam, name='vel_inlet')
vel_inlet.data[CUBA.VARIABLE] = CUBA.VELOCITY
vel_inlet.data[CUBA.VELOCITY] = (0, 0, 0.53)
pres_inlet = api.Neumann(foam, name='pres_inlet')
pres_inlet.data[CUBA.VARIABLE] = CUBA.PRESSURE

vel_outlet = api.Neumann(foam, name='vel_outlet')
vel_outlet.data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_outlet = api.Dirichlet(foam, name='pres_outlet')
pres_outlet.data[CUBA.VARIABLE] = CUBA.PRESSURE
pres_outlet.data[CUBA.PRESSURE] = 0.0

vel_walls = api.ShearStressPowerLawSlipVelocity(foam, 
                                                density = 250.0,
                                                linear_constant = 3.1e-3,
                                                power_law_index = 1.16,
                                                name='vel_walls')
vel_walls.data[CUBA.VARIABLE] = CUBA.VELOCITY
pres_walls = api.Neumann(name='pres_walls')
pres_walls.data[CUBA.VARIABLE] = CUBA.PRESSURE

inlet = api.Boundary(name='inlet', condition=[vel_inlet, pres_inlet])
walls = api.Boundary(name='walls', condition=[vel_walls, pres_walls])
outlet = api.Boundary(name='outlet', condition=[vel_outlet, pres_outlet])

cuds.add([inlet, walls, outlet])

end = time.time()
print "Time spend in boundary settings: ", end-start

start = time.time()
sim = Simulation(cuds, 'OpenFOAM', engine_interface=EngineInterface.FileIO)
end = time.time()
print "Time spend in Simulation initialization: ", end-start


# time setting
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=20,
                               size=1)
cuds.add([sim_time])

sim.run()

mesh_in_engine = cuds.get_by_name(mesh_name)


sm = api.MesoscopicStressModel(name='meso_stress_model')
cuds.add([sm])

cuds.remove('simulation_time')
sim_time = api.IntegrationTime(name='simulation_time',
                               current=0.0,
                               final=1,
                               size=1)
cuds.add([sim_time])


start = time.time()

number_of_outer_timesteps = 500
print "Working directory: ", mesh_in_engine.path
for time_i in range(number_of_outer_timesteps):
    # solve macroscopic scale
    print "Solve cfd"
    sim.run()
    print "Time: ", mesh_in_engine._time
    print "Mesoscale as analytic coupling"
    print " Update stress"

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
        hb_stress = cell.data[CUBA.STRESS_TENSOR]
        cell.data[CUBA.HOMOGENIZED_STRESS_TENSOR] = hb_stress
        updated_cells.append(cell)
    mesh_in_engine._update_cells(updated_cells)


end = time.time()
print "Time spend in run: ", end-start

start = time.time()
print "Working directory ", mesh_in_engine.path
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
