import os
import shutil
import sys

print("Post installation")

lib_dest = os.path.expandvars("$HOME/OpenFOAM/$USER-2.3.0/platforms/linux64GccDPOpt/lib")
bin_dest = os.path.expandvars("$HOME/OpenFOAM/$USER-2.3.0/platforms/linux64GccDPOpt/bin")

try:
    os.makedirs(lib_dest)
except OSError:
    # Already existent
    pass

try:
    os.makedirs(bin_dest)
except OSError:
    # Already existent
    pass

print("Installing extension files in home directory")

for f in ["driftFluxSimphonyFoam", "pimpleSimphonyFoam", "simpleSimphonyFoam"]:
    shutil.copy(os.path.join(sys.exec_prefix, "bin", f), bin_dest)

for f in [
    "libdriftFluxRelativeVelocityModelsII.so",
    "libdriftFluxRelativeVelocityModels.so",
    "libdriftFluxTransportModels.so",
    "libfoaminterface.so",
    "libincompressibleRASModelsSimphony.so",
    "libincompressibleTransportModelSimphony.so",
    "libincompressibleTurbulenceModelSimphony.so",
    "libshearStressPowerLawSlipVelocityIO.so",
    "libshearStressPowerLawSlipVelocity.so"]:

    shutil.copy(os.path.join(sys.exec_prefix, "lib", f), lib_dest)

print("Done")
