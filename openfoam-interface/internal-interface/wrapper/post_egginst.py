import os
import shutil
import sys

print("Post installation procedure. "
      "Installing openfoam to the user directory.")

lib_dest = os.environ.get("FOAM_USER_LIBBIN")
if lib_dest is None:
    lib_dest = os.path.expandvars(
        "$HOME/OpenFOAM/$USER-2.3.0/platforms/linux64GccDPOpt/lib")
    print("Could not obtain FOAM_USER_LIBBIN from environment. "
          "Using best guess {}".format(lib_dest))

bin_dest = os.environ.get("FOAM_USER_APPBIN")
if bin_dest is None:
    bin_dest = os.path.expandvars(
        "$HOME/OpenFOAM/$USER-2.3.0/platforms/linux64GccDPOpt/bin")
    print("Could not obtain FOAM_USER_APPBIN from environment. "
          "Using best guess {}".format(bin_dest))

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

print("Installing extension files")

# We have no way of reliably know which files have been created
# and are available. We can only manually list them.
for f in ["driftFluxSimphonyFoam", "pimpleSimphonyFoam", "simpleSimphonyFoam"]:
    # Put them from the egg deployed paths to the designated positions.
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
        "libshearStressPowerLawSlipVelocity.so"
        ]:

    # Same
    shutil.copy(os.path.join(sys.exec_prefix, "lib", f), lib_dest)

print("""Done.

Please note that the following entries must be added to your $HOME/.bashrc:

  source /opt/OpenFOAM-2.3.0/etc/bashrc
  export PATH=$PATH:{bin_dest}
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{lib_dest}

without these steps, importing simphonyfoaminterface at the python prompt
will fail to import the appropriate .so files. Verify the correct behavior
by starting python and performing

    >>> import simphonyfoaminterface

The prompt should reappear with no output.

""".format(bin_dest=bin_dest, lib_dest=lib_dest))
