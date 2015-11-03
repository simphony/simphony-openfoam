""" foam_variables

OpenFOAMvariable mappings

"""

from simphony.core.cuba import CUBA
# from .cuba_extension import CUBAExt

dataNameMap = {CUBA.PRESSURE: "p",
               # CUBAExt.DYNAMIC_PRESSURE: "p_rgh",
               CUBA.CONCENTRATION: "p_rgh",
               CUBA.VELOCITY: "U",
               CUBA.VOLUME_FRACTION: "alpha1",
               # CUBAExt.FLUX: "phi"
               CUBA.ANGULAR_VELOCITY: "phi"
               }

dataKeyMap = {"p": CUBA.PRESSURE,
              # "p_rgh": CUBAExt.DYNAMIC_PRESSURE,
              "p_rgh": CUBA.CONCENTRATION,
              "U": CUBA.VELOCITY,
              "alpha1": CUBA.VOLUME_FRACTION,
              # "phi": CUBAExt.FLUX
              "phi": CUBA.ANGULAR_VELOCITY
              }

dataTypeMap = {CUBA.PRESSURE: "scalar",
               # CUBAExt.DYNAMIC_PRESSURE: "scalar",
               CUBA.CONCENTRATION: "scalar",
               CUBA.VELOCITY: "vector",
               CUBA.VOLUME_FRACTION: "scalar",
               # CUBAExt.FLUX: "scalar"
               CUBA.ANGULAR_VELOCITY: "scalar"
               }

dataDimensionMap = {CUBA.PRESSURE: "[0, 2, -2, 0, 0, 0, 0]",
                    # CUBAExt.DYNAMIC_PRESSURE: "[0, 2, -2, 0, 0, 0, 0]",
                    CUBA.CONCENTRATION: "[0, 2, -2, 0, 0, 0, 0]",
                    CUBA.VELOCITY: "[0, 1, -1, 0, 0, 0, 0]",
                    CUBA.VOLUME_FRACTION: "[0, 0, 0, 0, 0, 0, 0]",
                    # CUBAExt.FLUX: "[1 0 -1 0 0 0 0]"
                    CUBA.ANGULAR_VELOCITY: "[1 0 -1 0 0 0 0]"
                    }

foamTypeMap = {"scalar": "volScalarField",
               "vector": "volVectorField"}
