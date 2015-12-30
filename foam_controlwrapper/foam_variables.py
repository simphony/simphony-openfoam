""" foam_variables

OpenFOAMvariable mappings

"""

from simphony.core.cuba import CUBA
# from .cuba_extension import CUBAExt

dataNameMap = {CUBA.PRESSURE: "p",
               # CUBAExt.DYNAMIC_PRESSURE: "p_rgh",
               CUBA.CONCENTRATION: "p_rgh",
               CUBA.VELOCITY: "U",
               CUBA.VOLUME_FRACTION: "alpha.phase1",
               # CUBAExt.FLUX: "phi"
               CUBA.DISTRIBUTION: "phi",
               # CUBAExt.STRESS_TENSOR: "Sigma"
               CUBA.ACCELERATION: "Sigma",
               # CUBAExt.MICROSCOPIC_STRESS_TENSOR: "S"
               CUBA.ANGULAR_ACCELERATION: "S",
                # CUBAExt.RELATIVE_VELOCITY: "Vdj"
               CUBA.ANGULAR_VELOCITY: "Vdj",
                # CUBAExt.DIFFUSION_VELOCITY: "Udm"
               CUBA.FORCE: "Udm"
              }

dataKeyMap = {"p": CUBA.PRESSURE,
              "p_rgh": CUBA.CONCENTRATION,
              "U": CUBA.VELOCITY,
              "alpha.phase1": CUBA.VOLUME_FRACTION,
              "phi": CUBA.DISTRIBUTION,
              "Sigma": CUBA.ACCELERATION,
              "S": CUBA.ANGULAR_ACCELERATION,
              "Vdj": CUBA.ANGULAR_VELOCITY,
              "Udm": CUBA.FORCE
              }

dataTypeMap = {CUBA.PRESSURE: "scalar",
               CUBA.CONCENTRATION: "scalar",
               CUBA.VELOCITY: "vector",
               CUBA.VOLUME_FRACTION: "scalar",
               CUBA.DISTRIBUTION: "surfacesScalar",
               CUBA.ACCELERATION: "tensor",
               CUBA.ANGULAR_ACCELERATION: "tensor",
               CUBA.ANGULAR_VELOCITY: "vector",
               CUBA.FORCE: "vector"
               }

dataDimensionMap = {CUBA.PRESSURE: "[0, 2, -2, 0, 0, 0, 0]",
                    CUBA.CONCENTRATION: "[0, 2, -2, 0, 0, 0, 0]",
                    CUBA.VELOCITY: "[0, 1, -1, 0, 0, 0, 0]",
                    CUBA.VOLUME_FRACTION: "[0, 0, 0, 0, 0, 0, 0]",
                    CUBA.DISTRIBUTION: "[1 0 -1 0 0 0 0]",
                    CUBA.ACCELERATION: "[1 -1 -2 0 0 0 0]",
                    CUBA.ANGULAR_ACCELERATION: "[1 -1 -2 0 0 0 0]",
                    CUBA.ANGULAR_VELOCITY: "[0, 1, -1, 0, 0, 0, 0]",
                    CUBA.FORCE: "[0, 1, -1, 0, 0, 0, 0]"
                    }

foamTypeMap = {"scalar": "volScalarField",
               "surfaceScalar": "surfaceScalarField",
               "vector": "volVectorField",
               "tensor": "volTensorField"}
cellDataTypes = ("scalar", "vector", "tensor")
