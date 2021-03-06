""" foam_variables

OpenFOAMvariable mappings

"""

from simphony.core.cuba import CUBA

dataNameMap = {CUBA.PRESSURE: "p",
               CUBA.DYNAMIC_PRESSURE: "p_rgh",
               CUBA.VELOCITY: "U",
               CUBA.VOLUME_FRACTION: "alpha",
               CUBA.FLUX: "phi",
               CUBA.HOMOGENIZED_STRESS_TENSOR: "Sigma",
               CUBA.STRAIN_TENSOR: "S",
               CUBA.RELATIVE_VELOCITY: "Vr",
               CUBA.DIFFUSION_VELOCITY: "Udm",
               CUBA.VOLUME_FRACTION_GRADIENT: "gradAlpha",
               CUBA.STRESS_TENSOR: "Stress",
               CUBA.MAGNITUDE: "Magnitude"
               }

dataKeyMap = {"p": CUBA.PRESSURE,
              "p_rgh": CUBA.DYNAMIC_PRESSURE,
              "U": CUBA.VELOCITY,
              "alpha": CUBA.VOLUME_FRACTION,
              "phi": CUBA.FLUX,
              "Sigma": CUBA.HOMOGENIZED_STRESS_TENSOR,
              "S": CUBA.STRAIN_TENSOR,
              "Vr": CUBA.RELATIVE_VELOCITY,
              "Udm": CUBA.DIFFUSION_VELOCITY,
              "gradAlpha": CUBA.VOLUME_FRACTION_GRADIENT,
              "Stress": CUBA.STRESS_TENSOR,
              "Magnitude": CUBA.MAGNITUDE
              }

dataTypeMap = {CUBA.PRESSURE: "scalar",
               CUBA.DYNAMIC_PRESSURE: "scalar",
               CUBA.VELOCITY: "vector",
               CUBA.VOLUME_FRACTION: "scalar",
               CUBA.FLUX: "surfaceScalar",
               CUBA.HOMOGENIZED_STRESS_TENSOR: "tensor",
               CUBA.STRAIN_TENSOR: "tensor",
               CUBA.RELATIVE_VELOCITY: "vector",
               CUBA.DIFFUSION_VELOCITY: "vector",
               CUBA.VOLUME_FRACTION_GRADIENT: "vector",
               CUBA.STRESS_TENSOR: "tensor",
               CUBA.MAGNITUDE: "scalar"
               }

dataDimensionMap = {CUBA.PRESSURE: [1, -1, -2, 0, 0, 0, 0],
                    CUBA.DYNAMIC_PRESSURE: [1, -1, -2, 0, 0, 0, 0],
                    CUBA.VELOCITY: [0, 1, -1, 0, 0, 0, 0],
                    CUBA.VOLUME_FRACTION: [0, 0, 0, 0, 0, 0, 0],
                    CUBA.FLUX: [1, 0, -1, 0, 0, 0, 0],
                    CUBA.HOMOGENIZED_STRESS_TENSOR: [1, -1, -2, 0, 0, 0, 0],
                    CUBA.STRAIN_TENSOR: [1, -1, -2, 0, 0, 0, 0],
                    CUBA.RELATIVE_VELOCITY: [0, 1, -1, 0, 0, 0, 0],
                    CUBA.DIFFUSION_VELOCITY: [0, 1, -1, 0, 0, 0, 0],
                    CUBA.VOLUME_FRACTION_GRADIENT: [0, -1, 0, 0, 0, 0, 0],
                    CUBA.STRESS_TENSOR: [1, -1, -2, 0, 0, 0, 0],
                    CUBA.MAGNITUDE: [0, -1, -1, 0, 0, 0, 0]
                    }

foamTypeMap = {"scalar": "volScalarField",
               "surfaceScalar": "surfaceScalarField",
               "vector": "volVectorField",
               "tensor": "volTensorField"}
cellDataTypes = ("scalar", "vector", "tensor")

phaseNames = ("phase1", "phase2")
