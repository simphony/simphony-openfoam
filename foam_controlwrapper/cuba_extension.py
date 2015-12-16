""" Provisional CUBA keywords specific for this revision

"""

from enum import IntEnum, unique


@unique
class CUBAExt(IntEnum):

    INCOMPRESSIBLE = 1
    COMPRESSIBLE = 2
    VOF = 3
    LAMINAR_MODEL = 4
    GE = 5
    PATCH_TYPE = 6
    PHASE_LIST = 7
    MAX_COURANT_NUMBER = 8
    SURFACE_TENSION = 9
    NUMBER_OF_CORES = 10
    DYNAMIC_PRESSURE = 11
    FLUX = 12