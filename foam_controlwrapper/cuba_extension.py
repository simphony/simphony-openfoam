from enum import IntEnum, unique


@unique
class CUBAExt(IntEnum):

    INCOMPRESSIBLE = 1
    COMPRESSIBLE = 2
    VOF = 3
    LAMINAR_MODEL = 4
    GE = 5
