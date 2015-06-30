""" foam_dicts

Managment of OpenFOAM dicts and Fields

"""

import simphonyfoaminterfaceII as foamface
from simphony.core.cuba import CUBA

dictionaryMaps = \
    {'pimpleFoam':
        {'transportProperties':
            {'transportModel': 'Newtonian',
             'nu nu     [ 0 2 -1 0 0 0 0 ]': '0.0'},
         'turbulenceProperties':
            {'simulationType': 'laminar'},
         'RASProperties':
            {'RASModel': 'laminar',
             'turbulence': 'off',
             'printCoeffs': 'on'},
         'controlDict':
            {'application': 'pimpleFoam',
             'startFrom': 'startTime',
             'startTime': '0',
             'stopAt': 'endTime',
             'endTime': '1',
             'deltaT': '1',
             'writeControl': 'timeStep',
             'writeInterval': '10000',
             'purgeWrite': '0',
             'writeFormat': 'ascii',
             'writePrecision': '6',
             'writeCompression': 'no',
             'timeFormat': 'general',
             'runTimeModifiable': 'yes'}}}

dictionaryTemplates = {'pimpleFoam':
                       {'transportProperties':
                        """
transportModel Newtonian;

nu nu     [ 0 2 -1 0 0 0 0 ]     0.0 ;

                        """,
                        'turbulenceProperties':
                        """
simulationType laminar;
                        """,
                        'RASProperties':
                        """
RASModel            laminar;

turbulence          off;

printCoeffs         on;
                        """,
                        'fvSolution':
                        """

solvers
{
"(p|pFinal)"
{
solver          GAMG;
tolerance       1e-07;
relTol          0.001;
smoother        GaussSeidel;
nPreSweeps      0;
nPostSweeps     2;
cacheAgglomeration true;
nCellsInCoarsestLevel 10;
agglomerator    faceAreaPair;
mergeLevels     1;
}

"(U|UFinal)"
{
solver          smoothSolver;
smoother        GaussSeidel;
nSweeps         2;
tolerance       1e-08;
relTol          0.01;
}

nuTilda
{
solver          smoothSolver;
smoother        GaussSeidel;
nSweeps         2;
tolerance       1e-08;
relTol          0.1;
}
}

PIMPLE
{
nOuterCorrectors    2;
nCorrectors         3;
nNonOrthogonalCorrectors 2;
pRefCell        0;
pRefValue       0;
}

relaxationFactors
{
fields
{
p               0.3;
}
equations
{
U               0.7;
nuTilda         0.7;
}
}
                        """,
                        'fvSchemes':
                        """
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)     Gauss upwind;
    div((nuEff*dev(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

fluxRequired
{
    default         no;
    pcorr           ;
    p               ;
}

                        """}}


def parse_map(mapContent):
    """ mapping to a string

    """
    result = ""
    for i, j in mapContent.iteritems():
        result += "%s\t%s; \n" % (i, j)

    return result


def modifyNumerics(mesh, SP):
    """ Modifies the numerical parameters of the simulation
    """
    fileContent = dictionaryTemplates['pimpleFoam']
    mapContent = dictionaryMaps['pimpleFoam']

    fvSchemesDict = fileContent['fvSchemes']
    fvSolutionDict = fileContent['fvSolution']

    nOfTimeSteps = SP[CUBA.NUMBER_OF_TIME_STEPS]
    deltaT = SP[CUBA.TIME_STEP]
    endTime = nOfTimeSteps*deltaT + mesh._time

    mapContent['controlDict']['startTime'] = str(mesh._time)
    mapContent['controlDict']['deltaT'] = str(deltaT)
    mapContent['controlDict']['endTime'] = str(endTime)

    controlDict = parse_map(mapContent['controlDict'])

    foamface.modifyNumerics(mesh.name, fvSchemesDict, fvSolutionDict,
                            controlDict)


def modifyFields(mesh, BC):
    """ Modifies the internal fields and boundary conditions
    """
    dimensionset = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    nCells = foamface.getCellCount(mesh.name)
    p_values = [0.0 for item in range(nCells)]
    U_values = [(0.0, 0.0, 0.0) for item in range(nCells)]
    for cell in mesh.iter_cells():
        p_values[mesh._uuidToFoamLabel[cell.uid]] = \
            cell.data[CUBA.PRESSURE]
        U_values[mesh._uuidToFoamLabel[cell.uid]] = \
            cell.data[CUBA.VELOCITY]
    foamface.setAllCellData(mesh.name, "p", dimensionset, p_values)
    foamface.setAllCellVectorData(mesh.name, "U", dimensionset, U_values)

#       Refresh boundary conditions
    velocityBCs = BC[CUBA.VELOCITY]
    myDict = ""
    for boundary in velocityBCs:
        myDict = myDict + str(boundary) + "\n{\n"
        if velocityBCs[boundary] == "zeroGradient":
            myDict = myDict + "\t type \t zeroGradient;\n"
            myDict = myDict + "\t value \t uniform (0 0 0);\n"
        elif velocityBCs[boundary] == "empty":
            myDict = myDict + "\t type \t empty;\n"
        else:
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform (" \
                            + str(velocityBCs[boundary][0]) + " " \
                            + str(velocityBCs[boundary][1]) + " " \
                            + str(velocityBCs[boundary][2]) + ");\n"
        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, "U", myDict)

    pressureBCs = BC[CUBA.PRESSURE]
    myDict = ""
    for boundary in pressureBCs:
        myDict = myDict + str(boundary) + "\n{\n"
        if pressureBCs[boundary] == "zeroGradient":
            myDict = myDict + "\t type \t zeroGradient;\n"
            myDict = myDict + "\t value \t uniform 0;\n"
        elif pressureBCs[boundary] == "fixedFluxPressure":
            myDict = myDict + "\t type \t fixedFluxPressure;\n"
            myDict = myDict + "\t value \t uniform 0;"
        elif pressureBCs[boundary] == "empty":
            myDict = myDict + "\t type \t empty;\n"
        else:
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform " \
                            + str(pressureBCs[boundary]) + ";\n"
        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, "p", myDict)
