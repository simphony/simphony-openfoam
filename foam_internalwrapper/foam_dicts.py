""" foam_dicts

Managment of OpenFOAM dicts and Fields

"""

import simphonyfoaminterface as foamface
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
             'runTimeModifiable': 'yes'}},
     'driftFluxSimphonyFoam':
         {'controlDict':
            {'application': 'driftFluxSimphonyFoam',
             'startFrom': 'startTime',
             'startTime': '0',
             'stopAt': 'endTime',
             'endTime': '0.1',
             'deltaT': '0.1',
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

                        """},
                       'driftFluxSimphonyFoam':
                       {'transportProperties':
                        """
phases (heavy light);

heavy
{
    transportModel  BinghamPlastic;

    "(plastic|BinghamPlastic)Coeffs"
    {
        coeff       0.00023143;
        exponent    179.26;

        BinghamCoeff    0.0005966;
        BinghamExponent 1050.8;
        BinghamOffset   0;

        muMax       10;
    }

    rho         1996;
}

light
{
    transportModel  Newtonian;

    nu              1.7871e-06;
    rho             996;
}

relativeVelocityModel simple;

"(simple|general)Coeffs"
{
    V0              (0 -0.002198 0);
    a               285.84;
    a1              0.1;
    residualAlpha   0;
}
                        """,
                        'fvSolution':
                        """
solvers
{
    "alpha.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 1;

        MULESCorr       yes;
        nLimiterIter    3;
        alphaApplyPrevCorr  yes;

        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0;
        minIter         1;
    }

    "alpha.*Diffusion"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-6;
        relTol          0;
        minIter         1;
    }

    p_rgh
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.01;
        smoother        GaussSeidel;
        cacheAgglomeration true;
        nCellsInCoarsestLevel 20;
        agglomerator    faceAreaPair;
        mergeLevels     1;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|k|epsilon)"
    {
        solver          PBiCG;
        preconditioner  DILU;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
        minIter         1;
    }

    "(U|k|epsilon)Final"
    {
        $k;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor no;
    nCorrectors     3;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        ".*"           1;
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
    default             Gauss linear;
    div(SigmaPrime)       Gauss linear;
    div((mu*dev(T(grad(U))))) Gauss linear;
    div(rhoPhi,U)       Gauss linearUpwind grad(U);
    div(tauDm)          Gauss linear;
    "div\(phi,alpha.*\)" Gauss vanLeer;
    "div\(phirb,alpha.*\)" Gauss linear;
    div(rhoPhi,k)       Gauss limitedLinear 1;
    div(rhoPhi,epsilon) Gauss limitedLinear 1;

    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
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
    p_rgh;
    "alpha.*";
}
                        """}}


def parse_map(mapContent):
    """ mapping to a string

    """
    result = ""
    for i, j in mapContent.iteritems():
        result += "%s\t%s; \n" % (i, j)

    return result


def modifyNumerics(mesh, SP, solver='pimpleFoam'):
    """ Modifies the numerical parameters of the simulation
    """
    fileContent = dictionaryTemplates[solver]
    mapContent = dictionaryMaps[solver]

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


def modifyFields(mesh, BC, solver='pimpleFoam'):
    """ Modifies the internal fields and boundary conditions
    """
    if solver=='pimpleFoam':
        name_pressure = 'p'
        ID_pressure = CUBA.PRESSURE
    elif solver == 'driftFluxSimphonyFoam':
        name_pressure = 'p_rgh'
        ID_pressure = CUBA.CONCENTRATION

    nCells = foamface.getCellCount(mesh.name)
    p_values = [0.0 for item in range(nCells)]
    U_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
    alpha_values = [0.0 for item in range(nCells)]
    for cell in mesh.iter_cells():
        p_values[mesh._uuidToFoamLabel[cell.uid]] = \
            cell.data[ID_pressure]
        U_values[mesh._uuidToFoamLabel[cell.uid]] = \
            list(cell.data[CUBA.VELOCITY])
        if solver == 'driftFluxSimphonyFoam':
            alpha_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.VOLUME_FRACTION]
    foamface.setAllCellData(mesh.name, name_pressure , p_values)
    foamface.setAllCellVectorData(mesh.name, "U", U_values)
    if solver == 'driftFluxSimphonyFoam':
        foamface.setAllCellData(mesh.name, "alpha1", alpha_values)

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
        elif velocityBCs[boundary] == "slip":
            myDict = myDict + "\t type \t slip;\n"
        else:
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform (" \
                            + str(velocityBCs[boundary][0]) + " " \
                            + str(velocityBCs[boundary][1]) + " " \
                            + str(velocityBCs[boundary][2]) + ");\n"
        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, "U", myDict)

    pressureBCs = BC[ID_pressure]
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

    foamface.setBC(mesh.name, name_pressure, myDict)
    
    
    if solver=='driftFluxSimphonyFoam':
        volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
        myDict = ""
        for boundary in volumeFractionBCs:
            myDict = myDict + str(boundary) + "\n{\n"
            if volumeFractionBCs[boundary] == "zeroGradient":
                myDict = myDict + "\t type \t zeroGradient;\n"
                myDict = myDict + "\t value \t uniform 0;\n"
            elif isinstance(volumeFractionBCs[boundary], tuple) and volumeFractionBCs[boundary][0] == "inletOutlet":
                myDict = myDict + "\t type \t inletOutlet;\n"
                myDict = myDict + "\t inletValue \t uniform %s;\n"%str(volumeFractionBCs[boundary][1])
            elif volumeFractionBCs[boundary] == "empty":
                 myDict = myDict + "\t type \t empty;\n"
            else:
                myDict = myDict + "\t type \t fixedValue;\n"
                myDict = myDict + "\t value \t uniform " \
                                + str(volumeFractionBCs[boundary]) + ";\n"
            myDict = myDict + "}\n"

        foamface.setBC(mesh.name, "alpha1", myDict)
        
