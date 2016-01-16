""" foam_dicts

Managment of OpenFOAM dicts and Fields

"""

import simphonyfoaminterface as foamface
from simphony.core.cuba import CUBA
from foam_controlwrapper.cuba_extension import CUBAExt

dictionaryMaps = \
    {'pimpleFoam':
        {'transportProperties':
            {'transportModel': 'Newtonian',
             'nu nu     [ 0 2 -1 0 0 0 0 ]': '0.0',
             'g': '( 0 0 0)'},
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
         {
            'controlDict':
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
                 'runTimeModifiable': 'yes'},
            'transportProperties':
                {'g': '( 0 0 0)',
                 'phases': '(phase1 phase2)',
                 'phase1':
                     {'transportModel': 'dummyViscosity',
                      'BinghamPlasticCoeffs':
                          {'coeff': '0.00023143',
                           'exponent': '179.26',
                           'BinghamCoeff': '0.0005966',
                           'BinghamExponent': '1050.8',
                           'BinghamOffset': '0',
                           'muMax': '10'},
                      'plasticCoeffs':
                          {'coeff': '0.00023143',
                           'exponent': '179.26',
                           'BinghamCoeff': '0.0005966',
                           'BinghamExponent': '1050.8',
                           'BinghamOffset': '0',
                           'muMax': '10'},
                      'nu': '1e-4',
                      'rho': '1996'},
                 'phase2':
                     {'transportModel': 'Newtonian',
                      'nu': '1.7871e-06',
                      'rho': '996'},
                 'stressModel': 'standard',
                 'relativeVelocityModel': 'simple',
                 'simpleCoeffs':
                     {'V0': '(0 -0.002198 0)',
                      'a': '285.84',
                      'a1': '0.1',
                      'residualAlpha': '0'},
                 'generalCoeffs':
                     {'V0': '(0 -0.002198 0)',
                      'a': '285.84',
                      'a1': '0.1',
                      'residualAlpha': '0'},
                 'fromMesoscaleCoeffs':
                     {}
                 }}}

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
    p
    {
        solver           GAMG;
        tolerance        1e-7;
        relTol           0.01;

        smoother         DICGaussSeidel;

        cacheAgglomeration true;
        nCellsInCoarsestLevel 10;
        agglomerator     faceAreaPair;
        mergeLevels      1;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    "(U|k|epsilon)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

    "(U|k|epsilon)Final"
    {
        $U;
        relTol          0;
    }
}

PIMPLE
{
    nNonOrthogonalCorrectors 0;
    nCorrectors         2;
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
    grad(p)         Gauss linear;
    grad(U)         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,R)      bounded Gauss upwind;
    div(R)          Gauss linear;
    div(phi,nuTilda) bounded Gauss upwind;
    div((nuEff*dev(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         none;
    laplacian(nuEff,U) Gauss linear corrected;
    laplacian(rAUf,p)  Gauss linear corrected;
    laplacian(DkEff,k) Gauss linear corrected;
    laplacian(DepsilonEff,epsilon) Gauss linear corrected;
    laplacian(DREff,R) Gauss linear corrected;
    laplacian(DnuTildaEff,nuTilda) Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
    interpolate(U)  linear;
}

snGradSchemes
{
    default         corrected;
}

fluxRequired
{
    default         no;
    p               ;
}

                        """},
                       'driftFluxSimphonyFoam':
                       {'transportProperties':
                        """
phases (phase1 phase2);

phase1
{
    transportModel  dummyViscosity;
    nu              1e-04;
    rho             1996;
}

phase2
{
    transportModel  Newtonian;

    nu              1.7871e-06;
    rho             996;
}

stressModel standard;
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
        if isinstance(j, dict):
            result += "%s{\n" % i
            result += parse_map(j)
            result += "}\n"
        else:
            result += "%s\t%s; \n" % (i, j)

    return result


def modifyNumerics(mesh, SP, SPExt, solver='pimpleFoam'):
    """ Modifies the numerical parameters of the simulation
    """
    fileContent = dictionaryTemplates[solver]
    mapContent = dictionaryMaps[solver]

    fvSchemesDict = fileContent['fvSchemes']
    fvSolutionDict = fileContent['fvSolution']
    transportPropertiesDict = fileContent['transportProperties']

    nOfTimeSteps = SP[CUBA.NUMBER_OF_TIME_STEPS]
    deltaT = SP[CUBA.TIME_STEP]
    endTime = nOfTimeSteps*deltaT + mesh._time

    mapContent['controlDict']['startTime'] = str(mesh._time)
    mapContent['controlDict']['deltaT'] = str(deltaT)
    mapContent['controlDict']['endTime'] = str(endTime)

    controlDict = parse_map(mapContent['controlDict'])

    if solver == 'pimpleFoam':
        control = mapContent['transportProperties']
        control['nu nu     [ 0 2 -1 0 0 0 0 ]'] = \
            SP[CUBA.DYNAMIC_VISCOSITY]/SP[CUBA.DENSITY]
        transportPropertiesDict = parse_map(control)
    elif solver == 'driftFluxSimphonyFoam':
        control = mapContent['transportProperties']
        density = SP[CUBA.DENSITY]
        viscosity = SP[CUBA.DYNAMIC_VISCOSITY]
        # phase1 uses mixtureViscosity models
        viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL][
            SPExt[CUBAExt.PHASE_LIST][0]]
        vmc = viscosity_model + 'Coeffs'
        viscosity_model_coeffs =\
            SPExt[CUBAExt.VISCOSITY_MODEL_COEFFS][viscosity_model]
        for coeff in viscosity_model_coeffs:
            control['phase1'][vmc][coeff] = viscosity_model_coeffs[coeff]
        control['phase1']['rho'] = density[SPExt[CUBAExt.PHASE_LIST][0]]

        viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL][
            SPExt[CUBAExt.PHASE_LIST][1]]
        if viscosity_model == 'Newtonian':
            control["phase2"]["transportModel"] = viscosity_model
            control['phase2']['nu'] = \
                viscosity[SPExt[CUBAExt.PHASE_LIST][1]] / \
                density[SPExt[CUBAExt.PHASE_LIST][1]]
            control['phase2']['rho'] = \
                density[SPExt[CUBAExt.PHASE_LIST][1]]
        else:
            control["phase2"]["transportModel"] = viscosity_model
            vmc = viscosity_model + 'Coeffs'
            viscosity_model_coeffs =\
                SPExt[CUBAExt.VISCOSITY_MODEL_COEFFS][viscosity_model]
            for coeff in viscosity_model_coeffs:
                control["phase2"][vmc][coeff] =\
                    viscosity_model_coeffs[coeff]

        if CUBAExt.STRESS_MODEL in SPExt:
                    control["stressModel"] =\
                        SPExt[CUBAExt.STRESS_MODEL]
        else:
            error_str = "Stress model not specified"
            raise ValueError(error_str)

        if CUBAExt.RELATIVE_VELOCITY_MODEL in SPExt:
            control["relativeVelocityModel"] =\
                SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL]
        else:
            error_str = "Relative velocity model not specified"
            raise ValueError(error_str)
        if SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL] == "simple"\
                or SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL] == "general":
            relVelModelCoeffs = SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL]+'Coeffs'
            if CUBAExt.RELATIVE_VELOCITY_MODEL_COEFFS in SPExt:
                coeffs = SPExt[CUBAExt.RELATIVE_VELOCITY_MODEL_COEFFS]
                V0 = coeffs["V0"]
                control[relVelModelCoeffs]["V0"] = "( " \
                    + str(V0[0]) + " " + str(V0[1]) + " " + str(V0[2]) + " )"
                control[relVelModelCoeffs]["a"] = coeffs["a"]
                control[relVelModelCoeffs]["a1"] = coeffs["a1"]
                control[relVelModelCoeffs]["residualAlpha"] =\
                    coeffs["residualAlpha"]
            if CUBAExt.EXTERNAL_BODY_FORCE_MODEL in SPExt:
                if SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL] == 'gravitation':
                    g = SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL_COEFFS]["g"]
                    control["g"] =\
                        "( " \
                        + str(g[0]) + " " + str(g[1]) + " " + str(g[2]) + " )"

        transportPropertiesDict = parse_map(control)

    foamface.modifyNumerics(mesh.name, fvSchemesDict, fvSolutionDict,
                            controlDict, transportPropertiesDict)


def modifyFields(mesh, BC, solver='pimpleFoam'):
    """ Modifies the internal fields and boundary conditions
    """
    if solver == 'pimpleFoam':
        name_pressure = 'p'
        ID_pressure = CUBA.PRESSURE
    elif solver == 'driftFluxSimphonyFoam':
        name_pressure = 'p_rgh'
        ID_pressure = CUBA.CONCENTRATION

    nCells = foamface.getCellCount(mesh.name)
    p_values = [0.0 for item in range(nCells)]
    U_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
    if solver == 'driftFluxSimphonyFoam':
        alpha_values = [0.0 for item in range(nCells)]
        vdj_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
        mu_sigma_values = [[0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0] for item in range(nCells)]
    for cell in mesh.iter_cells():
        p_values[mesh._uuidToFoamLabel[cell.uid]] = \
            cell.data[ID_pressure]
        U_values[mesh._uuidToFoamLabel[cell.uid]] = \
            list(cell.data[CUBA.VELOCITY])
        if solver == 'driftFluxSimphonyFoam':
            alpha_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.VOLUME_FRACTION]
            vdj_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.ANGULAR_VELOCITY]
            mu_sigma_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.ACCELERATION]
    foamface.setAllCellData(mesh.name, name_pressure, p_values)
    foamface.setAllCellVectorData(mesh.name, "U", U_values)
    if solver == 'driftFluxSimphonyFoam':
        foamface.setAllCellData(mesh.name, "alpha.phase1", alpha_values)

#       Refresh boundary conditions
    velocityBCs = BC[CUBA.VELOCITY]
    myDict = ""

    for boundary in velocityBCs:
        patch = velocityBCs[boundary]
        myDict = myDict + str(boundary) + "\n{\n"
        if patch == "zeroGradient":
            myDict = myDict + "\t type \t zeroGradient;\n"
            myDict = myDict + "\t value \t uniform (0 0 0);\n"
        elif patch == "empty":
            myDict = myDict + "\t type \t empty;\n"
        elif patch == "slip":
            myDict = myDict + "\t type \t slip;\n"
        elif isinstance(patch, tuple) and patch[0] == "pressureIOVelocity":
            myDict = myDict + "\t type \t pressureInletOutletVelocity;\n"
            myDict = myDict + "\t value \t uniform (" \
                + str(patch[1][0]) + " " \
                + str(patch[1][1]) + " " \
                + str(patch[1][2]) + ");\n"
        elif isinstance(patch, tuple) and patch[0] == "fixedValue":
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform (" \
                + str(patch[1][0]) + " " \
                + str(patch[1][1]) + " " \
                + str(patch[1][2]) + ");\n"
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
        elif isinstance(pressureBCs[boundary], tuple) and\
                pressureBCs[boundary][0] == "fixedValue":
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform " \
                            + str(pressureBCs[boundary][1]) + ";\n"
        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, name_pressure, myDict)

    if solver == 'driftFluxSimphonyFoam':
        volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
        myDict = ""
        for boundary in volumeFractionBCs:
            myDict = myDict + str(boundary) + "\n{\n"
            if volumeFractionBCs[boundary] == "zeroGradient":
                myDict = myDict + "\t type \t zeroGradient;\n"
                myDict = myDict + "\t value \t uniform 0;\n"
            elif volumeFractionBCs[boundary] == "empty":
                myDict = myDict + "\t type \t empty;\n"
            elif isinstance(volumeFractionBCs[boundary], tuple):
                if volumeFractionBCs[boundary][0] == "inletOutlet":
                    myDict = myDict + "\t type \t inletOutlet;\n"
                    myDict = myDict + "\t inletValue \t uniform " \
                        + str(volumeFractionBCs[boundary][1]) + ";\n"
                elif volumeFractionBCs[boundary][0] == "fixedValue":
                    myDict = myDict + "\t type \t fixedValue;\n"
                    myDict = myDict + "\t value \t uniform " \
                        + str(volumeFractionBCs[boundary][1]) + ";\n"
            myDict = myDict + "}\n"

        foamface.setBC(mesh.name, "alpha.phase1", myDict)
