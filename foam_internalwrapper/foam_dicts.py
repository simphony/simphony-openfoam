""" foam_dicts

Management of OpenFOAM dicts and Fields

"""
import os

import simphonyfoaminterface as foamface
from simphony.core.cuba import CUBA
from foam_controlwrapper.cuba_extension import CUBAExt
from foam_controlwrapper.foam_variables import (dataDimensionMap, dataKeyMap)

userLibs = """
("libshearStressPowerLawSlipVelocity.so"  "libtwoPhaseProperties.so")
"""
userLibsIO = """
("libshearStressPowerLawSlipVelocityIO.so" "libtwoPhaseProperties.so")
"""

dictionaryMaps = \
    {'pimpleFoam':
        {'transportProperties':
            {'transportModel': 'Newtonian',
             'nu nu     [ 0 2 -1 0 0 0 0 ]': '0.0',
             'CrossPowerLawCoeffs':
                 {
                     'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                     'nuInf': '           nuInf [ 0 2 -1 0 0 0 0 ] {}',
                     'm': '               m [ 0 0 1 0 0 0 0 ] {}',
                     'n': '               n [ 0 0 0 0 0 0 0 ] {}'
                     },
             'BirdCarreauCoeffs':
                 {
                     'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                     'nuInf': '           nuInf [ 0 2 -1 0 0 0 0 ] {}',
                     'k': '               k [ 0 0 1 0 0 0 0 ] {}',
                     'n': '               n [ 0 0 0 0 0 0 0 ] {}'
                     },
             'HerschelBulkleyCoeffs':
                 {
                     'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                     'tau0': '            tau0 [ 0 2 -2 0 0 0 0 ] {}',
                     'k': '               k [ 0 2 -1 0 0 0 0 ] {}',
                     'n': '               n [ 0 0 0 0 0 0 0 ] {}'
                     }
             },
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
              'runTimeModifiable': 'yes'}
         },
     'interFoam':
         {'transportProperties':
             {'phases': '(phase1 phase2)',
              'phase1':
                  {'transportModel': 'Newtonian',
                   'nu              nu [ 0 2 -1 0 0 0 0 ]': '0.0',
                   'rho             rho [ 1 -3 0 0 0 0 0 ]': '0.0',
                   'CrossPowerLawCoeffs':
                       {
                        'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                        'nuInf': '           nuInf [ 0 2 -1 0 0 0 0 ] {}',
                        'm': '               m [ 0 0 1 0 0 0 0 ] {}',
                        'n': '               n [ 0 0 0 0 0 0 0 ] {}'
                           },
                   'BirdCarreauCoeffs':
                       {
                        'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                        'nuInf': '           nuInf [ 0 2 -1 0 0 0 0 ] {}',
                        'k': '               k [ 0 0 1 0 0 0 0 ] {}',
                        'n': '               n [ 0 0 0 0 0 0 0 ] {}'
                           },
                   'HerschelBulkleyCoeffs':
                       {
                        'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                        'tau0': '            tau0 [ 0 2 -2 0 0 0 0 ] {}',
                        'k': '               k [ 0 2 -1 0 0 0 0 ] {}',
                        'n': '               n [ 0 0 0 0 0 0 0 ] {}'
                           }
                   },
              'phase2':
                  {'transportModel': 'Newtonian',
                   'nu              nu [ 0 2 -1 0 0 0 0 ]': '0.0',
                   'rho             rho [ 1 -3 0 0 0 0 0 ]': '0.0',
                   'CrossPowerLawCoeffs': {
                       'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                       'nuInf': '           nuInf [ 0 2 -1 0 0 0 0 ] {}',
                       'm': '               m [ 0 0 1 0 0 0 0 ] {}',
                       'n': '               n [ 0 0 0 0 0 0 0 ] {}'},
                   'BirdCarreauCoeffs': {
                       'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                       'nuInf': '           nuInf [ 0 2 -1 0 0 0 0 ] {}',
                       'k': '               k [ 0 0 1 0 0 0 0 ] {}',
                       'n': '               n [ 0 0 0 0 0 0 0 ] {}'},
                   'HerschelBulkleyCoeffs': {
                       'nu0': '             nu0 [ 0 2 -1 0 0 0 0 ] {}',
                       'tau0': '            tau0 [ 0 2 -2 0 0 0 0 ] {}',
                       'k': '               k [ 0 2 -1 0 0 0 0 ] {}',
                       'n': '               n [ 0 0 0 0 0 0 0 ] {}'}
                   },
              'sigma           sigma [ 1 0 -2 0 0 0 0 ]': '0.0'
              },
          'turbulenceProperties':
              {'simulationType': 'laminar'},
          'RASProperties':
              {'RASModel': 'laminar',
               'turbulence': 'off',
               'printCoeffs': 'on'},
          'controlDict':
              {'application': 'interFoam',
               'startFrom': 'startTime',
               'startTime': '0',
               'stopAt': 'endTime',
               'endTime': '1',
               'deltaT': '1',
               'writeControl': 'adjustableRunTime',
               'writeInterval': '1',
               'purgeWrite': '0',
               'writeFormat': 'ascii',
               'writePrecision': '6',
               'writeCompression': 'no',
               'timeFormat': 'general',
               'runTimeModifiable': 'yes',
               'adjustTimeStep': 'on',
               'maxCo': '1',
               'maxAlphaCo': '1',
               'maxDeltaT': '1'
               },
          'g':
              {'dimensions': '[0 1 -2 0 0 0 0]',
               'value':  '( 0 0 0 )'
               }
          },
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
               'runTimeModifiable': 'yes',
               'adjustTimeStep':  'on',
               'maxCo': '1',
               'maxDeltaT': '1'
               },
          'transportProperties':
              {'phases': '(phase1 phase2)',
               'phase1':
                   {'transportModel': 'dummyViscosity',
                    'BinghamPlasticCoeffs':
                        {'coeff': '0.0',
                         'exponent': '0.0',
                         'BinghamCoeff': '0.0',
                         'BinghamExponent': '0.0',
                         'BinghamOffset': '0',
                         'muMax': '0'},
                    'plasticCoeffs':
                        {'coeff': '0.0',
                         'exponent': '0.0',
                         'BinghamCoeff': '0.0',
                         'BinghamExponent': '0.0',
                         'BinghamOffset': '0',
                         'muMax': '0'},
                    'nu': '0.0',
                    'rho': '0.0'},
               'phase2':
                   {'transportModel': 'Newtonian',
                    'nu': '0.0',
                    'rho': '0.0'},
               'stressModel': 'standard',
               'relativeVelocityModel': 'simple',
               'simpleCoeffs':
                   {'V0': '(0 0.0 0)',
                    'a': '0.0',
                    'a1': '0.0',
                    'residualAlpha': '0'},
               'generalCoeffs':
                   {'V0': '(0 0.0 0)',
                    'a': '0.0',
                    'a1': '0.0',
                    'residualAlpha': '0'},
               'fromMesoscaleCoeffs':
                   {}
               },
          'g':
              {'dimensions': '[0 1 -2 0 0 0 0]',
               'value':  '( 0 0 0 )'
               }
          },
     'blockMesh':
         {'blockMeshDict':
              {'convertToMeters': 1,
               'vertices': (),
               'blocks': (),
               'edges': (),
               'boundary': (),
               'mergePatchPairs': ()},
          'controlDict':
              {'application': 'blockMesh',
               'startFrom': 'startTime',
               'startTime': '0',
               'stopAt': 'endTime',
               'endTime': '1',
               'deltaT': '1',
               'writeControl': 'timeStep',
               'writeInterval': '1'
               },
          'fvSchemes':
              {'divSchemes': {'default': 'none'},
               'gradSchemes': {'default': 'Gauss linear'},
               'laplacianSchemes': {'default': 'none'}
               },
          'fvSolution':
              {}
          }
     }

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
    div(phi,U)      bounded Gauss linearUpwind limited;
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
                       'interFoam':
                       {'transportProperties':
                            """
phases (phase1 phase2);


phase1
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 0.0;
    rho             rho [ 1 -3 0 0 0 0 0 ] 0.0;
 }
phase2
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 0.0;
    rho             rho [ 1 -3 0 0 0 0 0 ] 0.0;
 }
sigma           sigma [ 1 0 -2 0 0 0 0 ] 0.0;

                        """,
                        'g':
                        """
dimensions      [0 1 -2 0 0 0 0];
value           ( 0 0 0 );
                        """,
                        'fvSolution':
                        """

solvers
{
    "alpha.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 1;
        cAlpha          1;

        MULESCorr       yes;
        nLimiterIter    3;

        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }

    pcorr
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-5;
        relTol          0;
    }

    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-07;
        relTol          0.05;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    fields
    {
    }
    equations
    {
        ".*" 1;
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
    div(rhoPhi,U)  Gauss upwind;
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div((muEff*dev(T(grad(U))))) Gauss linear;
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
    pcorr;
    alpha.phase1;
}


                        """},
                       'driftFluxSimphonyFoam':
                       {'transportProperties':
                            """
phases (phase1 phase2);

phase1
{
    transportModel  dummyViscosity;
    nu              0.0;
    rho             0.0;
    BinghamPlasticCoeffs
    {
        coeff       0.0;
        exponent    0.0;

        BinghamCoeff    0.0;
        BinghamExponent 0.0;
        BinghamOffset   0;

        muMax       0;
    }
    plasticCoeffs
    {
        coeff       0.0;
        exponent    0.0;

        BinghamCoeff    0.0;
        BinghamExponent 0.0;
        BinghamOffset   0;

        muMax       0;
    }

}

phase2
{
    transportModel  Newtonian;

    nu              0.0;
    rho             0.0;
}

stressModel standard;
relativeVelocityModel simple;

simpleCoeffs
{
    V0              (0 0 0);
    a               0.0;
    a1              0.0;
    residualAlpha   0;
}
generalCoeffs
{
    V0              (0 0 0);
    a               0.0;
    a1              0.0;
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
                        """},
                       'blockMesh':
                       {'blockMeshDict':
                            """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (0 0 0)
    (30 0 0)
    (30 5 0)
    (0 5 0)
    (0 0 0.1)
    (30 0 0.1)
    (30 5 0.1)
    (0 5 0.1)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (500 20 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    walls
    {
        type wall;
        faces
        (
            (3 7 6 2)
            (1 5 4 0)
        );
    }
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (2 6 5 1)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //


"""}}

head = \
    """
FoamFile
{{
    version     {version};
    format      ascii;
    class       {foamclass};
    location    {location};
    object      {foamobject};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    """


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


def modifyNumerics(mesh, SP, SPExt, solver='pimpleFoam', io=False):
    """ Modifies the numerical parameters of the simulation

        Parameters
        ----------
        mesh : FoamMesh
            mesh to be added.
        SP : DataContainer
            System Parameters
        SPExt : dictionary
            Extension to System Parameters
        solver : str
            solver name
        io : bool
            whether or not to write dictionaries to disk

    """

    fileContent = dictionaryTemplates[solver]
    mapContent = get_dictionary_maps(solver, io)

    fvSchemesDict = fileContent['fvSchemes']
    fvSolutionDict = fileContent['fvSolution']
    transportPropertiesDict = fileContent['transportProperties']

    nOfTimeSteps = SP[CUBA.NUMBER_OF_TIME_STEPS]
    deltaT = SP[CUBA.TIME_STEP]
    interval = nOfTimeSteps*deltaT
    endTime = interval + mesh._time
    if CUBAExt.MAX_COURANT_NUMBER in SPExt:
        mapContent['controlDict']['maxCo'] = SPExt[CUBAExt.MAX_COURANT_NUMBER]
        mapContent['controlDict']['maxAlphaCo'] =\
            SPExt[CUBAExt.MAX_COURANT_NUMBER]
    mapContent['controlDict']['startTime'] = str(mesh._time)
    mapContent['controlDict']['deltaT'] = str(deltaT)
    mapContent['controlDict']['endTime'] = str(endTime)
    mapContent['controlDict']['maxDeltaT'] = str(deltaT)
    if io:
        mapContent['controlDict']['writeControl'] = 'adjustableRunTime'
        mapContent['controlDict']['writeInterval'] = str(interval)

    controlDict = parse_map(mapContent['controlDict'])

    if solver == 'pimpleFoam':
        control = mapContent['transportProperties']
        if CUBAExt.VISCOSITY_MODEL in SPExt:
            viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL]
            control['transportModel'] = viscosity_model
            vmc = viscosity_model + 'Coeffs'
            viscosity_model_coeffs =\
                SPExt[CUBAExt.VISCOSITY_MODEL_COEFFS][viscosity_model]
            for coeff in viscosity_model_coeffs:
                control[vmc][coeff] =\
                    control[vmc][coeff].format(viscosity_model_coeffs[coeff])
        else:
            control['nu nu     [ 0 2 -1 0 0 0 0 ]'] = \
                SP[CUBA.DYNAMIC_VISCOSITY]/SP[CUBA.DENSITY]

        transportPropertiesDict = parse_map(control)
    elif solver == 'interFoam':
        density = SP[CUBA.DENSITY]
        viscosity = SP[CUBA.DYNAMIC_VISCOSITY]
        control = mapContent['transportProperties']
        for i in range(2):
            if CUBAExt.VISCOSITY_MODEL in SPExt:
                viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL][
                    SPExt[CUBAExt.PHASE_LIST][i]]
            else:
                viscosity_model = 'Newtonian'
            phase_name = 'phase' + str(i + 1)
            if viscosity_model == 'Newtonian':
                control[phase_name]['nu              nu [ 0 2 -1 0 0 0 0 ]']\
                    = viscosity[SPExt[CUBAExt.PHASE_LIST][i]] / \
                    density[SPExt[CUBAExt.PHASE_LIST][i]]
                control[phase_name]['rho             rho [ 1 -3 0 0 0 0 0 ]']\
                    = density[SPExt[CUBAExt.PHASE_LIST][i]]
            else:
                control[phase_name]['transportModel'] = viscosity_model
                vmc = viscosity_model + 'Coeffs'
                viscosity_model_coeffs =\
                    SPExt[CUBAExt.VISCOSITY_MODEL_COEFFS][viscosity_model]
                for coeff in viscosity_model_coeffs:
                    control[phase_name][vmc][coeff] =\
                        control[phase_name][vmc][coeff].\
                        format(viscosity_model_coeffs[coeff])

        phases = SPExt[CUBAExt.PHASE_LIST]
        if CUBAExt.SURFACE_TENSION in SPExt:
            if phases in SPExt[CUBAExt.SURFACE_TENSION]:
                control['sigma           sigma [ 1 0 -2 0 0 0 0 ]'] = \
                    SPExt[CUBAExt.SURFACE_TENSION][phases]
            else:
                phases = SPExt[CUBAExt.PHASE_LIST]
                if phases in SPExt[CUBAExt.SURFACE_TENSION]:
                    control['sigma           sigma [ 1 0 -2 0 0 0 0 ]'] = \
                        SPExt[CUBAExt.SURFACE_TENSION][phases]
                else:
                    error_str = "Surface tension not specified"
                    raise ValueError(error_str)
        else:
            error_str = "Surface tension not specified"
            raise ValueError(error_str)
        transportPropertiesDict = parse_map(control)
    elif solver == 'driftFluxSimphonyFoam':
        control = mapContent['transportProperties']
        density = SP[CUBA.DENSITY]
        viscosity = SP[CUBA.DYNAMIC_VISCOSITY]
        # phase1 uses mixtureViscosity models
        viscosity_model = SPExt[CUBAExt.VISCOSITY_MODEL][
            SPExt[CUBAExt.PHASE_LIST][0]]
        if viscosity_model == 'Newtonian':
            viscosity_model = 'dummyViscosity'
        control["phase1"]["transportModel"] = viscosity_model
        if viscosity_model != 'slurry' and viscosity_model != 'dummyViscosity':
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

        transportPropertiesDict = parse_map(control)
    foamface.modifyNumerics(mesh.name, fvSchemesDict, fvSolutionDict,
                            controlDict, transportPropertiesDict, int(io))
    if io:
        # write numerics to case directory
        foamface.writeNumerics(mesh.name)

    if solver == 'driftFluxSimphonyFoam' or solver == 'interFoam':
        if io:
            control = mapContent['g']
            if CUBAExt.EXTERNAL_BODY_FORCE_MODEL in SPExt:
                if SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL] == 'gravitation':
                    g = SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL_COEFFS]['g']
                    control['value'] =\
                        "( " \
                        + str(g[0]) + " " + str(g[1]) + " " + str(g[2]) + " )"
                    gravityDict = parse_map(control)
                    foamface.modifyDictionary(mesh.name, 'g', gravityDict)
                    foamface.writeDictionary(mesh.name, 'g', True)
        else:
            if CUBAExt.EXTERNAL_BODY_FORCE_MODEL in SPExt:
                if SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL] == 'gravitation':
                    g = SPExt[CUBAExt.EXTERNAL_BODY_FORCE_MODEL_COEFFS]['g']
                    foamface.modifyUniformVectorField(mesh.name, 'g', list(g))


def check_boundary_names(bc_names, boundary_names, cuba):
    if not set(bc_names).issubset(set(boundary_names)):
        error_str = "Boundary name(s) used in boundary conditions "
        error_str += "for {} does not exist.\nUsed name(s) are: {}\n"
        error_str += "Boundary names defined in the mesh are: {}\n"
        raise ValueError(
            error_str.format(cuba.name,
                             list(set(bc_names).difference(
                                 set(boundary_names))),
                             boundary_names))


def modifyFields(mesh, BC, solver='pimpleFoam'):
    """ Modifies the internal fields and boundary conditions

        Parameters
        ----------
        mesh : FoamMesh
            mesh to be added.
        BC : DataContainer
            Boundary Conditions
        solver : str
            solver name

    """

    if solver == 'pimpleFoam':
        name_pressure = 'p'
        ID_pressure = CUBA.PRESSURE
    elif solver == 'driftFluxSimphonyFoam' or solver == 'interFoam':
        name_pressure = 'p_rgh'
        ID_pressure = CUBA.DYNAMIC_PRESSURE
        # check that boundarynames match with mesh boundary names
        volumeFractionBCs = BC[CUBA.VOLUME_FRACTION]
        check_boundary_names(volumeFractionBCs, mesh._boundaries.keys(),
                             CUBA.VOLUME_FRACTION)

    # check that boundarynames match with mesh boundary names
    velocityBCs = BC[CUBA.VELOCITY]
    pressureBCs = BC[ID_pressure]
    check_boundary_names(velocityBCs, mesh._boundaries.keys(), CUBA.VELOCITY)
    check_boundary_names(pressureBCs, mesh._boundaries.keys(), ID_pressure)

    nCells = foamface.getCellCount(mesh.name)
    p_values = [0.0 for item in range(nCells)]
    U_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
    if solver == 'driftFluxSimphonyFoam' or solver == 'interFoam':
        alpha_values = [0.0 for item in range(nCells)]
    if solver == 'driftFluxSimphonyFoam':
        vdj_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
        sigma_mu_values = [[0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0] for item in range(nCells)]
    for cell in mesh.iter_cells():
        p_values[mesh._uuidToFoamLabel[cell.uid]] = \
            cell.data[ID_pressure]
        U_values[mesh._uuidToFoamLabel[cell.uid]] = \
            list(cell.data[CUBA.VELOCITY])
        if solver == 'driftFluxSimphonyFoam' or solver == 'interFoam':
            alpha_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.VOLUME_FRACTION]
        if solver == 'driftFluxSimphonyFoam':
            vdj_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.RELATIVE_VELOCITY]
            sigma_mu_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.HOMOGENIZED_STRESS_TENSOR]
    foamface.setAllCellData(mesh.name, name_pressure, 0, p_values,
                            dataDimensionMap[dataKeyMap[name_pressure]])
    foamface.setAllCellVectorData(mesh.name, "U", 0, U_values,
                                  dataDimensionMap[dataKeyMap["U"]])
    if solver == 'driftFluxSimphonyFoam' or solver == 'interFoam':
        foamface.setAllCellData(mesh.name, "alpha.phase1", 0, alpha_values,
                                dataDimensionMap[dataKeyMap["alpha.phase1"]])
    if solver == 'driftFluxSimphonyFoam':
        foamface.setAllCellTensorData(mesh.name, "Sigma", 0, sigma_mu_values,
                                      dataDimensionMap[dataKeyMap["Sigma"]])

    # Refresh boundary conditions

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
        elif isinstance(patch, tuple) and patch[0] == "flowRate":
            myDict = myDict + "\t type \t flowRateInletVelocity;\n"
            myDict = myDict + "\t volumetricFlowRate \t" \
                + str(patch[1]) + ";\n"
        elif isinstance(patch, tuple) and patch[0] ==\
                "shearStressPowerLawSlipVelocity":
            myDict = myDict + "\t type \t shearStressPowerLawSlipVelocity;\n"
            parameters = patch[1]
            for key in parameters:
                myDict = myDict + "\t " + key + "\t " + str(parameters[key])\
                    + ";\n"
            myDict = myDict + "\t value \t uniform (0 0 0);\n"

        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, "U", myDict)

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

    if solver == 'driftFluxSimphonyFoam' or solver == 'interFoam':
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
                elif volumeFractionBCs[boundary][0] == "wettingAngle":
                    myDict = myDict + "\t type \t constantAlphaContactAngle;\n"
                    myDict = myDict + "\t theta0 \t" +\
                        str(volumeFractionBCs[boundary][1]) + ";\n"
                    myDict = myDict + "\t limit \t gradient;\n"
                    myDict = myDict + "\t value \t uniform 0;\n"
            myDict = myDict + "}\n"

        foamface.setBC(mesh.name, "alpha.phase1", myDict)


def get_foam_solver(CM):
    """ gives the name of the solver

        Parameters
        ----------
        CM : DataContainer
            Computational Methods

    """

    GE = CM[CUBAExt.GE]
    solver = "pimpleFoam"
    if CUBAExt.LAMINAR_MODEL in GE:
        if CUBAExt.VOF_MODEL in GE:
            solver = "interFoam"
        elif CUBAExt.MIXTURE_MODEL in GE:
            solver = "driftFluxSimphonyFoam"
        else:
            solver = "pimpleFoam"
        return solver
    else:
        error_str = "GE does not define supported solver: GE = {}"
        raise NotImplementedError(error_str.format(GE))


def write_dictionary(path, dictionary_name, dictionary_content):
    """write dictionary to file

        Parameters
        ----------
        path : str
            path to the file directory
        dictionary_name : str
            name of the dictionary
        dictionary_content : str
            content of the dictionary

    """

    heading = head.format(version='2.4', foamclass='dictionary',
                          location='system', foamobject=dictionary_name)

    foamface.writePathDictionary(path, dictionary_name, heading,
                                 dictionary_content)


def create_directories(case_directory):
    """Create default directories

    Parameters
    ----------
    case_directory : str
        directory to create directories on.

    """

    directories = ('constant', 'system', '0',
                   os.path.join('constant', 'polyMesh'))
    for dir in directories:
        directory = os.path.join(case_directory, dir)
        if not os.path.exists(directory):
            os.makedirs(directory)


def add_user_libs(map_content, io_library):
    """ Adds user libs to controlDict

    Parameters
    ----------
    map_content : dictionary
        specified solver file map
    io_library : bool
        use IO user libraries
    """

    controlDict = map_content['controlDict']
    if io_library:
        controlDict['libs'] = userLibsIO
    else:
        controlDict['libs'] = userLibs


def get_dictionary_maps(solver, io_library=True):
    """ returns dictionary map for specified solver

    Parameters
    ----------
    solver : str
        solver name
    io_library : bool
        use IO user libraries

    """
    mapContent = dictionaryMaps[solver]
    add_user_libs(mapContent, io_library)

    return mapContent
