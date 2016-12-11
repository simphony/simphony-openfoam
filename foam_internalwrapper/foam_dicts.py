""" foam_dicts

Management of OpenFOAM dicts and Fields

"""
import os

import simphonyfoaminterface as foamface
from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

from foam_controlwrapper.foam_variables import (dataDimensionMap, dataKeyMap,
                                                dataNameMap)

userLibs = """
("libshearStressPowerLawSlipVelocity.so"  "libtwoPhaseProperties.so")
"""
userLibsIO = """
("libshearStressPowerLawSlipVelocityIO.so" "libtwoPhaseProperties.so")
"""

dictionaryMaps = \
    {'simpleFoam':
        {'transportProperties':
            {'transportModel': 'Newtonian',
             'nu nu     [ 0 2 -1 0 0 0 0 ]': '0.0',
             'rho rho [ 1 -3 0 0 0 0 0 ]': '0.0',
             'stressModel': 'standard',
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
             {'application': 'simpleFoam',
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
              'runTimeModifiable': 'yes'
              }
         },
     'pimpleFoam':
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
              'runTimeModifiable': 'yes'
              }
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
     'driftFluxFoam':
         {'controlDict':
              {'application': 'driftFluxFoam',
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
                   {'dummy': '0'}
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

dictionaryTemplates = {'simpleFoam':
                       {'transportProperties':
                           """
transportModel Newtonian;

nu nu     [ 0 2 -1 0 0 0 0 ]     0.0 ;
rho  rho [ 1 -3 0 0 0 0 0 ]': '0.0';
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
        relTol           0.1;

        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     2;

        cacheAgglomeration true;
        nCellsInCoarsestLevel 10;
        agglomerator     faceAreaPair;
        mergeLevels      1;
    }


    "(U)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

}

SIMPLE
{
    nNonOrthogonalCorrectors 0;

    residualControl
    {
        p               1e-2;
        U               1e-3;
        "(k|epsilon|omega)" 1e-3;
    }
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
        k               0.7;
        epsilon         0.7;
        R               0.7;
        nuTilda         0.7;
    }
}
                        """,
                            'fvSchemes':
                                """
ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss upwind;
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,R)      bounded Gauss upwind;
    div(R)          Gauss linear;
    div(SigmaPrime)       Gauss linear;
    div(phi,nuTilda) bounded Gauss upwind;
    div((nu*dev(T(grad(U))))) Gauss linear;
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
    p               ;
}
                        """},
                       'pimpleFoam':
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
    div(SigmaPrime)       Gauss linear;
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
                       'driftFluxFoam':
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


def get_rheology_model_coeffs(rheology_model, cuds):
    material_uid = rheology_model.material
    material = cuds.get_by_uid(material_uid)
    density = material._data[CUBA.DENSITY]
    if isinstance(rheology_model, api.HerschelBulkleyModel):
        return {'name': 'HerschelBulkley',
                'nu0': str(rheology_model.initial_viscosity/density),
                'tau0': str(rheology_model.relaxation_time/density),
                'k': str(rheology_model.linear_constant/density),
                'n': str(rheology_model.power_law_index)}
    elif isinstance(rheology_model, api.BirdCarreauModel):
        return {'name': 'BirdCarreau',
                'nu0': str(rheology_model.initial_viscosity/density),
                'nuInf': str(rheology_model.maximum_viscosity/density),
                'k': str(rheology_model.linear_constant/density),
                'n': str(rheology_model.power_law_index)}
    elif isinstance(rheology_model, api.CrossPowerLawModel):
        return {'name': 'CrossPowerLaw',
                'nu0': str(rheology_model.initial_viscosity/density),
                'nuInf': str(rheology_model.maximum_viscosity/density),
                'm': str(rheology_model.linear_constant/density),
                'n': str(rheology_model.power_law_index)}
    elif isinstance(rheology_model, api.BinghamPlasticModel):
        # this is used in mixture model so no division by density
        return {'name': 'BinghamPlastic',
                'coeff': str(rheology_model.linear_constant[0]),
                'exponent': str(rheology_model.power_law_index[0]),
                'BinghamCoeff': str(rheology_model.linear_constant[1]),
                'BinghamOffset': str(rheology_model.power_law_index[1]),
                'muMax': str(rheology_model.maximum_viscosity)}

    error_str = "Rheology model {} not supported"
    raise ValueError(error_str.format(rheology_model.__class__.__name__))


def modifyNumerics(mesh, cuds, solver='pimpleFoam', io=False):
    """ Modifies the numerical parameters of the simulation

        Parameters
        ----------
        mesh : FoamMesh
            mesh to be added.
        cuds : CUDS
            CUDS
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

    for int_time in cuds.iter(api.IntegrationTime):
        current = int_time.current
        final = int_time.final
        deltaT = int_time.size
        break
    interval = final-current
    endTime = interval + mesh._time
    for solv_param in cuds.iter(api.SolverParameter):
        if CUBA.MAXIMUM_COURANT_NUMBER in solv_param._data:
            mapContent['controlDict']['maxCo'] =\
                solv_param._data[CUBA.MAXIMUM_COURANT_NUMBER]
            mapContent['controlDict']['maxAlphaCo'] =\
                solv_param._data[CUBA.MAXIMUM_COURANT_NUMBER]

    mapContent['controlDict']['startTime'] = str(mesh._time)
    mapContent['controlDict']['deltaT'] = str(deltaT)
    mapContent['controlDict']['endTime'] = str(endTime)
    mapContent['controlDict']['maxDeltaT'] = str(deltaT)
    if io:
        mapContent['controlDict']['writeControl'] = 'adjustableRunTime'
        mapContent['controlDict']['writeInterval'] = str(interval)

    controlDict = parse_map(mapContent['controlDict'])

    cfd = get_first(cuds.iter(api.Cfd))
    if is_empty(cuds.iter(api.Material)):
        error_str = "Materials missing on cuds"
        raise ValueError(error_str)

    materials = list(cuds.iter(api.Material))
    if solver == 'pimpleFoam' or solver == 'simpleFoam':
        control = mapContent['transportProperties']
        rheology_model = cfd.rheology_model
        material = materials[0]
        control['rho rho [ 1 -3 0 0 0 0 0 ]'] = material._data[CUBA.DENSITY]
        if rheology_model is not None and not \
                isinstance(rheology_model, api.NewtonianFluidModel):
            rheology_model_coeffs = get_rheology_model_coeffs(rheology_model,
                                                              cuds)
            rheology_model_name = rheology_model_coeffs['name']
            control['transportModel'] = rheology_model_name
            vmc = rheology_model_name + 'Coeffs'
            for coeff in rheology_model_coeffs:
                if coeff != 'name':
                    control[vmc][coeff] =\
                        control[vmc][coeff].format(
                        rheology_model_coeffs[coeff])
            control['nu nu     [ 0 2 -1 0 0 0 0 ]'] = \
                rheology_model.initial_viscosity
        else:
            control['nu nu     [ 0 2 -1 0 0 0 0 ]'] = \
                material._data[CUBA.DYNAMIC_VISCOSITY] / \
                material._data[CUBA.DENSITY]

        if not_empty(cuds.iter(api.StressModel)):
            stress_model = get_first(cuds.iter(api.StressModel))
            if isinstance(stress_model, api.StandardStressModel):
                control["stressModel"] = 'standard'
            elif isinstance(stress_model, api.MesoScopicStressModel):
                control["stressModel"] = 'fromMesoscale'
        else:
            control["stressModel"] = 'standard'

        transportPropertiesDict = parse_map(control)
    elif solver == 'interFoam':
        control = mapContent['transportProperties']
        for i in range(2):
            material = materials[i]
            density = material._data[CUBA.DENSITY]
            viscosity = material._data[CUBA.DYNAMIC_VISCOSITY]
            rheology_model = None
            if type(cfd.rheology_model) is list:
                for rm in cfd.rheology_model:
                    if rm._data[CUBA.MATERIAL] == material.uid:
                        rheology_model = rm
                        break
            elif CUBA.MATERIAL in cfd.rheology_model._data:
                if cfd.rheology_model._data[CUBA.MATERIAL] == material.uid:
                    rheology_model = cfd.rheology_model
            if rheology_model is not None:
                rheology_model_name = rheology_model.__class__.__name__
            else:
                rheology_model_name = 'Newtonian'
            phase_name = 'phase' + str(i + 1)
            if rheology_model_name == 'Newtonian':
                control[phase_name]['nu              nu [ 0 2 -1 0 0 0 0 ]']\
                    = material._data[CUBA.DYNAMIC_VISCOSITY] / \
                    material._data[CUBA.DENSITY]
                control[phase_name]['rho             rho [ 1 -3 0 0 0 0 0 ]']\
                    = material._data[CUBA.DENSITY]
            else:
                rheology_model_coeffs = get_rheology_model_coeffs(
                    rheology_model, cuds)
                vmc = rheology_model_coeffs['name']
                control[phase_name]['transportModel'] = vmc
                for coeff in rheology_model_coeffs:
                    if coeff != 'name':
                        control[phase_name][vmc][coeff] =\
                            control[phase_name][vmc][coeff].\
                            format(rheology_model_coeffs[coeff])

        if not_empty(cuds.iter(api.SurfaceTensionRelation)):
            surface_tension = get_first(cuds.iter(api.SurfaceTensionRelation))
            surface_tension = surface_tension.surface_tension
            control['sigma           sigma [ 1 0 -2 0 0 0 0 ]'] = \
                surface_tension
        else:
            error_str = "Surface tension not specified"
            raise ValueError(error_str)
        transportPropertiesDict = parse_map(control)
    elif solver == 'driftFluxFoam':
        control = mapContent['transportProperties']
        for i in range(2):
            mixture_model = get_first(cuds.iter(api.MixtureModel))
            disperse_material_uid = mixture_model.disperse
            if materials[i].uid == disperse_material_uid:
                phase_name = 'phase1'
            else:
                phase_name = 'phase2'

            density = materials[i]._data[CUBA.DENSITY]
            viscosity = materials[i]._data[CUBA.DYNAMIC_VISCOSITY]
            rheology_models = cfd.rheology_model
            rheology_model = None
            if rheology_models is not None:
                if type(rheology_models) is list:
                    for rm in rheology_models:
                        if rm._data[CUBA.MATERIAL] == materials[i].uid:
                            rheology_model = rm
                            break
                    if rheology_model is None:
                        error_str = "Rheology model not specified"
                        error_str += "for material {}"
                        raise ValueError(error_str.format(materials[i].name))

                elif CUBA.MATERIAL in rheology_models._data:
                    if rheology_models._data[CUBA.MATERIAL] == material.uid:
                        rheology_model = rheology_models
            if rheology_model is not None:
                rheology_model_coeffs = \
                    get_rheology_model_coeffs(rheology_model, cuds)
                rheology_model_name = rheology_model_coeffs['name']
                control[phase_name]["transportModel"] = rheology_model_name
                vmc = rheology_model_name + 'Coeffs'
                for coeff in rheology_model_coeffs:
                    if coeff != 'name':
                        control[phase_name][vmc][coeff] = \
                            rheology_model_coeffs[coeff]
                        control[phase_name]['rho'] = \
                            materials[0]._data[CUBA.DENSITY]
            else:
                if materials[i].uid == disperse_material_uid:
                    rheology_model_name = 'dummyViscosity'
                else:
                    rheology_model_name = 'Newtonian'
                control[phase_name]["transportModel"] = rheology_model_name
                control[phase_name]['nu'] = viscosity/density
                control[phase_name]['rho'] = density

        if not_empty(cuds.iter(api.StressModel)):
            stress_model = get_first(cuds.iter(api.StressModel))
            if isinstance(stress_model, api.StandardStressModel):
                control["stressModel"] = 'standard'
            elif isinstance(stress_model, api.MesoScopicStressModel):
                control["stressModel"] = 'fromMesoscale'
        else:
            error_str = "Stress model not specified"
            raise ValueError(error_str)

        if not_empty(cuds.iter(api.RelativeVelocityModel)):
            rel_vel_model = get_first(cuds.iter(api.RelativeVelocityModel))
            if isinstance(rel_vel_model, api.SimpleRelativeVelocityModel):
                rel_vel_model_name = 'simple'
            elif isinstance(rel_vel_model,
                            api.MesoScopicRelativeVelocityModel):
                rel_vel_model_name = 'fromMesoscale'
            else:
                error_str = "Relative velocity model not specified"
                raise ValueError(error_str)
            control["relativeVelocityModel"] = rel_vel_model_name

        else:
            error_str = "Relative velocity model not specified"
            raise ValueError(error_str)
        if isinstance(rel_vel_model, api.SimpleRelativeVelocityModel):
            relVelModelCoeffs = rel_vel_model_name + 'Coeffs'
            V0 = rel_vel_model.diffusion_velocity
            control[relVelModelCoeffs]["V0"] = "( " + str(V0[0]) + " " + \
                str(V0[1]) + " " + str(V0[2]) + " )"
            control[relVelModelCoeffs]["a"] = rel_vel_model.linear_constant

        transportPropertiesDict = parse_map(control)
    foamface.modifyNumerics(mesh.name, fvSchemesDict, fvSolutionDict,
                            controlDict, transportPropertiesDict, int(io))
    if io:
        # write numerics to case directory
        foamface.writeNumerics(mesh.name)

    if solver == 'driftFluxFoam' or solver == 'interFoam':
        if not_empty(cuds.iter(api.GravityModel)):
            gv = get_first(cuds.iter(api.GravityModel))
            if io:
                control = mapContent['g']
                if gv is not None:
                    g = gv.acceleration
                    control['value'] =\
                        "( " \
                        + str(g[0]) + " " + str(g[1]) + " " + str(g[2]) + " )"
                    gravityDict = parse_map(control)
                    foamface.modifyDictionary(mesh.name, 'g', gravityDict)
                    foamface.writeDictionary(mesh.name, 'g', True)
            else:
                if gv is not None:
                    g = gv.acceleration
                    foamface.modifyUniformVectorField(mesh.name, 'g', list(g))


def check_boundary_names(bc_names, boundary_names):
    if not set(bc_names).issubset(set(boundary_names)):
        error_str = "Boundary name(s) used in boundary conditions "
        error_str += "does not exist.\nUsed name(s) are: {}\n"
        error_str += "Boundary names defined in the mesh are: {}\n"
        raise ValueError(
            error_str.format(
                list(set(bc_names).difference(set(boundary_names))),
                boundary_names))


def get_foam_boundary_condition(condition):
    """ Return corresponding foam condition from Condition type
    """
    dyn_pres = condition._data[CUBA.VARIABLE] == CUBA.DYNAMIC_PRESSURE
    if isinstance(condition, api.Dirichlet):
        patch = []
        patch.append('fixedValue')
        patch.append(condition._data[condition._data[CUBA.VARIABLE]])
        return patch
    elif isinstance(condition, api.Neumann):
        # for dynamic pressure use fixedFluxPressure
        if dyn_pres:
            return 'fixedFluxPressure'
        else:
            return 'zeroGradient'
    elif isinstance(condition, api.Empty):
        return 'empty'
    elif isinstance(condition, api.ShearStressPowerLawSlipVelocity):
        patch = []
        patch.append('shearStressPowerLawSlipVelocity')
        patch.append({'rho': str(condition.density),
                      'beta': str(condition.linear_constant),
                      'n': str(condition.power_law_index)})
        return patch
    elif isinstance(condition, api.InletOutlet):
        if condition._data[CUBA.VARIABLE] == CUBA.VELOCITY:
            patch = []
            patch.append('pressureInletOutletVelocity')
            patch.append(condition._data[condition._data[CUBA.VARIABLE]])
            return patch
        elif condition._data[CUBA.VARIABLE] == CUBA.VOLUME_FRACTION:
            patch = []
            patch.append('inletOutlet')
            patch.append(condition._data[condition._data[CUBA.VARIABLE]])
            return patch
        else:
            error_str = "Boundary condition not supported:\n"
            error_str += " condition: {}\n"
            error_str += " variable: {}"
            raise ValueError(error_str.format(condition.__class__.__name__,
                                              condition._data[CUBA.VARIABLE]))

    elif isinstance(condition, api.SlipVelocity):
        return 'slip'
    elif isinstance(condition, api.WettingAngle):
        patch = []
        patch.append('wettingAngle')
        patch.append(condition.contact_angle)
        return patch
    else:
        patch = None
        return patch


def modifyFields(mesh, cuds, solver='pimpleFoam'):
    """ Modifies the internal fields and boundary conditions

        Parameters
        ----------
        mesh : FoamMesh
            mesh to be added.
        cuds : CUDS
            CUDS
        solver : str
            solver name

    """

    bc_names = []
    for boundary in cuds.iter(api.Boundary):
        bc_names.append(boundary.name)
    check_boundary_names(bc_names, mesh._boundaries.keys())

    if solver == 'pimpleFoam' or solver == 'simpleFoam':
        name_pressure = dataNameMap[CUBA.PRESSURE]
        ID_pressure = CUBA.PRESSURE
    elif solver == 'driftFluxFoam' or solver == 'interFoam':
        name_pressure = dataNameMap[CUBA.DYNAMIC_PRESSURE]
        ID_pressure = CUBA.DYNAMIC_PRESSURE

    # Refresh boundary conditions

    myDict = ""

    for boundary in cuds.iter(api.Boundary):
        patch = None
        for condition in boundary.condition:
            if condition.variable == CUBA.VELOCITY:
                patch = get_foam_boundary_condition(condition)
                break
        if patch is None:
            error_str = "Boundary condition not specified:\n"
            error_str += " boundary: {}\n"
            error_str += " variable: {}"
            raise ValueError(error_str.format(boundary.name,
                                              CUBA.VELOCITY))

        myDict = myDict + str(boundary.name) + "\n{\n"
        if patch == "zeroGradient":
            myDict = myDict + "\t type \t zeroGradient;\n"
            myDict = myDict + "\t value \t uniform (0 0 0);\n"
        elif patch == "empty":
            myDict = myDict + "\t type \t empty;\n"
        elif patch == "slip":
            myDict = myDict + "\t type \t slip;\n"
        elif isinstance(patch, list) and patch[0] == \
                "pressureInletOutletVelocity":
            myDict = myDict + "\t type \t pressureInletOutletVelocity;\n"
            myDict = myDict + "\t value \t uniform (" \
                + str(patch[1][0]) + " " \
                + str(patch[1][1]) + " " \
                + str(patch[1][2]) + ");\n"
        elif isinstance(patch, list) and patch[0] == "fixedValue":
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform (" \
                + str(patch[1][0]) + " " \
                + str(patch[1][1]) + " " \
                + str(patch[1][2]) + ");\n"
        elif isinstance(patch, list) and patch[0] == "flowRate":
            myDict = myDict + "\t type \t flowRateInletVelocity;\n"
            myDict = myDict + "\t volumetricFlowRate \t" \
                + str(patch[1]) + ";\n"
        elif isinstance(patch, list) and patch[0] ==\
                "shearStressPowerLawSlipVelocity":
            myDict = myDict + "\t type \t shearStressPowerLawSlipVelocity;\n"
            parameters = patch[1]
            for key in parameters:
                myDict = myDict + "\t " + key + "\t " + parameters[key]\
                    + ";\n"
            myDict = myDict + "\t value \t uniform (0 0 0);\n"

        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, dataNameMap[CUBA.VELOCITY], myDict)

    myDict = ""

    for boundary in cuds.iter(api.Boundary):
        patch = None
        for condition in boundary.condition:
            if condition.variable == ID_pressure:
                patch = get_foam_boundary_condition(condition)
                break
        if patch is None:
                error_str = "Boundary condition not specified:\n"
                error_str += " boundary: {}\n"
                error_str += " variable: {}"
                raise ValueError(error_str.format(boundary.name,
                                                  ID_pressure))
        myDict = myDict + str(boundary.name) + "\n{\n"
        if patch == "zeroGradient":
            myDict = myDict + "\t type \t zeroGradient;\n"
            myDict = myDict + "\t value \t uniform 0;\n"
        elif patch == "fixedFluxPressure":
            myDict = myDict + "\t type \t fixedFluxPressure;\n"
            myDict = myDict + "\t value \t uniform 0;"
        elif patch == "empty":
            myDict = myDict + "\t type \t empty;\n"
        elif isinstance(patch, list) and\
                patch[0] == "fixedValue":
            myDict = myDict + "\t type \t fixedValue;\n"
            myDict = myDict + "\t value \t uniform " \
                            + str(patch[1]) + ";\n"
        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, name_pressure, myDict)

    if solver == 'driftFluxFoam' or solver == 'interFoam':
        myDict = ""

        for boundary in cuds.iter(api.Boundary):
            patch = None
            for condition in boundary.condition:
                if condition.variable == CUBA.VOLUME_FRACTION:
                    patch = get_foam_boundary_condition(condition)
                    break
            if patch is None:
                error_str = "Boundary condition not specified:\n"
                error_str += " boundary: {}\n"
                error_str += " variable: {}"
                raise ValueError(error_str.format(boundary.name,
                                                  CUBA.VOLUME_FRACTION))

            myDict = myDict + str(boundary.name) + "\n{\n"
            if patch == "zeroGradient":
                myDict = myDict + "\t type \t zeroGradient;\n"
                myDict = myDict + "\t value \t uniform 0;\n"
            elif patch == "empty":
                myDict = myDict + "\t type \t empty;\n"
            elif isinstance(patch, list):
                if patch[0] == "inletOutlet":
                    myDict = myDict + "\t type \t inletOutlet;\n"
                    myDict = myDict + "\t inletValue \t uniform " \
                        + str(patch[1]) + ";\n"
                elif patch[0] == "fixedValue":
                    myDict = myDict + "\t type \t fixedValue;\n"
                    myDict = myDict + "\t value \t uniform " \
                        + str(patch[1]) + ";\n"
                elif patch[0] == "wettingAngle":
                    myDict = myDict + "\t type \t constantAlphaContactAngle;\n"
                    myDict = myDict + "\t theta0 \t" +\
                        str(patch[1]) + ";\n"
                    myDict = myDict + "\t limit \t gradient;\n"
                    myDict = myDict + "\t value \t uniform 0;\n"
            myDict = myDict + "}\n"
        foamface.setBC(mesh.name, dataNameMap[CUBA.VOLUME_FRACTION], myDict)

    nCells = foamface.getCellCount(mesh.name)
    p_values = [0.0 for item in range(nCells)]
    U_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
    if solver == 'driftFluxFoam' or solver == 'interFoam':
        alpha_values = [0.0 for item in range(nCells)]
    if solver == 'driftFluxFoam' or solver == 'simpleFoam':
        vdj_values = [[0.0, 0.0, 0.0] for item in range(nCells)]
        sigma_mu_values = [[0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0] for item in range(nCells)]
    for cell in mesh._iter_cells():
        p_values[mesh._uuidToFoamLabel[cell.uid]] = \
            cell.data[ID_pressure]
        U_values[mesh._uuidToFoamLabel[cell.uid]] = \
            list(cell.data[CUBA.VELOCITY])
        if solver == 'driftFluxFoam' or solver == 'interFoam':
            alpha_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.VOLUME_FRACTION]
        if solver == 'driftFluxFoam':
            vdj_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.RELATIVE_VELOCITY]
        if solver == 'driftFluxFoam' or solver == 'simpleFoam':
            sigma_mu_values[mesh._uuidToFoamLabel[cell.uid]] = \
                cell.data[CUBA.HOMOGENIZED_STRESS_TENSOR]
    foamface.setAllCellData(mesh.name, name_pressure, 0, p_values,
                            dataDimensionMap[dataKeyMap[name_pressure]])
    foamface.setAllCellVectorData(mesh.name, dataNameMap[CUBA.VELOCITY],
                                  0, U_values,
                                  dataDimensionMap[dataKeyMap[
                                      dataNameMap[CUBA.VELOCITY]]])
    if solver == 'driftFluxFoam' or solver == 'interFoam':
        foamface.setAllCellData(mesh.name, dataNameMap[CUBA.VOLUME_FRACTION],
                                0, alpha_values,
                                dataDimensionMap[dataKeyMap[
                                    dataNameMap[CUBA.VOLUME_FRACTION]]])
    if solver == 'driftFluxFoam' or solver == 'simpleFoam':
        foamface.setAllCellTensorData(mesh.name, dataNameMap[
                                      CUBA.HOMOGENIZED_STRESS_TENSOR],
                                      0, sigma_mu_values,
                                      dataDimensionMap[dataKeyMap[
                                          dataNameMap[
                                              CUBA.HOMOGENIZED_STRESS_TENSOR
                                              ]]])


def is_empty(generator):
    if generator is None:
        return True
    try:
        first = generator.next()
        return False
    except StopIteration:
        return True


def not_empty(generator):
    if generator is None:
        return False
    try:
        first = generator.next()
        return True
    except StopIteration:
        return False


def get_first(generator):
    if generator is None:
        return None
    try:
        first = generator.next()
        return first
    except StopIteration:
        return None


def get_simphony_io_solver(foam_solver):
    simphony_solvers = {'pimpleFoam': 'pimpleSimphonyFoam',
                        'simpleFoam': 'simpleSimphonyFoam',
                        'interFoam': 'interFoam',
                        'driftFluxFoam': 'driftFluxSimphonyFoam'}
    return simphony_solvers[foam_solver]


def get_foam_solver(cuds):
    """ gives the name of the solver

        Parameters
        ----------
        cuds : CUDS
            CUDS container

    """

    if not_empty(cuds.iter(api.SolverParameter)):
        for sp in cuds.iter(api.SolverParameter):
            if CUBA.STEADY_STATE in sp._data:
                solver = 'simpleFoam'
                return solver
    if is_empty(cuds.iter(api.Cfd)):
        error_str = "CFD physics model not present in the cuds"
        raise ValueError(error_str)
    cfd = get_first(cuds.iter(api.Cfd))
    laminar = isinstance(cfd.turbulence_model, api.LaminarFlowModel)
    mixture_model = not_empty(cuds.iter(api.MixtureModel))

    vof = not_empty(cuds.iter(api.FreeSurfaceModel))

    solver = "pimpleFoam"
    if laminar:
        if vof:
            solver = "interFoam"
        elif mixture_model:
            solver = "driftFluxFoam"
        return solver
    else:
        error_str = "Solver for the model not found:\n"
        error_str += " Laminar: {}\n"
        error_str += " Freesurface: {}\n"
        error_str += " Mixture model: {}"
        raise NotImplementedError(error_str.format(laminar,
                                                   vof, mixture_model))


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
