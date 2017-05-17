""" foam_dicts

Management of OpenFOAM dicts and Fields

"""
import os

import simphonyfoaminterface as foamface
from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

from foam_controlwrapper.foam_variables import (dataNameMap, phaseNames)

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
    material = cuds.get(material_uid)
    density = material.data[CUBA.DENSITY]
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

    for int_time in cuds.iter(item_type=CUBA.INTEGRATION_TIME):
        current = int_time.current
        final = int_time.final
        deltaT = int_time.size
        break

    interval = final-current
    endTime = interval + mesh._time
    for solv_param in cuds.iter(item_type=CUBA.SOLVER_PARAMETER):
        if CUBA.MAXIMUM_COURANT_NUMBER in solv_param.data:
            mapContent['controlDict']['maxCo'] =\
                solv_param.data[CUBA.MAXIMUM_COURANT_NUMBER]
            mapContent['controlDict']['maxAlphaCo'] =\
                solv_param.data[CUBA.MAXIMUM_COURANT_NUMBER]
        if CUBA.NUMBER_OF_PHYSICS_STATES in solv_param.data:
                interval = interval /\
                    solv_param.data[CUBA.NUMBER_OF_PHYSICS_STATES]

    mapContent['controlDict']['startTime'] = str(mesh._time)
    mapContent['controlDict']['deltaT'] = str(deltaT)
    mapContent['controlDict']['endTime'] = str(endTime)
    mapContent['controlDict']['maxDeltaT'] = str(deltaT)
    if io:
        mapContent['controlDict']['writeControl'] = 'adjustableRunTime'
        mapContent['controlDict']['writeInterval'] = str(interval)

    controlDict = parse_map(mapContent['controlDict'])

    cfd = get_first(cuds.iter(item_type=CUBA.CFD))
    if is_empty(cuds.iter(item_type=CUBA.MATERIAL)):
        error_str = "Materials missing on cuds"
        raise ValueError(error_str)

    materials = list(cuds.iter(item_type=CUBA.MATERIAL))
    if solver == 'pimpleFoam' or solver == 'simpleFoam':
        control = mapContent['transportProperties']
        rheology_model = cfd.rheology_model
        material = materials[0]
        control['rho rho [ 1 -3 0 0 0 0 0 ]'] = material.data[CUBA.DENSITY]
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
                material.data[CUBA.DYNAMIC_VISCOSITY] / \
                material.data[CUBA.DENSITY]

        if not_empty(cuds.iter(item_type=CUBA.STRESS_MODEL)):
            stress_model = get_first(cuds.iter(item_type=CUBA.STRESS_MODEL))
            if isinstance(stress_model, api.StandardStressModel):
                control["stressModel"] = 'standard'
            elif isinstance(stress_model, api.MesoscopicStressModel):
                control["stressModel"] = 'fromMesoscale'
        else:
            control["stressModel"] = 'standard'

        transportPropertiesDict = parse_map(control)
    elif solver == 'interFoam':

        control = mapContent['transportProperties']

        for phase_name, material in mesh._foamPhaseNameToMaterial.iteritems():
            density = material.data[CUBA.DENSITY]
            viscosity = material.data[CUBA.DYNAMIC_VISCOSITY]
            rheology_model = None
            if type(cfd.rheology_model) is list:
                for rm in cfd.rheology_model:
                    if rm.data[CUBA.MATERIAL] == material.uid:
                        rheology_model = rm
                        break
            elif CUBA.MATERIAL in cfd.rheology_model.data:
                if cfd.rheology_model.data[CUBA.MATERIAL] == material.uid:
                    rheology_model = cfd.rheology_model
            if rheology_model is not None:
                rheology_model_name = rheology_model.__class__.__name__
            else:
                rheology_model_name = 'Newtonian'
            if rheology_model_name == 'Newtonian':
                control[phase_name]['nu              nu [ 0 2 -1 0 0 0 0 ]']\
                    = material.data[CUBA.DYNAMIC_VISCOSITY] / \
                    material.data[CUBA.DENSITY]
                control[phase_name]['rho             rho [ 1 -3 0 0 0 0 0 ]']\
                    = material.data[CUBA.DENSITY]
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

        if not_empty(cuds.iter(item_type=CUBA.SURFACE_TENSION_RELATION)):
            surface_tension = get_first(cuds.iter(
                item_type=CUBA.SURFACE_TENSION_RELATION))
            surface_tension = surface_tension.surface_tension
            control['sigma           sigma [ 1 0 -2 0 0 0 0 ]'] = \
                surface_tension
        else:
            error_str = "Surface tension not specified"
            raise ValueError(error_str)
        transportPropertiesDict = parse_map(control)
    elif solver == 'driftFluxFoam':
        control = mapContent['transportProperties']
        mixture_model = get_first(cuds.iter(item_type=CUBA.MIXTURE_MODEL))
        disperse_material_uid = mixture_model.disperse
        for pname, material in mesh._foamPhaseNameToMaterial.iteritems():
            # always disperse phase name phaseNames[0]
            if disperse_material_uid == material.uid:
                phase_name = phaseNames[0]
            else:
                phase_name = phaseNames[1]
            density = material.data[CUBA.DENSITY]
            viscosity = material.data[CUBA.DYNAMIC_VISCOSITY]
            rheology_models = cfd.rheology_model
            rheology_model = None
            if rheology_models is not None:
                if type(rheology_models) is list:
                    for rm in rheology_models:
                        if rm.data[CUBA.MATERIAL] == material.uid:
                            rheology_model = rm
                            break
                    if rheology_model is None:
                        error_str = "Rheology model not specified"
                        error_str += "for material {}"
                        raise ValueError(error_str.format(material.name))

                elif CUBA.MATERIAL in rheology_models.data:
                    if rheology_models.data[CUBA.MATERIAL] == material.uid:
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
                            materials[0].data[CUBA.DENSITY]
            else:
                if material.uid == disperse_material_uid:
                    rheology_model_name = 'dummyViscosity'
                else:
                    rheology_model_name = 'Newtonian'
                control[phase_name]["transportModel"] = rheology_model_name
                control[phase_name]['nu'] = viscosity/density
                control[phase_name]['rho'] = density

        if not_empty(cuds.iter(item_type=CUBA.STRESS_MODEL)):
            stress_model = get_first(cuds.iter(item_type=CUBA.STRESS_MODEL))
            if isinstance(stress_model, api.StandardStressModel):
                control["stressModel"] = 'standard'
            elif isinstance(stress_model, api.MesoscopicStressModel):
                control["stressModel"] = 'fromMesoscale'
        else:
            error_str = "Stress model not specified"
            raise ValueError(error_str)

        if not_empty(cuds.iter(item_type=CUBA.RELATIVE_VELOCITY_MODEL)):
            rel_vel_model = get_first(cuds.iter(
                item_type=CUBA.RELATIVE_VELOCITY_MODEL))
            if isinstance(rel_vel_model, api.SimpleRelativeVelocityModel):
                rel_vel_model_name = 'simple'
            elif isinstance(rel_vel_model,
                            api.MesoscopicRelativeVelocityModel):
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
        if not_empty(cuds.iter(item_type=CUBA.GRAVITY_MODEL)):
            gv = get_first(cuds.iter(item_type=CUBA.GRAVITY_MODEL))
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


def get_condition_variable(condition, solver):
    """Returns condition variable as CUBA key.
    If not defined raise an error
    """

    if hasattr(condition, 'variables'):
        variable = CUBA[condition.variables[0].split('.')[1]]
    else:
        if CUBA.VARIABLE not in condition.data:
            error_str = "condition.data[CUBA.VARIABLE] must be defined\n"
            error_str += " condition: {}\n"
            error_str += " name: {}"
            raise ValueError(error_str.format(condition.__class__.__name__,
                                              condition.name))
        else:
            variable = condition.data[CUBA.VARIABLE]

    if (solver == 'driftFluxFoam' or solver == 'interFoam') and \
            variable == CUBA.PRESSURE:
        condition.data[CUBA.DYNAMIC_PRESSURE] = condition.data[CUBA.PRESSURE]
        return CUBA.DYNAMIC_PRESSURE
    else:
        return variable


def get_foam_boundary_condition(condition, phaseNameToMaterial, solver):
    """ Return corresponding foam condition from Condition type
    """
    if isinstance(condition, api.ShearStressPowerLawSlipVelocity):
        patch = []
        patch.append('shearStressPowerLawSlipVelocity')
        patch.append({'rho': str(condition.density),
                      'beta': str(condition.linear_constant),
                      'n': str(condition.power_law_index)})
        return patch
    elif isinstance(condition, api.MixedCondition):
        variable = get_condition_variable(condition, solver)
        if isinstance(condition, api.InletOutletVelocity):
            patch = []
            patch.append('pressureInletOutletVelocity')
            patch.append(condition.data[variable])
            return patch
        elif isinstance(condition, api.TotalPressureCondition):
            patch = []
            patch.append('totalPressure')
            patch.append(condition.data[variable])
            return patch
        elif isinstance(condition, api.InletOutletVolumeFraction):
            patch = []
            patch.append('inletOutlet')
            if condition.material == phaseNameToMaterial[phaseNames[0]]:
                vf = condition.data[CUBA.VOLUME_FRACTION]
            else:
                vf = 1 - condition.data[CUBA.VOLUME_FRACTION]
            patch.append(vf)
            return patch
        else:
            error_str = "Boundary condition not supported:\n"
            error_str += " condition: {}\n"
            error_str += " variable: {}"
            raise ValueError(error_str.format(condition.__class__.__name__,
                                              variable))

    elif isinstance(condition, api.FreeSlipVelocity):
        return 'slip'
    elif isinstance(condition, api.WettingAngle):
        patch = []
        patch.append('wettingAngle')
        patch.append(condition.contact_angle)
        return patch
    elif isinstance(condition, api.Dirichlet):
        patch = []
        patch.append('fixedValue')
        variable = get_condition_variable(condition, solver)
        if variable == CUBA.VOLUME_FRACTION:
            if condition.material == phaseNameToMaterial[phaseNames[0]]:
                vf = condition.data[CUBA.VOLUME_FRACTION]
            else:
                vf = 1 - condition.data[CUBA.VOLUME_FRACTION]
            patch.append(vf)
        else:
            patch.append(condition.data[variable])
        return patch
    elif isinstance(condition, api.Neumann):
        # for dynamic pressure use fixedFluxPressure
        variable = get_condition_variable(condition, solver)
        dyn_pres = variable == CUBA.DYNAMIC_PRESSURE
        if dyn_pres:
            return 'fixedFluxPressure'
        else:
            return 'zeroGradient'
    elif isinstance(condition, api.EmptyCondition):
        return 'empty'
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
    for boundary in cuds.iter(item_type=CUBA.BOUNDARY):
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

    for boundary in cuds.iter(item_type=CUBA.BOUNDARY):
        patch = None
        for condition in boundary.condition:
            variable = get_condition_variable(condition, solver)
            if variable == CUBA.VELOCITY:
                patch = get_foam_boundary_condition(
                    condition, mesh._foamPhaseNameToMaterial, solver)
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

    for boundary in cuds.iter(item_type=CUBA.BOUNDARY):
        patch = None
        for condition in boundary.condition:
            variable = get_condition_variable(condition, solver)
            if variable == ID_pressure:
                patch = get_foam_boundary_condition(
                    condition, mesh._foamPhaseNameToMaterial, solver)
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
        elif isinstance(patch, list) and\
                patch[0] == "totalPressure":
            myDict = myDict + "\t type \t totalPressure;\n"
            myDict = myDict + "\t p0 \t uniform " + str(patch[1]) + ";\n"
            myDict = myDict + "\t U  \t U;\n"
            myDict = myDict + "\t phi \t phi;\n"
            myDict = myDict + "\t rho \t rho;\n"
            myDict = myDict + "\t psi \t none;\n"
            myDict = myDict + "\t gamma \t 1;\n"
            myDict = myDict + "\t value \t uniform 0;"
        myDict = myDict + "}\n"

    foamface.setBC(mesh.name, name_pressure, myDict)

    if solver == 'driftFluxFoam' or solver == 'interFoam':
        myDict = ""

        for boundary in cuds.iter(item_type=CUBA.BOUNDARY):
            patch = None
            for condition in boundary.condition:
                variable = get_condition_variable(condition, solver)
                if variable == CUBA.VOLUME_FRACTION or \
                        variable == CUBA.CONTACT_ANGLE:
                    patch = get_foam_boundary_condition(
                        condition, mesh._foamPhaseNameToMaterial, solver)
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
        foamface.setBC(mesh.name,
                       dataNameMap[CUBA.VOLUME_FRACTION] + '.' + phaseNames[0],
                       myDict)


def is_empty(generator):
    if generator is None:
        return True
    try:
        generator.next()
        return False
    except StopIteration:
        return True


def not_empty(generator):
    if generator is None:
        return False
    try:
        generator.next()
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

    if not_empty(cuds.iter(item_type=CUBA.SOLVER_PARAMETER)):
        for sp in cuds.iter(item_type=CUBA.SOLVER_PARAMETER):
            if CUBA.STEADY_STATE in sp.data:
                solver = 'simpleFoam'
                return solver
    if is_empty(cuds.iter(item_type=CUBA.CFD)):
        error_str = "CFD physics model not present in the cuds"
        raise ValueError(error_str)
    cfd = get_first(cuds.iter(item_type=CUBA.CFD))
    laminar = isinstance(cfd.turbulence_model, api.LaminarFlowModel)
    mixture_model = not_empty(cuds.iter(item_type=CUBA.MIXTURE_MODEL))

    vof = not_empty(cuds.iter(item_type=CUBA.FREE_SURFACE_MODEL))

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
