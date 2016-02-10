""" foam_templates

Templates for OpenFOAM -files

"""
import os
from cuba_extension import CUBAExt

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

blockMeshDict = """
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
        type patch;
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


"""

dictionaryTemplates = {'pimpleFoam':
                       {os.path.join('constant', 'transportProperties'):
                        """
transportModel Newtonian;

nu nu     [ 0 2 -1 0 0 0 0 ]     0.0 ;
                       """,
                        os.path.join('constant', 'turbulenceProperties'):
                        """
simulationType laminar;
                        """,
                        os.path.join('constant', 'RASProperties'):
                        """
RASModel            laminar;

turbulence          off;

printCoeffs         on;
                        """,
                        os.path.join('system', 'controlDict'):
                        """
application     pimpleFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         1;

deltaT          0.0001;

writeControl    adjustableRunTime;

writeInterval   0.001;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable yes;


                        """,
                        os.path.join('system', 'decomposeParDict'):
                        """
numberOfSubdomains 4;

method          scotch;

                        """,
                        os.path.join('system', 'fvSolution'):
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
                        os.path.join('system', 'fvSchemes'):
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
                       'simpleFoam':
                       {os.path.join('constant', 'transportProperties'):
                        """
transportModel Newtonian;

nu nu     [ 0 2 -1 0 0 0 0 ]     0.0 ;

                        """,
                        os.path.join('constant', 'turbulenceProperties'):
                        """
simulationType laminar;
                        """,
                        os.path.join('constant', 'RASProperties'):
                        """
RASModel            laminar;

turbulence          off;

printCoeffs         on;
                        """,
                        os.path.join('system', 'controlDict'):
                        """
application simpleFoam;

startFrom startTime;

startTime 0;

stopAt endTime;

endTime 1000;

deltaT 1;

writeControl runTime;

writeInterval 1000;

purgeWrite 0;

writeFormat ascii;

writePrecision 6;

writeCompression no;
timeFormat general;

timePrecision 6;

runTimeModifiable yes;
                        """,
                        os.path.join('system', 'decomposeParDict'):
                        """
numberOfSubdomains 4;

method          scotch;

                        """,
                        os.path.join('system', 'fvSolution'):
                        """

solvers
{
p
{
solver          GAMG;
tolerance       1e-06;
relTol          0.1;
smoother        GaussSeidel;
nPreSweeps      0;
nPostSweeps     2;
cacheAgglomeration true;
nCellsInCoarsestLevel 10;
agglomerator    faceAreaPair;
mergeLevels     1;
}

U
{
solver          smoothSolver;
smoother        GaussSeidel;
nSweeps         2;
tolerance       1e-08;
relTol          0.1;
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

SIMPLE
{
nNonOrthogonalCorrectors 0;
pRefCell        0;
pRefValue       0;

residualControl
{
p               1e-5;
U               1e-5;
nuTilda         1e-5;
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
nuTilda         0.7;
}
}
                        """,
                        os.path.join('system', 'fvSchemes'):
                        """
ddtSchemes
{
default         steadyState;
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
div(phi,nuTilda) bounded Gauss linearUpwind grad(nuTilda);
div((nuEff*dev(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
default         none;
laplacian(nuEff,U) Gauss linear corrected;
laplacian((1|A(U)),p) Gauss linear corrected;
laplacian(DnuTildaEff,nuTilda) Gauss linear corrected;
laplacian(1,p)  Gauss linear corrected;
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
                       {os.path.join('constant', 'transportProperties'):
                        """
phases (phase1 phase2);


phase1
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1e-06;
    rho             rho [ 1 -3 0 0 0 0 0 ] 1000;
 }
phase2
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1.48e-05;
    rho             rho [ 1 -3 0 0 0 0 0 ] 1;
 }
sigma           sigma [ 1 0 -2 0 0 0 0 ] 0.07;

                        """,
                        os.path.join('constant', 'turbulenceProperties'):
                        """
simulationType  laminar;

                        """,
                        os.path.join('constant', 'RASProperties'):
                        """
RASModel        laminar;

turbulence      off;

printCoeffs     on;

                        """,
                        os.path.join('constant', 'g'):
                        """
dimensions      [0 1 -2 0 0 0 0];
value           ( 0 0 0 );
                        """,
                        os.path.join('system', 'controlDict'):
                        """
application     interFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1;

deltaT          0.001;

writeControl    adjustableRunTime;

writeInterval   0.05;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression uncompressed;

timeFormat      general;

timePrecision   6;

runTimeModifiable yes;

adjustTimeStep  off;

maxCo           1;
maxAlphaCo      1;

maxDeltaT       1;

                        """,
                        os.path.join('system', 'decomposeParDict'):
                        """
numberOfSubdomains 4;

method          scotch;

                        """,
                        os.path.join('system', 'fvSolution'):
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
                        os.path.join('system', 'fvSchemes'):
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
    div(rhoPhi,U)  Gauss linearUpwind grad(U);
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
                           {os.path.join('constant', 'transportProperties'):
                               """

phases (phase1 phase2);


phase1
{
    transportModel  BinghamPlastic;
    BinghamPlasticCoeffs
    {
        coeff       0.00023143;
        exponent    179.26;

        BinghamCoeff    0.0005966;
        BinghamExponent 1050.8;
        BinghamOffset   0;

        muMax       10;
    }
    plasticCoeffs
    {
        coeff       0.00023143;
        exponent    179.26;

        BinghamCoeff    0.0005966;
        BinghamExponent 1050.8;
        BinghamOffset   0;

        muMax       10;
    }
     rho             1000;
 }
phase2
{
    transportModel  Newtonian;
    nu              1.48e-05;
    rho             1;
 }

stressModel standard;
relativeVelocityModel simple;

simpleCoeffs
{
    V0              (0 0 0);
    a               285.84;
    a1              0.1;
    residualAlpha   0;
}
generalCoeffs
{
    V0              (0 0 0);
    a               285.84;
    a1              0.1;
    residualAlpha   0;
}
fromMesoscaleCoeffs
{

}
                        """,
                            os.path.join('constant', 'turbulenceProperties'):
                                """
simulationType  laminar;

                        """,
                            os.path.join('constant', 'g'):
                                """
dimensions      [0 1 -2 0 0 0 0];
value           ( 0 0 0 );
                        """,
                            os.path.join('system', 'controlDict'):
                                """
application     driftFluxSimphonyFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         0.1;

deltaT          0.1;

writeControl    adjustableRunTime;

writeInterval   10000;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression no;

timeFormat      general;

timePrecision   6;

runTimeModifiable yes;

adjustTimeStep  on;

maxCo           5;

maxDeltaT       1;
                        """,
                            os.path.join('system', 'decomposeParDict'):
                                """
numberOfSubdomains 4;

method          scotch;

                        """,
                            os.path.join('system', 'fvSolution'):
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
                            os.path.join('system', 'fvSchemes'):
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
    default             none;
    div(tauPrime)       Gauss linear;
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

                        """}
                       }


dataTemplates = {'pimpleFoam': {'U': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}
                    """, 'p': """
dimensions [ 0 2 -2 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """}, 'simpleFoam': {'U': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}
                    """, 'p': """
dimensions [ 0 2 -2 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """}, 'interFoam': {'U': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}

                    """,
                                        'alpha.phase1': """
dimensions [ 0 0 0 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """,
                                        'p_rgh': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform 0;

boundaryField
{

}

                    """}, 'driftFluxSimphonyFoam': {'U': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}

                    """, 'p_rgh': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform 0;

boundaryField
{

}

                    """, 'alpha.phase1': """
dimensions [ 0 0 0 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """, 'Vdj': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}
                    """, 'Sigma': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{

}
                    """, 'gradAlpha1': """
dimensions [0 -1 0 0 0 0 0];

internalField uniform (0 0 0);

boundaryField
{

}
                    """, 'Stress': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{

}
                    """, 'S': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{

}

                    """}, 'generalFoam': {'U': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}

                    """, 'p_rgh': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform 0;

boundaryField
{

}

                    """, 'p': """
dimensions [ 0 2 -2 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}

                    """, 'alpha.phase1': """
dimensions [ 0 0 0 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """, 'Vdj': """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}
                    """, 'Sigma': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{

}
                    """, 'gradAlpha1': """
dimensions [0 -1 0 0 0 0 0];

internalField uniform (0 0 0);

boundaryField
{

}
                    """, 'Stress': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{

}
                    """, 'S': """
dimensions [1 -1 -2 0 0 0 0];

internalField uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{

}

                    """}


                 }

multiphase_solvers = ("interFoam", "driftFluxSimphonyFoam")


def get_foam_solver(CM):
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
