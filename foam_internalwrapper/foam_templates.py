""" foam_templates

Templates for OpenFOAM -files

"""
import os

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

dictionaryMaps = {'pimpleFoam':
                       {'transportProperties':
                          {'transportModel':'Newtonian',
                           'nu nu     [ 0 2 -1 0 0 0 0 ]': '0.0'},
                       'turbulenceProperties':
                          {'simulationType':'laminar'},
                       'RASProperties':
                          {'RASModel':'laminar',
                          'turbulence':'off',
                          'printCoeffs':'on'},
                       'controlDict':
                          {'application':'pimpleFoam',
                           'startFrom':'startTime',
                           'startTime':'0',
                           'stopAt':'endTime',
                           'endTime':'1',
                           'deltaT':'1',
                           'writeControl':'timeStep',
                           'writeInterval':'10000',
                           'purgeWrite':'0',
                           'writeFormat':'ascii',
                           'writePrecision':'6',
                           'writeCompression':'no',
                           'timeFormat':'general',
                           'runTimeModifiable':'yes'}
                       }
                      }

                          
                        
dictionaryTemplates = {'simpleFoam':
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

writeControl timeStep;

writeInterval 1000;

purgeWrite 0;

writeFormat ascii;

writePrecision 6;

writeCompression no;
timeFormat general;

runTimeModifiable yes;
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
			'pimpleFoam':
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

endTime 1;

deltaT 1;

writeControl timeStep;

writeInterval 10000;

purgeWrite 0;

writeFormat ascii;

writePrecision 6;

writeCompression no;
timeFormat general;

runTimeModifiable yes;
                        """,
                        os.path.join('system', 'fvSolution'):
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
                       'interFoam':
                       {os.path.join('constant', 'transportProperties'):
                        """
twoPhase
{
    transportModel  twoPhase;
    phase1          phase1;
    phase2          phase2;
}

phase1
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1e-06;
}
phase2
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1.48e-05;
}

                        """,
                        os.path.join('constant', 'turbulenceProperties'):
                        """
simulationType  laminar;

                        """,
                        os.path.join('constant', 'RASProperties'):
                        """
RASModel        laminar;

turbulence      on;

printCoeffs     on;

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

adjustTimeStep  on;

maxCo           0.2;
maxAlphaCo      0.2;

maxDeltaT       1;

                        """,
                        os.path.join('system', 'fvSolution'):
                        """
solvers
{
    pcorr
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-10;
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
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-07;
        relTol          0;
    }

    "(U|k|omega)"
    {
        solver          PBiCG;
        preconditioner  DILU;
        tolerance       1e-06;
        relTol          0;
    }

    "(U|k|omega)Final"
    {
        solver          PBiCG;
        preconditioner  DILU;
        tolerance       1e-08;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor no;
    nCorrectors     3;
    nNonOrthogonalCorrectors 0;
    nAlphaCorr      1;
    nAlphaSubCycles 4;
    cAlpha          2;
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
    div(rho*phi,U)  Gauss linear;
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss interfaceCompression;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,R)      Gauss upwind;
    div(R)          Gauss linear;
    div(phi,nuTilda) Gauss upwind;
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
    alpha;
}

                        """}}

vectorTemplates = {'simpleFoam':
                   {os.path.join('0', 'U'):
                    """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}
                    """},
                    'pimpleFoam':
                   {os.path.join('0', 'U'):
                    """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}
                    """},
                   'interFoam':
                   {os.path.join('0', 'U'):
                    """
dimensions [ 0 1 -1 0 0 0 0 ];

internalField uniform (0 0 0);

boundaryField
{

}

                    """}}

scalarTemplates = {'simpleFoam':
                   {os.path.join('0', 'p'):
                    """
dimensions [ 0 2 -2 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """},
                    'pimpleFoam':
                   {os.path.join('0', 'p'):
                    """
dimensions [ 0 2 -2 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}
                    """},
                   'interFoam':
                   {os.path.join('0', 'p_rgh'):
                    """
dimensions [ 0 2 -2 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}

                    """,
                    os.path.join('0', 'alpha1'):
                    """
dimensions [ 0 0 0 0 0 0 0 ];

internalField uniform 0;

boundaryField
{

}

                    """}}
