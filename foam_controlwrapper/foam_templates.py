""" foam_templates

Templates for OpenFOAm -files

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

timePrecision 6;

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
                        """}}
vectorTemplates = {'simpleFoam':
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
                    """}}
