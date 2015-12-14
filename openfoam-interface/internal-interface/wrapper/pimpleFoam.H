/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2015 SimPhoNy -project
   \\/     M anipulation  |
-------------------------------------------------------------------------------


License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


\*---------------------------------------------------------------------------*/

#include "initContinuityErrs.H"
pimpleControl pimple(mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "\nStarting time loop\n" << endl;

std::string pName = "p";
std::string UName = "U";
volScalarField& p = find_Data<scalar>(runTime,pName);
volVectorField& U = find_Data<vector>(runTime,UName);
surfaceScalarField& phi = const_cast<surfaceScalarField&>(runTime.db().lookupObject<surfaceScalarField>(word("phi"))); 


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U,phi,IStringStream
    (
        "transportModel  Newtonian;"
        "nu              nu [ 0 2 -1 0 0 0 0 ] 0.01;"
    )()
);

IStringStream RASPropertiesDict
(
    "RASModel        laminar;"
    "turbulence      off;"
    "printCoeffs     off;"
);


autoPtr<incompressible::RASModel> turbulence
(
    new incompressible::RASModels::laminar
    (
        U, 
        phi, 
        laminarTransport, 
        RASPropertiesDict()
    )
);



Info<< "\nStarting time loop\n" << endl;

U.storePrevIter();
p.storePrevIter();

while (runTime.run())
{
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"

    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        // Solve the Momentum equation

        tmp<fvVectorMatrix> UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + turbulence->divDevReff(U)
        // ==
        //    fvOptions(U)
        );

        UEqn().relax();

        //fvOptions.constrain(UEqn());

        volScalarField rAU(1.0/UEqn().A());

        if (pimple.momentumPredictor())
        {
            solve(UEqn() == -fvc::grad(p));

            //fvOptions.correct(U);
        }

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
            
            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn().H();
            
            if (pimple.nCorrPISO() <= 1)
            {
                UEqn.clear();
            }

            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + rAUf*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            //fvOptions.relativeFlux(phiHbyA);

            // Non-orthogonal pressure corrector loop
            while (pimple.correctNonOrthogonal())
            {
                // Pressure corrector
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
                );
            
                pEqn.setReference(pRefCell, pRefValue);
            
                pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
            
                if (pimple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();
            
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();


            //fvOptions.correct(U);
        }

        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }

    //runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}

Info<< "End\n" << endl;



// ************************************************************************* //
