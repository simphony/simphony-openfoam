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
simpleControl simple(mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "\nStarting time loop\n" << endl;

std::string pName = "p";
std::string UName = "U";
volScalarField& p = find_Data<scalar>(runTime,pName);
volVectorField& U = find_Data<vector>(runTime,UName);
surfaceScalarField& phi = const_cast<surfaceScalarField&>(runTime.db().lookupObject<surfaceScalarField>(word("phi"))); 


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

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


U.storePrevIter();
p.storePrevIter();

while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // --- Pressure-velocity SIMPLE corrector
    {
		// Momentum predictor

		tmp<fvVectorMatrix> UEqn
		(
			fvm::div(phi, U)
		  + turbulence->divDevReff(U)
		  //==
		  //fvOptions(U)
		);

		UEqn().relax();

		//fvOptions.constrain(UEqn());

		solve(UEqn() == -fvc::grad(p));

		//fvOptions.correct(U);

		{
            volScalarField rAU(1.0/UEqn().A());
            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn().H();
            UEqn.clear();

            surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf());
            adjustPhi(phiHbyA, U, p);

            //fvOptions.relativeFlux(phiHbyA);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            //fvOptions.correct(U);
        }
	}

    turbulence->correct();

    //runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}

Info<< "End\n" << endl;



// ************************************************************************* //
