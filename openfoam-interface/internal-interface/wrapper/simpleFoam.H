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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "simpleControl/simpleControlMod.H"

#include "initContinuityErrs.H"

simpleControlMod simple(mesh);

std::string pName = "p";
std::string UName = "U";
std::string Sigma_muName = "Sigma";
std::string SName = "S";
std::string StressName = "Stress";


volScalarField& p = find_Data<scalar>(mesh,pName);

volVectorField& U = find_Data<vector>(mesh,UName);
surfaceScalarField& phi = const_cast<surfaceScalarField&>(mesh.lookupObject<surfaceScalarField>(word("phi"))); 
volTensorField& Sigma_mu = find_Data<tensor>(mesh,Sigma_muName);
volTensorField& S = find_Data<tensor>(mesh,SName);
volTensorField& Stress = find_Data<tensor>(mesh,StressName);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

dictionary& TPDict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word("transportProperties")));
std::string TPString = dictionary_to_string(TPDict);

singlePhaseTransportModel laminarTransport(U,phi,IStringStream
    (
     TPString
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


word stressModelType(TPDict.lookup("stressModel"));
Info<< "Stressmodel is "<< stressModelType<<"\n" << endl;

dimensionedScalar rho(TPDict.lookup("rho"));

Info<< "Creating field SigmaPrime\n" << endl;
volTensorField SigmaPrime
(
    IOobject
    (
        "SigmaPrime",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Sigma_mu - Stress
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "\nStarting time loop\n" << endl;

//scale pressure
dimensionedScalar dimlessRho = dimensionedScalar("dimlessRho",dimensionSet(0,0,0,0,0,0,0),rho.value());

p=p/dimlessRho;

U.storePrevIter();
p.storePrevIter();

scalar cumulativeContErr = 0;

while (simple.loop())
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // --- Pressure-velocity SIMPLE corrector
    {
    // Momentum predictor

      if(stressModelType != word("fromMesoscale")){
	// do SigmaPrime = Sigma_mu-Sigma
	//microscale tensor is selected as the same in macroscale 
	// in order to obtain SigmaPrime=0
	Sigma_mu = Stress;
      }

      SigmaPrime = Sigma_mu - Stress;

     tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvc::div(SigmaPrime)/rho
	+ turbulence->divDevReff(U)
	//      ==
	//        fvOptions(U)
    );

    UEqn().relax();

    //    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    //    fvOptions.correct(U);
 
    volScalarField rAU(1.0/UEqn().A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();
    UEqn.clear();

    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf());

    //    fvOptions.makeRelative(phiHbyA);

    adjustPhi(phiHbyA, U, p);

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
    //    fvOptions.correct(U);

    }

    turbulence->correct();

    S = fvc::grad(U) + dev(T(fvc::grad(U)));
    Stress = laminarTransport.nu()*rho*S;


    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;

  }
p=p*dimlessRho;

Info<< "End\n" << endl;



// ************************************************************************* //
