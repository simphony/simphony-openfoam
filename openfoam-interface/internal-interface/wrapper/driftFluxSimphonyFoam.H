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

pimpleControl pimple(mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "\nStarting time loop\n" << endl;

std::string p_rghName = "p_rgh";
std::string UName = "U";
std::string Sigma_muName = "Sigma";
std::string SName = "S";
std::string gradAlpha1Name = "gradAlpha1";
std::string StressName = "Stress";
std::string VrName = "Vr";
std::string alphaName = "alpha.phase1";

volScalarField& p_rgh = find_Data<scalar>(mesh,p_rghName);
volVectorField& U = find_Data<vector>(mesh,UName);
volTensorField& Sigma_mu = find_Data<tensor>(mesh,Sigma_muName);
volTensorField& S = find_Data<tensor>(mesh,SName);
volVectorField& gradAlpha1 = find_Data<vector>(mesh,gradAlpha1Name);
volTensorField& Stress = find_Data<tensor>(mesh,StressName);
volVectorField& Vr = find_Data<vector>(mesh,VrName);

volScalarField& alpha1 = find_Data<scalar>(mesh,alphaName);

surfaceScalarField& phi = const_cast<surfaceScalarField&>(mesh.lookupObject<surfaceScalarField>(word("phi"))); 


IOdictionary& transportProperties = const_cast<IOdictionary&>(mesh.lookupObject<IOdictionary>(word("transportProperties"))); 



incompressibleTwoPhaseInteractingMixture mixture(U, phi, transportProperties);
// Relative Velocity
autoPtr<relativeVelocityModel> UdmModelPtr
(
    relativeVelocityModel::New
    (
        transportProperties,
        mixture
    )
);


relativeVelocityModel& UdmModel(UdmModelPtr());


word stressModelType(transportProperties.lookup("stressModel"));
Info<< "Stressmodel is "<< stressModelType<<"\n" << endl;


//volScalarField& alpha1(mixture.alpha1());

volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rhod();
const dimensionedScalar& rho2 = mixture.rhoc();

Info<<"Density phase 1 "<<rho1<<endl;
Info<<"Density phase 2 "<<rho2<<endl;

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
    Sigma_mu - (mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U)))
);

// Mixture density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.rho()
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// g 
uniformDimensionedVectorField& g = const_cast<uniformDimensionedVectorField&>(mesh.lookupObject<uniformDimensionedVectorField>(word("g"))); 

Info<<"g "<<g<<endl; 
Info<< "Calculating field (g.h)f\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("gh", g & mesh.Cf());

volScalarField p
(
    IOobject
    (
        "p",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

// MULES Correction
tmp<surfaceScalarField> tphiAlphaCorr0;

Info<< "\nStarting time loop\n" << endl;


U.storePrevIter();
p_rgh.storePrevIter();
alpha1.storePrevIter();


scalar cumulativeContErr = 0;


while (runTime.run())
{
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"

    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;


    #include "modelingMicroScale.H"


    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        #include "alphaControls.H"

        UdmModel.correct();

	//	Info <<UdmModel.Udm()<<endl;

        #include "alphaEqnSubCycle.H"

        mixture.correct();

        #include "UEqn.H"


        // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }

    }

    #include "calcStrain.H"

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

}

Info<< "End\n" << endl;


// ************************************************************************* //
