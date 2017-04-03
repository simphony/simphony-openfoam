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
//std::string Sigma_muName = "Sigma_mu";
std::string Sigma_muName = "Sigma";

Info<< "Creating field p_rgh\n" << endl;
new volScalarField
(
    IOobject
    (
        "p_rgh",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p_rgh", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0)
);

Info<< "Creating field U\n" << endl;
volVectorField* U = new volVectorField
(
    IOobject
    (
        "U",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);

Info<< "Creating field Vr\n" << endl;
new volVectorField
(
    IOobject
    (
        "Vr",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Vr", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);



Info<< "Creating field Sigma" << endl;
new volTensorField
(
    IOobject
    (
        Sigma_muName,
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor(Sigma_muName, dimensionSet(1, -1, -2, 0, 0, 0, 0), tensor::zero)
);

Info<< "Calculating face flux field phi\n" << endl;

new surfaceScalarField
(
    IOobject
    (
        "phi",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(*U) & mesh.Sf()
);

Info<< "Creating field alpha.phase1\n" << endl;
new volScalarField
(
    IOobject
    (
        "alpha.phase1",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha.phase1", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
);

new volTensorField
(
    IOobject
    (
        "S",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("S", dimensionSet(0, 0, -1, 0, 0, 0, 0), tensor::zero)
);

new volVectorField
(
 IOobject
 (
  "gradAlpha1",
  runTimes[name]->timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("gradAlpha1", dimensionSet(0, -1, 0, 0, 0, 0, 0), vector::zero)
 );

new volTensorField
(
 IOobject
 (
  "Stress",
  runTimes[name]->timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedTensor("Stress", dimensionSet(1, -1, -2, 0, 0, 0, 0), tensor::zero)
 );

Info<< "Creating transportProperties\n" << endl;
IStringStream TP
    (
    "phases (phase1 phase2);\n"
    "phase1\n"
    "{\n"
    "    transportModel  dummyViscosity;\n"
    "    nu              0.0;\n"
    "    rho             0.0;\n"
    "}\n"
    "phase2\n"
    "{\n"
    "    transportModel  Newtonian;\n"
    "    nu              0.0;\n"
    "    rho             0.0;\n"
    "}\n"
    "relativeVelocityModel simple;\n"
    "simpleCoeffs\n"
    "{\n"
    "    V0              (0 0 0);\n"
    "    a               0.0;\n"
    "    a1              0.0;\n"
    "    residualAlpha   0;\n"
    "}\n"
    "generalCoeffs\n"
    "{\n"
    "    V0              (0 0 0);\n"
    "    a               0.0;\n"
    "    a1              0.0;\n"
    "    residualAlpha   0;\n"
    "}\n"
    );

new IOdictionary
(
    IOobject
    (
        "transportProperties",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    TP()
);

IStringStream RaP
    (
     "RASModel            laminar;"
     "turbulence            off;"
    );

new IOdictionary
(
    IOobject
    (
        "RASProperties",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    RaP()
);


/*
IStringStream grav
    (
     "dimensions      [0 1 -2 0 0 0 0];\n"
     "value           ( 0 0 0 );\n"
    );
*/
new uniformDimensionedVectorField
(
 IOobject
 (
  "g",
  runTimes[name]->constant(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 dimensionedVector("g", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector::zero)
 );


