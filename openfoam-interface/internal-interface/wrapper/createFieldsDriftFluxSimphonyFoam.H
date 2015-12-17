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
    dimensionedScalar("p", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0)
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

Info<< "Creating field Sigma_mu" << endl;
new volTensorField
(
    IOobject
    (
        "Sigma_mu",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("Sigma_mu", dimensionSet(1, -1, -2, 0, 0, 0, 0), tensor::zero)
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

Info<< "Creating field alpha1\n" << endl;
new volScalarField
(
    IOobject
    (
        "alpha1",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha1", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
);


Info<< "Creating transportProperties\n" << endl;
IStringStream TP
    (
    "phases (heavy light);\n"
    "heavy\n"
    "{\n"
    "    transportModel  dummyViscosity;\n"
    "    nu              1e-04;\n"
    "    rho             1996;\n"
    "}\n"
    "light\n"
    "{\n"
    "    transportModel  Newtonian;\n"
    "    nu              1.7871e-06;\n"
    "    rho             996;\n"
    "}\n"
    "relativeVelocityModel simple;\n"
    "\"(simple|general)Coeffs\"\n"
    "{\n"
    "    V0              (0 -0.002198 0);\n"
    "    a               285.84;\n"
    "    a1              0.1;\n"
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