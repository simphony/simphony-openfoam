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

new volScalarField
(
    IOobject
    (
        "p_rgh",
        runTimes[name]->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("p_rgh", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0)
);

new volVectorField
(
    IOobject
    (
        "U",
        runTimes[name]->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);

new volVectorField
(
    IOobject
    (
        "Vr",
        runTimes[name]->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Vr", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);



new volTensorField
(
    IOobject
    (
        Sigma_muName,
        runTimes[name]->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor(Sigma_muName, dimensionSet(1, -1, -2, 0, 0, 0, 0), tensor::zero)
);


new volScalarField
(
    IOobject
    (
        "alpha.phase1",
        runTimes[name]->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
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
  IOobject::READ_IF_PRESENT,
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
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedTensor("Stress", dimensionSet(1, -1, -2, 0, 0, 0, 0), tensor::zero)
 );

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
        runTimes[name]->constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    TP()
);

IStringStream TuP
    (
     "simulationType laminar;"
    );

new IOdictionary
(
    IOobject
    (
        "turbulenceProperties",
        runTimes[name]->constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    TuP()
);

IStringStream RaP
    (
     "RASModel            laminar;"
    );

new IOdictionary
(
    IOobject
    (
        "RASProperties",
        runTimes[name]->constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    RaP()
);

IStringStream deC
    (
     "numberOfSubdomains 2;"
     "method          scotch;"
    );

new IOdictionary
(
    IOobject
    (
        "decomposeParDict",
        runTimes[name]->system(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    deC()
);

IStringStream grav
    (
     "dimensions      [0 1 -2 0 0 0 0];"
     "value           ( 0 0 0 );"
    );


new IOdictionary
(
    IOobject
    (
        "g",
        runTimes[name]->constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    grav()
);
