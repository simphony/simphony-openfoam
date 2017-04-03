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

new volScalarField
(
    IOobject
    (
        "p",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::MUST_WRITE
    ),
    mesh,
    dimensionedScalar("p", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0)
);

new volVectorField
(
    IOobject
    (
        "U",
        runTimes[name]->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::MUST_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);

