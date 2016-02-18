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

#include "incompressibleTwoPhaseInteractingMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseInteractingMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseInteractingMixture::
incompressibleTwoPhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),
    muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            subDict(phase1Name_),
            U,
            phi
        )
    ),

    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rhod_("rho", dimDensity, muModel_->viscosityProperties().lookup("rho")),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties().lookup("rho")),
    dd_
    (
        "d",
        dimLength,
        muModel_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    alphaMax_(muModel_->viscosityProperties().lookupOrDefault("alphaMax", 1.0)),

    U_(U),
    phi_(phi),

    mu_
    (
        IOobject
        (
            "mu",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}

Foam::incompressibleTwoPhaseInteractingMixture::
incompressibleTwoPhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    Istream& Is
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
		Is
    ),
    twoPhaseMixture(U.mesh(), *this, true),

    muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            subDict(phase1Name_),
            U,
            phi
        )
    ),

    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rhod_("rho", dimDensity, muModel_->viscosityProperties().lookup("rho")),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties().lookup("rho")),
    dd_
    (
        "d",
        dimLength,
        muModel_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    alphaMax_(muModel_->viscosityProperties().lookupOrDefault("alphaMax", 1.0)),

    U_(U),
    phi_(phi),

    mu_
    (
        IOobject
        (
            "mu",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}

Foam::incompressibleTwoPhaseInteractingMixture::
incompressibleTwoPhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    IOdictionary& dict
)
:
    IOdictionary
    (
        dict
    ),
    twoPhaseMixture(U.mesh(), *this, true),

    muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            subDict(phase1Name_),
            U,
            phi
        )
    ),

    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rhod_("rho", dimDensity, muModel_->viscosityProperties().lookup("rho")),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties().lookup("rho")),
    dd_
    (
        "d",
        dimLength,
        muModel_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    alphaMax_(muModel_->viscosityProperties().lookupOrDefault("alphaMax", 1.0)),

    U_(U),
    phi_(phi),

    mu_
    (
        IOobject
        (
            "mu",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::incompressibleTwoPhaseInteractingMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            muModel_().read(subDict(phase1Name_))
         && nucModel_().read(subDict(phase2Name_))
        )
        {
            muModel_->viscosityProperties().lookup("rho") >> rhod_;
            nucModel_->viscosityProperties().lookup("rho") >> rhoc_;

            dd_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel_->viscosityProperties().lookupOrDefault("d", 0)
            );

            alphaMax_ =
                muModel_->viscosityProperties().lookupOrDefault
                (
                    "alphaMax",
                    1.0
                );

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
