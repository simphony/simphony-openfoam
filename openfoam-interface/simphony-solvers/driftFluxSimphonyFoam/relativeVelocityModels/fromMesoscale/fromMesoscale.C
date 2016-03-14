/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "fromMesoscale.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(fromMesoscale, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, fromMesoscale, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::fromMesoscale::fromMesoscale
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
  relativeVelocityModel(dict, mixture),
  Vr_(const_cast< GeometricField < vector, fvPatchField, volMesh > &>(mixture.U().mesh().lookupObject< GeometricField <vector, fvPatchField, volMesh> >(word("Vr"))))
  /*
	Vr_
	(
        IOobject
        (
            "Vr",
            mixture.U().time().timeName(),
            mixture.U().mesh(),
	    IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mixture.U().mesh()
	)*/
{}

Foam::relativeVelocityModels::fromMesoscale::fromMesoscale
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const volVectorField& Vr
)
:
    relativeVelocityModel(dict, mixture),
	Vr_
	(
		Vr
    )
 
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::fromMesoscale::~fromMesoscale()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::fromMesoscale::correct()
{
    Udm_ =
        (rhoc_/rho())
      *Vr_;
}


// ************************************************************************* //
