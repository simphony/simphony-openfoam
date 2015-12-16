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

Class
    Foam::fromMesoscale

Description
    fromMesoscale relative velocity model

SourceFiles
    fromMesoscale.C

\*---------------------------------------------------------------------------*/

#ifndef fromMesoscale_H
#define fromMesoscale_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class fromMesoscale Declaration
\*---------------------------------------------------------------------------*/

class fromMesoscale
:
    public relativeVelocityModel
{
    // Private data

        //- Relative velocity
	volVectorField Vdj_;


public:

    //- Runtime type information
    TypeName("fromMesoscale");


    // Constructors

        //- Construct from components
        fromMesoscale
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );

        //- Construct from components
        fromMesoscale
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
			const volVectorField& Vdj
        );

    //- Destructor
    ~fromMesoscale();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
