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

#include "twoPhaseMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMixture::twoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    phase2Name_(wordList(dict.lookup("phases"))[1]),

    alpha1_
    (
        const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("alpha.phase1"))
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_
    )
{

}

Foam::twoPhaseMixture::twoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict,
    bool not_read
)
:
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    phase2Name_(wordList(dict.lookup("phases"))[1]),
    alpha1_
    (
        const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("alpha.phase1"))
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_
    )
{
}




// ************************************************************************* //
