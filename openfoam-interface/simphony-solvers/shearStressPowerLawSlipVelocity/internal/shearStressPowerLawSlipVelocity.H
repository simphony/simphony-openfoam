/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::shearStressPowerLawSlipVelocityFvPatchVectorField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a velocity condition for cases with
    slip velocity relative to wall shear stress.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        U            | velociy field name      | no          | U
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            shearStressPowerLawSlipVelocity;
        beta            3.1e-3;
	rho             250.0;
        n               1.16;
        value           uniform (0 0 0);
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchVectorField

SourceFiles
    shearStressPowerLawSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef shearStressPowerLawSlipVelocityFvPatchVectorField_H
#define shearStressPowerLawSlipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class shearStressPowerLawSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class shearStressPowerLawSlipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

       // power law multiplier for slip velocity ratio to wall shear stress
       scalar beta_;
       // density
       scalar rho_;
       // power law power 
       scalar n_;

public:

    //- Runtime type information
    TypeName("shearStressPowerLawSlipVelocity");

    // Constructors

        //- Construct from patch and internal field
        shearStressPowerLawSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        shearStressPowerLawSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given shearStressPowerLawSlipVelocityFvPatchVectorField
        //  onto a new patch
        shearStressPowerLawSlipVelocityFvPatchVectorField
        (
            const shearStressPowerLawSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        shearStressPowerLawSlipVelocityFvPatchVectorField
        (
            const shearStressPowerLawSlipVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new shearStressPowerLawSlipVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        shearStressPowerLawSlipVelocityFvPatchVectorField
        (
            const shearStressPowerLawSlipVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new shearStressPowerLawSlipVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        const scalar& beta() const
        {
	  return beta_;
        }
        const scalar& rho() const
        {
	  return rho_;
        }
        const scalar& n() const
        {
	  return n_;
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
