/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "shearStressPowerLawSlipVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "laminar.H"
#include "simphonyInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::
shearStressPowerLawSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(p, iF),
  beta_(0),
  rho_(0),
  n_(1)
{}


Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::
shearStressPowerLawSlipVelocityFvPatchVectorField
(
    const shearStressPowerLawSlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
  fixedValueFvPatchVectorField(ptf, p, iF, mapper),
  beta_(ptf.beta_),
  rho_(ptf.rho_),
  n_(ptf.n_)
{}


Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::
shearStressPowerLawSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
  fixedValueFvPatchVectorField(p, iF),
  beta_(readScalar(dict.lookup("beta"))),
  rho_(readScalar(dict.lookup("rho"))),
  n_(readScalar(dict.lookup("n")))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::
shearStressPowerLawSlipVelocityFvPatchVectorField
(
    const shearStressPowerLawSlipVelocityFvPatchVectorField& mwvpvf
)
:
  fixedValueFvPatchVectorField(mwvpvf),
  beta_(mwvpvf.beta_),
  rho_(mwvpvf.rho_),
  n_(mwvpvf.n_)
{}


Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::
shearStressPowerLawSlipVelocityFvPatchVectorField
(
    const shearStressPowerLawSlipVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(mwvpvf, iF),
  beta_(mwvpvf.beta_),
  rho_(mwvpvf.rho_),
  n_(mwvpvf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = dimensionedInternalField().mesh();

    const fvPatch& p = patch();

    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    
    const surfaceScalarField& phi = mesh.lookupObject<surfaceScalarField>("phi");
    dictionary& TPDict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word("transportProperties")));

    std::string TPString = dictionary_to_string(TPDict);

    singlePhaseTransportModel laminarTransport(U,phi,IStringStream
					       (
						TPString
						)()
					       );

 
    autoPtr<incompressible::RASModel> model
      (
       new incompressible::RASModels::laminar(U, phi, laminarTransport,
					      IStringStream
					      (
					       "RASModel            laminar;"
					       "turbulence            off;" 
					       )())
       );

 
    volSymmTensorField devReff(model->devReff());

    word stressModelType(TPDict.lookup("stressModel"));
    if(stressModelType != word("fromMesoscale")){
      const volTensorField& SigmaPrime = mesh.lookupObject<volTensorField>("SigmaPrime");
      dimensionedScalar rho(TPDict.lookup("rho"));

      devReff = devReff - dev(twoSymm(SigmaPrime/rho)); 
    }
   

    // compute the slip velocity

    vectorField Up = p.Cf()*0.0;

    label patchI = p.index();
    const scalarField& magSf = p.magSf();
    const vectorField& Sf = p.Sf();

    forAll(p, facei)
      {
	vector tauw = ((Sf[facei]/magSf[facei]) & devReff.boundaryField()[patchI][facei]);
	scalar magTauw = mag(tauw);
	scalar magPowLaw = pow(rho_*magTauw, n_);
        // cell velocity direction
	label celli = p.faceCells()[facei];
	scalar Umag = mag(U[celli]);
	scalar coef = beta_*magPowLaw/(Umag + SMALL);
	Up[facei] = coef*U[celli];
      }


    // project direction to wall tangent
    scalarField phip
      (
       p.patchField<surfaceScalarField, scalar>(phi)
       );

    const vectorField n(p.nf());
    tmp<scalarField> Un = phip/(magSf + VSMALL);

    vectorField::operator=(Up + n*(Un - (n & Up)));
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::shearStressPowerLawSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
  fixedValueFvPatchVectorField::write(os);
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        shearStressPowerLawSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
