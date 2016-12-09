/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //scale pressure
    dictionary& TPDict = const_cast<dictionary&>(mesh.lookupObject<dictionary>(word("transportProperties")));

     dimensionedScalar rho(TPDict.lookup("rho"));
     dimensionedScalar dimlessRho = dimensionedScalar("dimlessRho",dimensionSet(0,0,0,0,0,0,0),rho.value());
     p=p/dimlessRho;

     word stressModelType(TPDict.lookup("stressModel"));
     Info<< "Stressmodel is "<< stressModelType<<"\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	if(stressModelType != word("fromMesoscale")){
	  //microscale tensor is selected as the same in macroscale 
	  // in order to obtain SigmaPrime=0
	  Sigma = Stress;
	}
	SigmaPrime = Sigma - Stress;


        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();


	S = fvc::grad(U) + dev(T(fvc::grad(U)));
	Stress = turbulence->nuEff()*rho*S;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    p = p*dimlessRho;
    p.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
