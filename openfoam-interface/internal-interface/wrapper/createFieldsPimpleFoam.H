    Info<< "Creating field p\n" << endl;
    new volScalarField
    (
        IOobject
        (
            "p",
            runTimes[name]->timeName(),
            *runTimes[name],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("p", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0)
    );

    Info<< "Creating field U\n" << endl;
    volVectorField* U = new volVectorField
    (
        IOobject
        (
            "U",
            runTimes[name]->timeName(),
            *runTimes[name],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
    );

    Info<< "Calculating face flux field phi\n" << endl;

    new surfaceScalarField
    (
        IOobject
        (
            "phi",
            runTimes[name]->timeName(),
            *runTimes[name],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(*U) & mesh.Sf()
    );

