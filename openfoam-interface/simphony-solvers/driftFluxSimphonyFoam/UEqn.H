    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + fvc::div(UdmModel.tauDm())
      + fvc::div(tauPrime)
      //+ turbulence->divDevRhoReff(U)
      - fvm::laplacian(mixture.mu(), U)
      - fvc::div(mixture.mu()*dev(T(fvc::grad(U))))
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
