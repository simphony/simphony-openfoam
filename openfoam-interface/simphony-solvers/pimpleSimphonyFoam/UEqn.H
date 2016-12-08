// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
 fvm::ddt(U)
 + fvm::div(phi, U)
 - fvc::div(SigmaPrime)/rho
 + turbulence->divDevReff(U)
 // ==
 //  fvOptions(U)
);

UEqn().relax();

//fvOptions.constrain(UEqn());

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    //    fvOptions.correct(U);
}
