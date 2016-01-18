S = fvc::grad(U) + T(fvc::grad(U));
gradAlpha1 = fvc::grad(alpha1);
Stress = mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U));
