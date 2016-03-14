S = fvc::grad(U) + T(fvc::grad(U));
Stress = mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U));
gradAlpha1 = fvc::grad(alpha1);
volScalarField betac((1-alpha1)*mixture.rhoc());
volScalarField betad(alpha1*mixture.rhod());
Vr = (betac+betad)/betac*UdmModel.Udm();
