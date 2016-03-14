if(stressModelType != word("fromMesoscale")){
    // do SigmaPrime = Sigma_mu-Sigma
    //microscale tensor is selected as the same in macroscale in order to obtain SigmaPrime=0
    Sigma_mu = (mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U)));
}
//calculating SigmaPrime
SigmaPrime = Sigma_mu - (mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U)));


