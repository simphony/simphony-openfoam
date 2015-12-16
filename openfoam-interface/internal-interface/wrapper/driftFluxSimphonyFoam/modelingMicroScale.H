if(1){
    //Re-read Sigma_mu and do SigmaPrime = Sigma_mu-Sigma

    //reading from microscale (meanwhile this tensor is selected as the same in macroscale in order to obtain SigmaPrime=0)
    Sigma_mu = (mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U)));

    //calculating SigmaPrime
    SigmaPrime = Sigma_mu - (mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U)));

}
