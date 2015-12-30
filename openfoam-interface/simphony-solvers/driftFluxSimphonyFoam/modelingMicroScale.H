if(1){
    //Re-read tau and do tauprime = tau_m-tau

    //reading from microscale (meanwhile this tensor is selected as the same in macroscale in order to obtain tauPrime=0)
    Sigma = (mixture.mu()*fvc::grad(U) + mixture.mu()*T(fvc::grad(U)));

    //calculating tauPrime
    tauPrime = Sigma - (mixture.mu()*fvc::grad(U) - mixture.mu()*T(fvc::grad(U)));

}
