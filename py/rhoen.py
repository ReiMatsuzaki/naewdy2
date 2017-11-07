def rhoen(dvr, ):
    """
    compute density of electron and nuclear wave function.
    .   rho(r,R,t) = <Psi(r,R,t)|delta(r)|Psi(r,R,t)>_r
    .   delta(r) = \sum_{ij}E_{ij}phi^*_i(r,R)
    where E_{ij} is excitation operator. Time dependent wave function forms
    .   Psi(r,R,t) = sum_{an} C_{an}(t) g_a(R) Psi_n(r;R)
    .              = sum_{anI} A_{an}(t) g_a(R) C_{In} Phi_I(r;R)
    finally, density becomes
    .   rho(r,R,t) = sum_{abnmIJ} A^*_{an}A_{bm} g_a^*(R)g_b(R) phi_i(r,R)phi_j(r,R)
    .                             C_{In}C_{Jm} a_{IJ}^{ij}(R)
    """
    
