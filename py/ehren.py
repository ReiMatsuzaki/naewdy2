import numpy as np

from naewdy2.math import eigh_sort

dot = np.dot
exp = np.exp
tr = np.transpose

class Ehren1d(object):
    """
    1 dimensional semiclassical Ehrenfest theory dynamics.
    Ehrenfest force is given as
    .    F_k = -\sum_IJ C_I^* C_J [dH_IJ/dR_k + \sum_K {X^k_IK H_KJ^k - H_IK X^k_KJ} + i \sum_l dot{R}_l[dX^l_IJ/dR_k - dX^k_IJ/dR_l] ]
    In 1 dimensional case, 
    .    F = -\sum_IJ C_I^* C_J [ dH_IJ/dR + \sum_K {X_IK H_KJ - H_IK X_KJ} ]
    """
    def __init__(self, dt, m, q0, p0, c0):
        self.dt = dt
        self.m = m
        self.q = q0
        self.p = p0
        self.c = c0
        
    def update(self, hIJ, dhIJ, xIJ):
        """
        hIJ  : matrix(float)
        dhIJ : matrix(float)
        xIJ  : matrix(float)
        """

        # -- classical part --
        f = -dot(self.c.conj(),
                 dot(dhIJ, self.c) +
                 dot(xIJ, dot(hIJ, self.c)) -
                 dot(hIJ, dot(xIJ, self.c)) )
        dp = f
        dq = self.p/self.m

        # -- quantum part --
        heff = hIJ - 1.0j*dq*xIJ
        (e, u) = eigh_sort(heff)
        uH = tr(u.conj())

        # -- update variables --
        self.c = dot(uH, self.c)
        self.c = exp(-1.0j*e*self.dt) * self.c
        self.c = dot(u, self.c)
        
        self.p += self.dt * f.real
        self.q += self.dt * self.p/self.m
        

        
        

        
