import numpy as np

from naewdy2.math import uni_inte

tr = np.transpose
dot = np.dot

def idx(nel, a, i):
    return a*nel+i
    
class PBE1d1p(object):
    def __init__(self, dt, m, dvr, c0, uni_inte_opts={}):
        self.dt = dt
        self.dvr = dvr      # DVR object
        (self.nnuc, self.nel) = c0.shape
        if(dvr.num != self.nnuc):
            raise RuntimeError("dvr size and given c0 is mismatch")

        n = self.nel*self.nnuc
        self.c = np.reshape(c0, n)
        norm = np.sqrt(dot(self.c.conj(), self.c))
        self.c /= norm

        self.uni_inte_opts = uni_inte_opts

        self.tmat = -1/(2*m)*self.dvr.dmat(2)

    def update(self, dq, hIJq, xIJq):
        """
        hIJq : matrix(float)
        .      Electronic Hamiltonian evaluated at q
        xIJq : matrix(float)
        .      Derivative coupling evaluated at q
        """

        h = self.calc_hmat(dq,hIJq,xIJq)
        self.c = uni_inte(h, self.dt, self.c, self.uni_inte_opts)
        
    def calc_hmat(self, dq, hIJq, xIJq):
        """
        <gA.PhiI|H-i.Dt|gB.PhiJ> = d_AB H^{el}_IJ  +  d_IJ T_AB  -i d_AB X_IJ dot{Q}
        """
        n = self.nel*self.dvr.num
        h = np.zeros((n,n), dtype=complex)
                
        elid = np.identity(self.nel)
        for a in range(self.nnuc):
            a0 = self.idx(a,0)
            a1 = self.idx(a+1,0)
            h[a0:a1,a0:a1] = hIJq -1.0j*dq*xIJq
            for b in range(self.nnuc):
                b0 = self.idx(b,0)
                b1 = self.idx(b+1,0)
                h[a0:a1,b0:b1] += self.tmat[a,b]*elid

        return h
        
    def c_as_matrix(self):
        return np.reshape(self.c, (self.nnuc, self.nel))
        
    def idx(self,a,i):
        return idx(self.nel,a,i)
        
