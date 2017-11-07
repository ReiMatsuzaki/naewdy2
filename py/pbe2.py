import numpy as np
from naewdy2.math import uni_inte

tr = np.transpose
dot = np.dot

class PBE1d1p(object):
    def __init__(self, dt, m, dvr, c0, uni_inte_opts={}):
        self.dt = dt
        self.dvr = dvr
        (self.nnuc, self.nel) = c0.shape
        if(dvr.num != self.nnuc):
            raise RuntimeError("dvr size and given c0 is mismatch")        

        self.c = np.reshape(c0, (self.nel*self.nnuc))
        n = self.nel*self.nnuc
        self.n = np.reshape(c0, n)
        norm = np.sqrt(dot(self.c.conj(), self.c))
        self.c /= norm

        self.uni_inte_opts = uni_inte_opts

        self.tmat = -1/(2*m)*self.dvr.dmat(2)

    def update(self, dRn, hIJn, xIJn, vnn_An, vne_An):
        """
        dRn :
        .      dot{R}_n. time derivative of classical path Rn(t)
        hIJn : 
        .      electronic Hamiltonian at claasial path Rn(t)
        xIJn : 
        .      Derivative coupling evaluated at classical path Rn(t)
        vnn_An :
        .     {V_nn(RA) - V_nn(Rn)}_A
        vne_An :
        .     {<Phi_I(Rn)|V_ne(RA)-V_ne(Rn)|Phi_J(Rn)>}_{AIJ}
        """
        h = self.calc_hmat(dRn, hIJn, xIJn, vnn_An, vne_An)
        self.c = uni_inte(h, self.dt, self.c, self.uni_inte_opts)

    def calc_hmat(self, dRn, hIJn, xIJn, vnn_An, vne_An):
        n = self.nel*self.dvr.num
        h = np.zeros((n,n), dtype=complex)

        elid = np.identity(self.nel)
        for a in range(self.nnuc):
            a0 = self.idx(a,0)
            a1 = self.idx(a+1,0)
            h[a0:a1,a0:a1] = hIJn + vne_An[a,:,:] + elid*vnn_An[a]
            for b in range(self.nnuc):
                b0 = self.idx(b,0)
                b1 = self.idx(b+1,0)
                h[a0:a1,a0:a1] += self.tmat[a,b]*elid

        return h

    def c_as_matrix(self):
        return np.reshape(self.c, (self.nnuc, self.nel))
        
    def idx(self,a,i):
        return a*self.nel+i
        
    
