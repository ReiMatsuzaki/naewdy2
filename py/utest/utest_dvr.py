import unittest
import numpy as np
import math

import matplotlib.pyplot as plt

from naewdy2.math import *
from naewdy2.dvr import *


class TestDVR(unittest.TestCase):    
    def test_at(self):
        dx = 0.001
        x  = 0.3
        dvr = ExpDVR(3, -4.0, 4.0)
        cs = [0.2, 0.3, 0.4, 0.5, 1.1, 1.2, 1.3]
        y0 = dvr.at(cs, x,    nd=0)
        yp = dvr.at(cs, x+dx, nd=0)
        ym = dvr.at(cs, x-dx, nd=0)
        y1 = dvr.at(cs, x,    nd=1)
        y2 = dvr.at(cs, x,    nd=2)

        self.assertAlmostEqual((yp-ym)/(2*dx), y1, 6)
        self.assertAlmostEqual((yp+ym-2*y0)/(dx*dx), y2, 6)
    
    def test_dmat(self):
        # == calculation ==
        x0 = -2.0
        xN = +2.0        
        dvr = ExpDVR(2, x0, xN)
        d1 = dvr.dmat(1)
        d2 = dvr.dmat(2)
        
        # == check basic property ==
        self.assertAlmostEqual(d1[0,1], -d1[1,0], 7)
        self.assertAlmostEqual(0.0,     d1[0,0])
        self.assertAlmostEqual(0.0,     d1[1,1])
        self.assertAlmostEqual(d2[0,1], d2[1,0], 7)

        # == check value with numerical integration ==
        xs = np.linspace(x0, xN, 1000)
        dx = xs[1]-xs[0]
        i = 0; j = 1
        c = np.zeros(dvr.num); c[i] = 1.0
        y0s = dvr.at(c, xs, nd=0)
        c = np.zeros(dvr.num); c[j] = 1.0
        y1s = dvr.at(c, xs, nd=1)
        y2s = dvr.at(c, xs, nd=2)
        
        self.assertAlmostEqual(np.dot(y0s, y1s)*dx, d1[i,j])
        self.assertAlmostEqual(np.dot(y0s, y2s)*dx, d2[i,j])

    def test_fit(self):
        dvr = ExpDVR(12, -3.0, 3.0)
        gauss = lambda x: 1.3*np.exp(-1.2*(x-0.1)**2)
        fit_c = dvr.fit(gauss)
        x = 0.2
        self.assertAlmostEqual(dvr.at(fit_c, x), gauss(x), 6)
    
    def test_harmo(self):

        m = 2000.0
        k = 0.5
        w = np.sqrt(k/m)

        def harmo(n, x):
            mw = sqrt(k*m)
            z = sqrt(mw)*x
            a = sqrt(1.0/(math.factorial(n)*2**n) * sqrt(mw/pi))
            if(n==0):
                h = 1.0
            elif(n==1):
                h = 2*z
            elif(n==2):
                h = 4*z**2-2
            return a * h * exp(-z*z/2)
        
        dvr = ExpDVR(30, -3.0, 3.0)
        d2 = dvr.dmat(2)
        v  = np.diag([k/2*x*x for x in dvr.xs])
        h = -1/(2*m)*d2 + v
        (e,u) = eigh_sort(h)

        # -- check orthogonal --
        self.assertAlmostEqual(1.0, dot(u[:,0], u[:,0]))
        self.assertAlmostEqual(0.0, dot(u[:,0], u[:,1]))

        for n in range(3):
            # -- check eigan values --
            self.assertAlmostEqual(w*(0.5+n), e[n])            

            # -- check eigen functions --
            x = 0.2
            y0 = harmo(n,x)
            y1 = dvr.at(u[:,n], x, nd=0)
            if(y0*y1.real<0):
                y1 = -y1
            self.assertAlmostEqual(y0, y1, 5, 'n={0}'.format(n))
            
    def test_krylov(self):

        m = 1.0
        w = 1.0
        a = m*w/2
        x0 = 1.0
        p0 = 0.0
        dvr = ExpDVR(128, -3.0, 3.0)
        n = dvr.num
        dt = 1.0

        c0 = np.zeros(n, dtype=complex)
        for a in range(len(dvr.xs)):
            c0[a] = gauss(a=1.0, zeta=a, r0=x0, p0=p0, x=dvr.xs[a])
            #c0 = np.ones(n, dtype=complex)
        c0 /= norm(c0)

        h = np.zeros((n,n))
        h[:,:] = np.diag([m*w*w/2*x*x for x in dvr.xs]) -1/(2*m)*dvr.dmat(2)

        c_diag = uni_inte(h, dt, c0, opts={"inte":"diag"})
        c_kry10  = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":10})
        c_kry20  = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":20})
        #c_kry30  = uni_inte(h, dt, c0, krylov_num=30, inte="krylov")
        c_kry_full  = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":n})

        print np.sum(abs(h-np.transpose(h.conj())))
        self.assertAlmostEqual(1.0, norm(c_diag))
        self.assertAlmostEqual(1.0, norm(c_kry10))

        print np.sum(abs(c_diag-c_kry10))/len(c_diag)
        print np.sum(abs(c_diag-c_kry20))/len(c_diag)
        #print np.sum(abs(c_diag-c_kry30))/len(c_diag)
        print np.sum(abs(c_diag-c_kry_full))/len(c_diag)
        
    def test_sigma_H(self):
        m = 1.2
        w = 1.0
        a = m*w/2
        x0 = 1.0
        p0 = 0.0
        dvr = ExpDVR(3, -3.0, 3.0)
        n = dvr.num
        dt = 1.0
        nel = 2

        D1 = dvr.dmat(1)
        D2 = dvr.dmat(2)

        c0 = np.zeros(n, dtype=complex)
        for a in range(len(dvr.xs)):
            c0[a] = gauss(a=1.0, zeta=a, r0=x0, p0=p0, x=dvr.xs[a])
        c0 /= norm(c0)

        hel = np.zeros((n,nel,nel))
        hel[:,0,0] = [m*w*w/2*x*x     for x in dvr.xs]
        hel[:,1,1] = [m*w*w/2*x*x+1.0 for x in dvr.xs]
        hel[:,0,1] = [0.1*x           for x in dvr.xs]
        hel[:,1,0] = hel[:,0,1]

        xij = np.zeros((n,nel,nel))
        xij[:,1,0] = [0.01*np.exp(-(x-1.0)**2) for x in dvr.xs]
        xij[:,0,1] = -xij[:,1,0]

        #c = np.array([(jn+ie)/1000.0 for jn in range(n) for ie in range(nel)])
        c = np.array([(jn+ie)/1000.0 for ie in range(nel) for jn in range(n) ])
        Hc_sigma = sigma_H(hel, xij, D1, D2, m, c)
        Hc_build = np.dot(build_H(hel, xij, D1, D2, m), c)

        self.assertAlmostEqual(0.0, np.sum(abs(Hc_sigma-Hc_build))/len(c))
            
if __name__=='__main__':
    unittest.main()
