import unittest
import numpy as np
import math
from datetime import datetime

import matplotlib.pyplot as plt

from naewdy2.math import *
from naewdy2.dvr import *


class TestDVR(unittest.TestCase):    
    def _test_at(self):
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
    
    def _test_dmat(self):
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

    def _test_fit(self):
        print "fit begin"
        dvr = ExpDVR(12, -3.0, 3.0)
        g = lambda x:1.3*np.exp(-1.2*(x-0.1)**2)
        gs = g(dvr.xs)
        c = dvr.fit(gs)
        x = 0.2
        self.assertAlmostEqual(dvr.at(c, x), g(x), 6)
        print "fit end"
    
    def _test_harmo(self):

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
        print np.min(e)

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
        
                        
    def _test_krylov(self):

        m = 1.0
        w = 1.0
        a = m*w/2
        x0 = 1.0
        p0 = 0.0
        dvr = ExpDVR(128, -5.0, 5.0)
        n = dvr.num
        dt = 1.0

        g0 = np.exp(-a*(dvr.xs-x0)**2 + 1j*p0*(dvr.xs-x0))
        c0 = dvr.fit(g0)
        c0 /= norm(c0)

        h = np.zeros((n,n))
        h[:,:] = np.diag([m*w*w/2*x*x for x in dvr.xs]) -1/(2*m)*dvr.dmat(2)

        pd.DataFrame({"val":dvr.xs}).to_csv("out/xs.csv", index=None)
        y0 = dvr.at(c0)
        pd.DataFrame({"re":y0.real, "im":y0.imag}).to_csv("out/psi0.csv", index=None)

        (c_diag, od) = uni_inte(h, dt, c0, opts={"inte":"diag"})
        (c_kry10, o10) = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":10})
        (c_kry20, o20) = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":20})
        (c_kry60, o60) = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":60})
        (c_kry_full, of)  = uni_inte(h, dt, c0, opts={"inte":"krylov", "krylov_num":n})

        self.assertAlmostEqual(1.0, norm(c_diag))
        self.assertAlmostEqual(1.0, norm(c_kry10))

        print 
        print np.sum(abs(c_diag-c_kry10))/len(c_diag)
        print np.sum(abs(c_diag-c_kry20))/len(c_diag)
        print np.sum(abs(c_diag-c_kry60))/len(c_diag)
        print np.sum(abs(c_diag-c_kry_full))/len(c_diag)

    def test_sigma_H(self):
        m = 1.2
        w = 1.0
        a = m*w/2
        x0 = 1.0
        p0 = 0.0
        dvr = ExpDVR(256, -5.0, 5.0)
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

        d0 = datetime.now()
        Hc_sigma = sigma_H(hel, xij, D1, D2, m)(c)
        d1 = datetime.now()
        print "sigma", d1-d0

        d0 = datetime.now()
        H = build_H(hel, xij, D1, D2, m)
        d1 = datetime.now()
        print "build", d1-d0

        d0 = datetime.now()
        Hc_build = np.dot(H, c)
        d1 = datetime.now()
        print "dot", d1-d0

        self.assertAlmostEqual(0.0, np.sum(abs(Hc_sigma-Hc_build))/len(c))

    def _test_krylov_2state(self):
        m = 1.2
        w = 1.0
        a = m*w/2
        x0 = 1.0
        p0 = 0.0
        
        dvr = ExpDVR(256, -5.0, 5.0)
        n = dvr.num
        dt = 1.0
        nstate = 2

        D1 = dvr.dmat(1)
        D2 = dvr.dmat(2)        

        hel = np.zeros((n,nstate,nstate))
        hel[:,0,0] = [m*w*w/2*x*x     for x in dvr.xs]
        hel[:,1,1] = [m*w*w/2*x*x+1.0 for x in dvr.xs]
        hel[:,0,1] = [0.1*x           for x in dvr.xs]
        hel[:,1,0] = hel[:,0,1]

        xij = np.zeros((n,nstate,nstate))
        xij[:,1,0] = [0.01*np.exp(-(x-1.0)**2) for x in dvr.xs]
        xij[:,0,1] = -xij[:,1,0]

        g0 = np.exp(-a*(dvr.xs-x0)**2 + 1j*p0*(dvr.xs-x0))
        c0 = np.reshape(np.transpose(np.array([dvr.fit(g0), np.zeros(n)])), n*nstate)

        
        Hc = sigma_H(hel, xij, D1, D2, m)
        print "build begin", datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        H  = build_H(hel, xij, D1, D2, m)
        print "build end", datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        print
        print "krylov begin", datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        (c1, o1) = uni_inte_krylov(Hc, dt, 10, c0)
        print "krylov end", datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        print
        print "diag begin", datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        (c2, o2) = uni_inte(H, dt, c0, opts={"inte":"diag"})
        print "diag end", datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        print
        print "krylov 2state"
        print np.sum(abs(c1-c2))/len(c1)
        
if __name__=='__main__':
    unittest.main()
    
        