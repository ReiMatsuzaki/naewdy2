import unittest
import time
import os

import pandas as pd

from naewdy2.nsh import *
from naewdy2.math import ijv2mat

class TestCase(unittest.TestCase):
    def assertMatProp(self,mattype,a,msg=""):
        (n,m) = a.shape
        mattype_list = [mattype]
        if(mattype=="overlap"):
            mattype_list.append("hermite")

        if "overlap" in mattype_list:
            self.assertEqual(n,m,msg=msg+"\noverlap must be squared")
            for i in range(n):
                self.assertAlmostEqual(1.0, a[i,i],
                                       msg="""
{2}
diagonal element need to unity
i      : {0}
a[i,i] : {1}
""".format(i, a[i,i], msg))

        if "hermite" in mattype_list:
            for i in range(n):
                for j in range(i):
                    self.assertAlmostEqual(a[i,j], a[j,i], msg="""
{4}
hermiticity is broken
i  : {0}
j  : {1}
a[i,j] : {2}
a[j,i] : {3}
""".format(i,j,a[i,j],a[j,i],msg))
        
    def assertMatEqual(self,a,b,prec=7,msg=""):
        self.assertEqual(a.shape,b.shape, msg="""{2}                         
shape of a and b are different
a: {0}
b: {1}
        """.format(a.shape, b.shape, msg))
        
        (n,m) = a.shape
        for i in range(n):
            for j in range(m):
                self.assertAlmostEqual(a[i,j], b[i,j], prec, msg="""{4}
Matrix a and b is largely different at (i,j).
(i,j) = ({0},{1})
a[i,j] = {2}
b[i,j] = {3}
""".format(i,j,a[i,j],b[i,j],msg))

class TestNsh(TestCase):
    def _test_read(self):
        with open("out/0/out/nshel.json") as f:
            j = json.load(f)
            gs = nshel2gto(j)
            print gtomat(gs, prim_dw(2))
                        
    def _test_first(self):
        g1 = GTO([1.0], [1.0], [0.0,0.0,0.0], [0,0,0], True)
        g2 = GTO([1.3], [1.0], [0.0,0.0,0.0], [0,0,0], True)
        gs = [g1, g2]
        s = gtomat(gs, op_s())

    def _test_s_d(self):
        """
        g1 = GTO([18.73,2.83,0.64], [0.215,0.36,0.42],
        [0.0,0.0,0.0], [0,0,0], True)
        g2 = GTO([0.2], [1.0],
        [0.0,0.0,4.0], [2,0,0], True)
        print g2
        """
        g1 = GTO([1.1], [1.0],
                 [0.0,0.0,0.0], [1,0,0], True)
        g2 = GTO([0.2], [1.0],
                 [0.0,0.0,4.0], [1,0,1], True)
        #print g1
        #print g2
        #print gtomat([g1,g2], prim_s)
        #print gtoele(g1, prim_s, g2)

    def _test_dz(self):
        
        g1 = GTO([1.5, 0.3], [1.0, 0.8], [0.0,0.0,0.0], [2,0,0], True)
    
        ex = [1.3, 0.5]
        cs = [0.8, 0.4]
        r  = np.array([0.0, 0.0, 1.1])
        dx = 0.001
        dr = np.array([0.0, 0.0, dx])
        n  = [0,0,2]
        g0 = GTO(ex, cs, r,    n, True)
        gp = GTO(ex, cs, r+dr, n, True)
        gm = GTO(ex, cs, r-dr, n, True)

        op = op_s()
        ref = (gtoele(g1,op_s(),gp) - gtoele(g1,op_s(),gm))/(2*dx)
        calc = gtoele(g1, op_dw(2), g0)
        self.assertAlmostEqual(ref, calc, 5)

    def _test_gammainc(self):
        from scipy import integrate
        import numpy as np

        t = np.linspace(0,1.0,220)
        for (m,z) in [(0,1.0), (0,2.0), (1,1.1), (2,1.2), (4,1.1), (1,0.0)]:
            y = t**(2*m) * np.exp(-z*t*t)
            ref = integrate.simps(y, t)
            calc = mole_gammainc(m, z)
            self.assertAlmostEqual(ref, calc,
                                   msg="""
ref = {0}
calc = {1}
|ref-calc| = {4}
(m,z) = ({2},{3})""".format(ref,calc,m,z,abs(ref-calc)))
        
    def _test_na(self):
        r0 = [0.0,0.0,0.0]
        g0 = GTO([1.2],      [1.0],      r0, [0,0,0], True)
        g1 = GTO([1.5, 0.3], [1.0, 0.8], r0, [2,0,0], True)

        c = np.array([1.0,0.0,0.0])
        #print "na:"
        #print gtoele(g1, prim_na(c), g0)
            
class TestMatEle(TestCase):
    def setUp(self):
        with open("gms/out/nshel.json") as f:
            j = json.load(f)
            self.gs = nshel2gto(j, True)
            self.c = j["c"]
            self.zan = j["zan"]
        
    def _test_smat(self):
        
        calc = gtomat(self.gs, op_s())
        df = pd.read_csv("gms/out/s.csv")
        ref  = ijv2mat(df)
        
        self.assertMatEqual(ref, calc)
            
    def _test_tmat(self):
        calc = gtomat(self.gs, op_t())
        df = pd.read_csv("gms/out/t.csv")
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc)
        #self.assertAlmostEqual(0.0, np.max(abs(calc-ref)))
        
    def _test_hmat(self):
        
        print "basis size:", len(self.gs)
        
        t0 = time.time()        
        calc = gtomat(self.gs, op_t())
        t1 = time.time()
        for ia in range(3):
            z = self.zan[ia]
            wc = [self.c[ir][ia] for ir in range(3)]
            calc += z*gtomat(self.gs, op_na(wc))
        t2 = time.time()        
        df = pd.read_csv("gms/out/h.csv")
        ref = ijv2mat(df)
        self.assertMatEqual(ref,calc)
        t3 = time.time()

        print "calc tmat:", t1-t0
        print "calc vmat:", t2-t1
        print "check mat:", t3-t2    

    def test_nshel(self):

        nucs = Nucs()
        ia1 = nucs.add_atom([0.0,0.0,0.0], 1, 1.0)
        ia2 = nucs.add_atom([1.0,0.0,0.0], 2, 0.3)
        nshel = Nshel(nucs)
        nshel.add_shel("s", [1.1],   {0:[1.0]}, ia1)
        nshel.add_shel("p", [1.3],   {1:[1.0]}, ia2)
        nshel.add_shel("dxx", [1.2], {2: [1.0]}, ia2)
        
        nshel.setup(True)
        gs = nshel.to_gtos()
        
        calc = nshel.smat()
        ref = gtomat(gs, op_s())
        self.assertMatProp("overlap", calc)
        self.assertMatEqual(ref, calc)

        calc = nshel.tmat()
        ref = gtomat(gs, op_t())
        self.assertMatProp("hermite", calc)
        self.assertMatEqual(ref, calc)

        calc = nshel.vmat()
        ref = gtomat(gs, op_na(nucs.ws[0])) + 0.3*gtomat(gs, op_na(nucs.ws[1]))
        self.assertMatProp("hermite", calc)
        self.assertMatEqual(ref, calc, msg="test_nshel.Check Nuclear Attraction")   
        
    def test_nshel_h2(self):
        out = "gms/h2/out"
        with open(os.path.join(out, "nshel.json")) as f:
            j = json.load(f)
            nshel = nshel_load(j)
            nshel.setup()
            
        calc = nshel.smat()
        df = pd.read_csv(os.path.join(out, "s.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_nshel_h2")

        calc = nshel.tmat()
        df = pd.read_csv(os.path.join(out, "t.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_nshel_h2. T matrix")

        calc = nshel.tmat() + nshel.vmat()
        df = pd.read_csv(os.path.join(out, "h.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_nshel_h2. H core matrix")        
        
    def test_nshel_gms(self):

        with open("gms/out/nshel.json") as f:
            j = json.load(f)
            nshel = nshel_load(j)
            nshel.setup(True)
            gs = nshel2gto(j, True)

        calc = nshel.smat()
        self.assertMatProp("overlap", calc)
        df = pd.read_csv("gms/out/s.csv")
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc)

        calc = nshel.tmat()
        self.assertMatProp("hermite", calc)
        df = pd.read_csv("gms/out/t.csv")
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc)

        calc = nshel.tmat() + nshel.vmat()
        self.assertMatProp("hermite", calc)
        df = pd.read_csv("gms/out/h.csv")
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_neshl_gms. H core")        

        
if __name__ == '__main__':
    unittest.main()
        
