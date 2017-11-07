import unittest
import numpy as np
import math
import matplotlib.pyplot as plt
from naewdy2.ehren import *

class TestEhren(unittest.TestCase):    
    def test_first(self):
        dt = 0.1
        m = 1.0
        x0 = 0.0
        p0 = 1.0
        nt = 10
        c0 = np.ones(1)
        
        ehren = Ehren1d(dt, m, x0, p0, c0)
        hIJ = np.zeros((1,1))
        dhIJ = np.zeros((1,1))
        xIJ  = np.zeros((1,1))

        ts = []
        xs = []
        for it in range(nt):
            t = it*dt
            x = ehren.q
            print t,x

            ts.append(t)
            xs.append(x)
            
            hIJ[0,0] = x*x
            dhIJ[0,0] = 2*x
            ehren.update(hIJ,dhIJ,xIJ)

        plt.plot(ts, xs)
        plt.show()
#        self.assertAlmostEqual((yp-ym)/(2*dx), y1, 6)

if __name__=='__main__':
    unittest.main()
