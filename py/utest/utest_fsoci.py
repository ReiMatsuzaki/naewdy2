import unittest
from numpy import sqrt
from naewdy2.fsoci import *

class TestFsoci(unittest.TestCase):
    def test_sign_ai(self):
        self.assertAlmostEqual(0, sign_ai([1,2,3], 4))
        self.assertAlmostEqual(1, sign_ai([1,2,3], 3))
        self.assertAlmostEqual(-1, sign_ai([1,2,3], 2))
        self.assertAlmostEqual(1, sign_ai([1,2,3], 1))

    def test_aiaj(self):
        self.assertAlmostEqual(1,  aiaj([1,2,3], 1, 1, [1,2,3]))
        self.assertAlmostEqual(0,  aiaj([1,2,3], 4, 1, [1,2,3]))
        self.assertAlmostEqual(1,  aiaj([1,2,4], 4, 3, [1,2,3]))
        self.assertAlmostEqual(-1, aiaj([1,3,4], 4, 2, [1,2,3]))

    def test_eij(self):
        self.assertAlmostEqual(sqrt(2.0),
                               eij([1,2,3], [1,2,3],
                                   1, 1,
                                   [1,2,3], [1,2,3]))
        
if __name__ == '__main__':
    unittest.main()
