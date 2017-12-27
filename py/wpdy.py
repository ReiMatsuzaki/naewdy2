import numpy as np

class WavePacketDyanmics(object):
    def __init__(nx):
        self.xs = np.zeros(nx)
        self.vs = np.zeros(nx)
        self.cvs = np.zeros(nx)
        self.fs = np.zeros(nx, dtype=complex)
