import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import kevrekidis as kv


class ifNet(obj):
    def __init__(self, X0, ):
        self.X = X0
        self.t = t0
        self.Xhist = [X0]
    
    def diff(self, X, Vleak= -65, Vreset= -65, dt=0.05, gLeak=0.02,
            Vthresh= -50, C=0.4, t0=0, V0= -65, I=1.0, tref=5,
            Vi=-75, tauSyn=10, synInc=0.2):
        V = X[0:2]
        g = X[2:4]  # For a real network, this should be a matrix, and there
                     #  should be a sum down there.
        I = np.array([1, 1])
        delta = (V > Vthresh).astype(float)
        
        dV = (I - gLeak*(V - Vleak) - g*(V - Vi)) / C
        dg = -g/tauSyn + synInc * delta
        
        dX = np.array(X.shape)
        return np.hstack((dV, dg))
    
    
    def step(self):
        


def part2():
    X0 = np.array([-65, -65, .1, .1])
    X, T = kv.dynSys.integration.integrate(X0, diff)

    
if __name__ == "__main__":
    part2()
    plt.show()
    
