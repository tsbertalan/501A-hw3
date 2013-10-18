import numpy as np
from integrateAndFire import ifNeuron, ifNetwork

def twoNet(IDC1, IDC2, Irand, tmax=500, synInc=0.2):
    '''
    A network of two mutually-inhibiting integrate-and-fire neurons.
    Immediately performs integration upon construction. This is convenient.
    '''
    def functor_In(IDCn):
        def In():
            if Irand != 0:
                return np.random.normal(loc=IDCn, scale=Irand)
            else:
                return IDCn
        return In
    I1 = functor_In(IDC1)
    I2 = functor_In(IDC2)
    
    n1 = ifNeuron(I=I1, label=r"$n_1$", synInc=synInc)
    n2 = ifNeuron(I=I2, label=r"$n_2$", synInc=synInc)
    n1.addInhibitor(n2)
    n2.addInhibitor(n1)
    
    net = ifNetwork([n1, n2])
    net.integrate(tmax)
    return net
