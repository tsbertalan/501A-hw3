import matplotlib.pyplot as plt
import numpy as np

from hw3p1 import fa, save

def ratePlot(T=0.3, b=0.5):
    r = lambda u: .5 + .5 * np.tanh((u - T)/b)
    f, a = fa()
    
    u = np.arange(T-b*5, T+b*5, .01)
    
    a.plot(u, r(u), 'k')
    a.set_xticks([T-b, T, T+b])
    a.set_xticklabels([r'$\theta-\beta$', r'$\theta$', r'$\theta+\beta$'])
    a.set_xlim([-b*5, b*5])
    
    a.axvline(T, color='black')
    a.axvline(T-b, color='black')
    a.axvline(T+b, color='black')
    a.axhline(.5, color='black')
    
    a.set_xlabel(r'$u(t)$')
    a.set_ylabel(r'$r(u(t))$')
    
    save(f, 'r_of_u')
    
if __name__ == "__main__":
    ratePlot()

    plt.show()
    
