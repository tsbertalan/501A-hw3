import matplotlib.pyplot as plt
import numpy as np

from hw3p1 import fa, save

def ratePlot(T=0.3, b=0.5):
    r = lambda u: .5 + .5 * np.tanh((u - T)/b)
    f, a = fa()
    
    u = np.arange(-b*5, b*5, .01)
    
    a.plot(u, r(u), 'k')
    a.set_xticks([-b*4, b, b*4])
    a.set_xticklabels([r'$-4\beta$', r'$\beta$', r'$4\beta$'])
    a.set_xlim([-b*5, b*5])
    
    a.axvline(b, color="black")
    
    a.set_xlabel(r'$u(t)$')
    a.set_ylabel(r'$r(u(t))$')
    
    save(f, 'r_of_u')
    
if __name__ == "__main__":
    ratePlot()

    plt.show()
    
