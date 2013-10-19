import matplotlib.pyplot as plt
import numpy as np

from hw3_1 import fa

def ratePlot(T=0.3, b=0.5):
    r = lambda u: .5 + .5 * np.tanh((u - T)/b)
    
    f, a = fa()     
    
if __name__ == "__main__":
    ratePlot()

    plt.show()
    
