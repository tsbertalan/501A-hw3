import matplotlib.pyplot as plt
import numpy as np

from integrateAndFire import ifNeuron, ifNetwork, waitingTimes
from twoNet import twoNet

def fa(**kwargs):
    '''Create a figure with a single axis. kwargs are passed to figure().'''
    f = plt.figure(**kwargs)
    a = f.add_subplot(1, 1, 1)
    return f, a


def part1A(tmax=1000):
    '''
    Clearly I=0.3 is a bifurcation point for the default parameters
    below this threshold, the steady-state I value is a balance if inward
    current, I, and outward (leak) current. Or perhaps I have my directions
    confused.
    '''
    I = 0.290
    d = ifNeuron(I=I, label=r"$I=%.3f$" % I)
    I = 0.300
    a = ifNeuron(I=I, label=r"$I=%.3f$" % I)
    I = 0.301
    b = ifNeuron(I=I, label=r"$I=%.3f$" % I)
    I = 0.305
    c = ifNeuron(I=I, label=r"$I=%.3f$" % I)
    for n in d, a, b, c:
        n.integrate(tmax)
        for i in range(len(n.Vhist)):
            if n.Vhist[i] > n.Vthresh:
                n.Vhist[i] = 10
    fig, ax = showStateHistories([d, a, b, c], showg=False)
    ax.legend(loc="best")
    save(fig, '501a-hw3-part1A', enum=False, verbose=True)
        
        
def part1B():
    
    n = ifNeuron(dt=0.1)
    Ihist = [0]
    thist_I = [0]
    
    def Ioft(t):
        return t * 4.0 / 12500.0
    def rampI():
        I = Ioft(n.t)
        Ihist.append(I)
        thist_I.append(n.t)
        return I
    n.I_ = rampI
    
    fig, ax = fa()
    
    n.integrate(tmax=12500)
    spikeDelays = []
    
    Vhist = n.getVHist()
    thist = n.getTHist()
    
    Vhist[Vhist > n.Vthresh] = 10
    
    spikeDelays = waitingTimes(np.array(n.spikeTimes))
    delayTimes = np.array(n.spikeTimes[1:])
    currents = Ioft(delayTimes)
    
    ax.scatter(currents, 1 / spikeDelays * 1000, label='simulated')
    def f(I):
        tcurve = -n.C/n.gLeak * np.log(1 - n.gLeak/I * (n.Vthresh - n.Vleak))
        return 1 / (n.tref + tcurve)
    ax.plot(currents, f(currents) * 1000, label='calculated', color="red")
    
    ax.legend(loc='best')
    ax.set_ylabel(r'spiking frequency $[\mathrm{Hz}]$')
    ax.set_xlabel(r'applied current $[\mathrm{nA}]$')
    save(fig, "part1B", enum=False, verbose=True)


def part1C_ramp():
    Irands = [4, 8, 10, 0, 1.0, 2.0]
    fig = plt.figure(figsize=(11,8.5))
    A = [fig.add_subplot(2, 3, i+1) for i in range(6)]
    for i, Irand in enumerate(Irands):
        ax = A[i]
        n = ifNeuron(dt=0.1)
        Ihist = [0]
        thist_I = [0]
        
        def Ioft(t):
            return t * 4.0 / 12500.0
        def Ifunc():
            I = Ioft(n.t)
            Ihist.append(I)
            thist_I.append(n.t)
            # np.random.normal(), without args, returns a single, standard-normal, number. So,
            return I + Irand * np.random.normal()
            # ... is equivalent to np.random.normal(loc=I, scale=Irand) 
        n.I_ = Ifunc
        
        n.integrate(tmax=12500)
        spikeDelays = []
        
        Vhist = n.getVHist()
        thist = n.getTHist()
        
        Vhist[Vhist > n.Vthresh] = 10
        
        spikeDelays = np.array(n.spikeTimes[1:]) - np.array(n.spikeTimes[:-1])
        delayTimes = np.array(n.spikeTimes[1:])
        currents = Ioft(delayTimes)
        
        ax.scatter(currents, 1 / spikeDelays * 1000, label='simulated')
        
        def f(I):
            tcurve = -n.C/n.gLeak * np.log(1 - n.gLeak/I * (n.Vthresh - n.Vleak))
            return 1 / (n.tref + tcurve)
        ax.plot(currents, f(currents) * 1000, label='calculated', color="red")
        
        ax.set_title(r'$I_\mathrm{rand}=%.1f$ $[\mathrm{nA}]$' % Irand)

        ax.set_ylim((0, 200))
        ax.set_xlim((0, 4))

        if i is 3:
            ax.legend(loc='best')
            ax.set_ylabel(r'spiking frequency $[\mathrm{Hz}]$')
            ax.set_xlabel(r'applied current $[\mathrm{nA}]$')
        else:
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            
        save(fig, 'part1C_ramp', enum=False)


def part2_timeCourse(Irand=6, tmax=200, IDC2=1):
    
    # Even base-currents
    net = twoNet(1, IDC2, Irand, tmax=tmax)
    n1, n2 = net.neurons

    V1 = n1.getVHist()
    t  = n1.getTHist()
    V2 = n2.getVHist()
    
    fig, axes = showStateHistories([n1, n2])
    
    fig.suptitle(  "n1: %d spikes" % len(n1.spikeTimes) + '; '
                 + "n2: %d spikes" % len(n2.spikeTimes) + '\n'
                 + r"$I_\mathrm{rand}=%.2f$" % Irand)
    fig.subplots_adjust(top=0.9)
    save(fig, "part2_timecourse-Irand%.2f" % Irand, enum=False, verbose=True)
    

def part2_ratioData(ntrials=20, IDC1=1.0, IDC2max=12.0, dI=.1, Irand=1.0):
    '''Current-ratio / spike-ratio relationship.'''
    from time import time
    start = time()
    
    curRats_ = []
    spiRats_ = []
    delays = []
    for trial in range(ntrials):
        print
        curTime = time() - start
        delays.append(curTime)
        if len(delays) > 1:
            remain = np.mean(waitingTimes(delays)) * (ntrials - trial)
        else:
            remain = np.nan
        print "trial", trial+1, "of", ntrials, ",",
        print "%.1f min remaining" % (remain / 60.),
        curRats = []
        spiRats = []
        for IDC2 in np.arange(IDC1, IDC2max+dI, dI):
            net = twoNet(IDC1, IDC2, Irand, tmax=200)
            n1, n2 = net.neurons
            n1spikes = len(n1.spikeTimes)
            n2spikes = len(n2.spikeTimes)
            curRats.append(IDC2 / IDC1)
            if abs(IDC2 % 1.0) < dI:
                print int(IDC2),# ":", n1spikes, n2spikes
            spiRats.append(n2spikes / float(n1spikes))
        curRats_.append(np.array(curRats))
        spiRats_.append(np.array(spiRats))
    
    curRats = np.vstack(curRats_)
    spiRats = np.vstack(spiRats_)
    fname = "driveRatios-tri_%d-IDC1_%.1f-IDC2max_%.1f-dI_%.1f-Irand_%.1f.npz" % (
                                             ntrials, IDC1, IDC2max, dI, Irand)
    np.savez_compressed(fname, curRats=curRats, spiRats=spiRats,
                        curRats_=curRats_, spiRats_=spiRats_,
                        ntrials=ntrials, IDC1=IDC1, IDC2max=IDC2max, dI=dI)
    print
    print 'saving', fname
    return fname
    
def part2_ratioShow(fname):
    '''Show the saved data from part2_ratioData.'''
    fig, ax = fa()
    
    data = np.load(fname)
    curRats = data['curRats']
    spiRats = data['spiRats']
    curRats_ = data['curRats_']
    spiRats_ = data['spiRats_']
    IDC1 = data['IDC1']
    IDC2max = data['IDC2max']

    for c, s in zip(curRats_, spiRats_):
        ax.plot(c, s, color=(0,0,0,.01)) # black; 1% opacity
    # average across four trials
    curRats = np.mean(curRats, axis=0)
    window_len = 32
    curRats = smooth(curRats, window_len=window_len)
    spiRats = np.mean(spiRats, axis=0)
    spiRats = smooth(spiRats, window_len=window_len)
    
    ax.set_xlim((IDC1, IDC2max))
    ax.plot(curRats, spiRats, color="black")
    ax.set_ylabel(r"$nSpikes_2:nSpikes_1$")
    ax.set_xlabel(r"$I_\mathrm{DC2}:I_\mathrm{DC1}$")
    print curRats.shape
    print curRats[::window_len].shape
#     for x in curRats[::window_len]:
#         ax.axvline(x=x)
    ax.set_xscale('log')
    save(fig, fname, verbose=True, enum=False)


def part3():
    '''Input current is not noisy, but syanptic conductances are Poisson-
    distributed.'''
    m = 40
    def synInc():
        return 0.2 * np.random.poisson(lam=m)
    net = twoNet(1, 1, 0, synInc=synInc, tmax=1000)
    n1, n2 = net.neurons
    
    # cause the lower-g neuron to have ... a higher g.
    if n1.inhibWeights[0] < n2.inhibWeights[0]:
        n1.inhibWeights[0] = 0.7
    else:
        n2.inhibWeights[0] = 0.7
    
    net.integrate(tmax=2000)
    
    fig, axes = showStateHistories([n1, n2])
    save(fig, 'part3', enum=False, verbose=True)
    
    
def showStateHistories(listOfNeurons, showg=True):
    fig = plt.figure()
    if showg:
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
    else:
        ax1 = fig.add_subplot(1, 1, 1)
    for n in listOfNeurons:
        T = n.getTHist()
        ax1.plot(T, n.getVHist(), label=n.label)
        if showg:
            ax2.plot(T[1:], n.getGHist()[:, 0], label=n.label)
    
    xlabel = r"$t$ $[\mathrm{msec}]$"

    ax1.set_ylabel(r"$V$ $[mV]$")
    ax1.legend(loc="best")
    ax1.set_xlim((T.min(), T.max()))

    if showg:
        ax1.set_xticklabels([])
        ax2.set_ylabel(r"$g_{i,j}$ $[\mathrm{\mu S}]$")
        ax2.legend(loc="best")
        ax2.set_xlim((T.min(), T.max()))
        ax2.set_xlabel(xlabel)
    else:
        ax1.set_xlabel(xlabel)
    
    if showg:
        return fig, [ax1, ax2]
    else:
        return fig, ax1


def save(fig, filename, enum=True, ext=".pdf", verbose=False):
    if enum:
        from time import time
        filename += "-%d" % time()
    filename += ext
    if verbose:
        print "saving", filename
    fig.savefig(filename)


def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=numpy.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]    
    
    
if __name__ == "__main__":
    part1A()
    part1B()
    part1C_ramp()
    part2_timeCourse(Irand=0.0, tmax=400)
    part2_timeCourse(Irand=0.8, tmax=400)
    part2_timeCourse(Irand=3.0, tmax=400)
    
#     fname = part2_ratioData(ntrials=2, IDC2max=2.0, dI=0.25)
#     fname = part2_ratioData(ntrials=44, IDC2max=24.0, dI=0.025)
#     fname = "driveRatios-tri20-IDC11.000000-IDC2max12.000000-dI0.100000.npz"
    fname = "driveRatios-tri_44-IDC1_1.0-IDC2max_24.0-dI_0.0-Irand_1.0.npz"
    part2_ratioShow(fname)
    part3()

    plt.show()
    
