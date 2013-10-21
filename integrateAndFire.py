import numpy as np

class ifNeuron(object):
    
    def __init__(self, Vleak=-65, Vreset=-65, dt=0.05, gLeak=0.02, Vthresh=-50,
                 C=0.4, t0=0, V0=-65, I=1.0, tref=5, tauSyn=10, synInc=0.2,
                 Vinhibsyn=-75, label=""):
        '''An integrate-and-fire neuron, with Euler-integration. Can be combined
        with other inhibitory neurons into a network (see ifNetwork).
        
        Optional Arguments
        ==================
        Vleak=-65 : scalar
            [mV] resting potential
        Vreset=-65 : scalar
            [mv] potential to which the neuron is reset after a spike.
        dt=0.05 : scalar
            [msec] timestep for Euler-integration
        gLeak=0.02 : scalar
            [uS] leak conductance
        Vthresh=-50 : scalar
            [mV] threshold potential for resetting
        C=0.4 : scalar
            [nF] membrane capacitance
        t0=0 : scalar
            [mS] start time for Euler integration
        V0=-65 : scalar
            [mV] initial membrane potential
        I=1.0 : scalar or callable
            [nA] tonic membrane current. A zero-argument, scalar-valued
            function can also be passed.
        tref=5 : scalar
            [msec] refractory time
        tauSyn=10 : scalar
            [msec] synaptic conductance decay time constant
        synInc=0.2 : scalar or callable
            [mS/msec] synaptic increment per spike received. A zero-argument,
            scalar-valued function can also be passed.
        Vinhibsyn=-75 : scalar
            [mV] Nernst potential for the inhibitory synaptic current
        label="" : string
            a label describing this neuron'''
        self.Vleak = Vleak
        self.Vreset = Vreset
        self.dt = dt
        self.gLeak = gLeak
        self.Vthresh = Vthresh
        self.C = C
        self.t = t0
        self.V = V0
        self.I_ = I
        self.Vhist = [V0]
        self.thist = [t0]
        self.tref = tref
        self.lastSpike = t0 - Vthresh
        self.spikeTimes = []

        self.inhibitors = []
        self.inhibWeights = np.array([])
        self.excitors = []  # unused
        self.excitWeights = np.array([])  # unused
        self.Vinhibsyn = Vinhibsyn
        
        self.Ghist = []  # There is no Ghist0, like there is a V0 and a t0,
        #  since the length of this vector depends on how many inhibitors are
        #  added after instantiation. 
        
        self.tauSyn = tauSyn
        self.synInc_ = synInc
        
        self.newInhibWeights = self.newV = None
        
        self.label = label
        
    def synInc(self):
        '''If the supplied synInc object is a function, return its return value;
        else just return it.'''
        if callable(self.synInc_):
            return self.synInc_()
        else:
            return self.synInc_
        
    def addInhibitor(self, other, g0=0.0):
        '''Add another neuron to this neuron's list of inhibitors.
        
        Arguments
        =========
        other (neuron)
            The other neuron afferent to, and inhibiting, this neuron.
            
        Optional Arguments
        ==================
        g0=0.0 (float)
            The initial synaptic conductance for this pair.'''
        self.inhibitors.append(other)
        weights = list(self.inhibWeights)
        weights.append(g0)
        self.inhibWeights = np.array(weights)
        
    def Ifunc(self):
        '''If the supplied I object is a function, return its return value;
        else just return it.'''
        if callable(self.I_):
            return self.I_()
        else:
            return self.I_ 
    
    def saveNextState(self):
        '''We'll do Euler-integration for now.
        This does everything *except* updating the neuron's voltages and its
        vector of inhibitory connection weights. That way,
        neurons in a ifNetwork can get eachothers' last-step voltages for use in
        their own calculations.'''
        V = self.V
        self.t += self.dt
        self.thist.append(self.t)
        if V > self.Vthresh:  # spike
            V = self.Vreset
            self.lastSpike = self.t
            self.spikeTimes.append(self.t)
        else:
            if self.t > self.lastSpike + self.tref:  # done refracting
                V += self.dVdt() * self.dt
        # no refractory period was mentioned for the synaptic conductances.
        if len(self.inhibitors) > 0:
            self.newInhibWeights = self.inhibWeights + self.dgdt() * self.dt
        self.newV = V

    def overwriteCurrentState(self):
        '''Actually save the temporary potentials and connection weights to the
        member variables used by dVdt and dgdt. For use after all neurons in the
        network have executed saveNextState.'''
        self._setV(self.newV)
        if len(self.inhibitors) > 0:
            G = self.newInhibWeights
            self._setG(G) 

    def step(self):
        '''Perform an Euler step. For use in single-neuron networks only
        (actually, just use the integrate() method).'''
        self.saveNextState()
        self.overwriteCurrentState()
        
    def dVdt(self):
        '''The current rate of potential change.'''
        leak = self.Ifunc() - self.gLeak * (self.V - self.Vleak)
        
        inhib = 0.0
        for i in range(len(self.inhibitors)):
            inhib -= self.inhibWeights[i] * (self.V - self.Vinhibsyn)
        # Show the relative contributions of leak and inhibitory components:
#         print "abs(leak):abs(inhib)", np.abs(leak) / np.abs(inhib)
        return (leak + inhib) / self.C
    
    def dgdt(self):
        '''The current vector of rates of synaptic conductance change.'''
        thresholds = np.array([n.Vthresh for n in self.inhibitors])
        voltages = np.array([n.V for n in self.inhibitors])
        spiking = (voltages > thresholds).astype(float)
        return -self.inhibWeights / self.tauSyn + self.synInc() * spiking 
    
    def integrate(self, tmax):
        '''Euler-integrate a single neuron forward in time. Don't use for networks.
        
        Arguments
        =========
        tmax (float)
            The time (in msec) at which point integration should stop, as an absolute
            number, not a difference from the neuron's current saved time.''' 
        while self.t < tmax:
            self.step()
    
    def _setV(self, V):
        '''Pre-allocating a numpy array for histories would be good,
        but then we'd have to keep track of our current index in that array,
        and we'd still have to extend it if we integrated past our expected
        maximum time. Still, this remains a potential source of speedup if
        any is later needed.'''
        self.V = V
        self.Vhist.append(V)
        
    def _setG(self, inhibWeights):
        '''inhibWeights is a vector of weights. Inhibitory only for now.'''
        self.inhibWeights = inhibWeights
        self.Ghist.append(inhibWeights)
        
    def _setI(self, I):
        self.I_ = I
    
    def getVHist(self):
        '''Returns a 1D numpy array of the potential history.'''
        return np.array(self.Vhist)
    
    def getTHist(self):
        '''Returns a 1D numpy array of the time history.'''
        return np.array(self.thist)
    
    def getGHist(self):
        '''Returns a 2D numpy array of the vector-of-synaptic-conductances history.'''
        return np.vstack(self.Ghist)

    
class ifNetwork(object):
    
    def __init__(self, neuronList):


        '''A simple container of several ifNeurons, to facilitate their integration.
        
        Arguments
        =========
        neuronList (list of ifNeuron)
            A list of neurons, which should already have called addInhibitor on
            eachother as appropriate.
            Note that the network's start time (before integration) will be taken
            from the current time of the first neuron in the supplied list.'''
        self.neurons = neuronList
        # Ensure that the neurons are integrated through this class rather
        #  than individually (which would be dumb):
        for n in self.neurons:
            n.step = n.integrate = None
        self.t = self.neurons[0].t  # Choice of index is arbitrary.
        
    def step(self):
        '''First has all neurons calculate their next state, then has them all 
        save that state. This is done in two steps so they can reference
        eachothers' current state as they're calculating their next states.'''
        for n in self.neurons:
            n.saveNextState()
            
        for n in self.neurons:
            n.overwriteCurrentState()
            
    def integrate(self, tmax):
        '''This class's raison d'etre. Integrate the network of neurons forward
        in time.
        
        Arguments
        =========
        tmax (float)
            The time (in msec) at which point integration should stop, as an absolute
            number, not a difference from the network's current saved time.'''
        while self.t < tmax:
            self.step()
            self.t = self.neurons[0].t


def waitingTimes(tlist):
    '''Returns an ndarray vector 1 shorter than that given, containing
    the differences between consecutive entries.'''
    if len(tlist) > 0:
        T = np.array(tlist)
        return T[1:] - T[:-1]
    else:
        return np.nan


