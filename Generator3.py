from pylon import Generator
from numpy import cos, sin, pi, array, zeros, ones, conj, array,arccos

# General constant
a=cos(2*pi/3)+1j*sin(2*pi/3)
Apos=array([1.0, a*a, a])# creating positive sequence three 120-shifted unity voltage.

class Generator3(Generator):
    """Three phase Generator Model
    """
    def __init__(self, bus3, type='PQ',S=array([complex(0,0),complex(0,0),complex(0,0)])):
        #Generator.__init__(self,bus)
        self.E=ones(3)*Apos
        self.Esch=ones(3)*Apos
        self.S=S # output power of generator
        self.Srated=1.0 # 3 phase rated MVA
        self.I=conj(self.S/self.E)
        self.terminalI=zeros(3)
        self.gena=Generator(bus3.Busa)
        self.genb=Generator(bus3.Busb)
        self.genc=Generator(bus3.Busc)
        self.type=type # this shows how the generator is modeled and controled.
        self.bus=bus3
        self.bus.add_generator(self)
        self.Qmax=(self.Srated*sin(arccos(0.99))) # 3-phase reactive power upper limit
        self.Qmin=-(self.Srated*sin(arccos(0.99))) # 3-phase reactive power lower limit

    def calculateTerminalCurrent(self):
        """
        Calculate injected terminal current of three phase generator.
        """
        self.Sphase=ones(3)*self.S/3 # assume that the specified power is divided
                                     # evenly among three phases.
        self.terminalI=conj(self.Sphase/self.bus.E)
        return self.terminalI
    def updateS(self,S):
        self.S=S
        self.bus.sourceS=S
