from pylon.util import _Named
from numpy import array,conj,zeros
from numpy import dot
from Bus3 import Bus3
from constants import A

D=array([[1.0,-1.0,0.0],[0.0,1.0,-1.0],[-1.0,0.0,1.0]])
DI=D.transpose().copy()
class noConnectionSpecified(Exception):
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg


class Load3(_Named):
	
    def __init__(self,Sa=0.0,Sb=0.0,Sc=0.0,connection='delta',type='CS', bus=Bus3()):
        """
        Load3 is a CONSTANT-POWER three-phase load model.

        """
        self.S=array([Sa,Sb,Sc]) # complex power in kVA
        self.connection=connection
        self.type=type
        self.bus=bus
        self.E=self.bus.E
        self.Ell=dot(D,self.E)
        self.Iterm=0
        self.calculateLineCurrent()
    def calculateLineCurrent(self):
        """
        Calculate the terminal current of this load.

        """
        self.E=self.bus.E
        self.Ell=dot(D,self.E)
        if self.connection=='delta':
            self.Ill=conj(self.S/self.Ell) # delta current
            self.Iterm= dot(DI,self.Ill) # delta current to line current
        elif self.connection=='wye':
            self.Iln=conj(self.S/self.E)
            self.Iterm=self.Iln
        else:
            raise noConnectionSpecified
        
class Load3CZ(Load3):
    def __init__(self,Z,connection='delta',bus=Bus3()):
        """
        Load3 is a CONSTANT-IMPEDANCE three-phase load model based on Load3.

        """
        self.Za=Z[0]
        self.Zb=Z[1]
        self.Zc=Z[2]
        self.Z=array([self.Za, self.Zb, self.Zc])
        self.connection=connection
        self.bus=bus
        self.E=bus.E
    def calculateLineCurrent(self):
        """
        Calculate the terminal current of this load.

        """
        self.E=self.bus.E
        self.Ell=dot(D,self.E)
        if self.connection=='delta':
            self.Ill=self.Ell/self.Z # delta current
            self.Iterm= dot(DI,self.Ill) # delta current to line current
        elif self.connection=='wye':
            self.Iln=self.E/self.Z
            self.Iterm=self.Iln
        else:
            raise noConnectionSpecified
        
class Load3CI(Load3):
    def __init__(self,I=zeros(3),connection='delta',bus=Bus3()):
        """
        Load3 is a CONSTANT-CURRENT three-phase load model based on Load3.

        """
        self.Ia=I[0]
        self.Ib=I[1]
        self.Ic=I[2]
        self.I=array([self.Ia, self.Ib, self.Ic])
        self.connection=connection
        self.bus=bus
        self.E=bus.E
    def calculateLineCurrent(self):
        """
        Calculate the terminal current of this load.

        """
        self.E=self.bus.E
        self.Ell=dot(D,self.E)
        if self.connection=='delta':
            self.Ill=self.I # delta current
            self.Iterm= dot(DI,self.Ill) # delta current to line current
        elif self.connection=='wye':
            self.Iln=self.I
            self.Iterm=self.Iln
        else:
            raise noConnectionSpecified
        
        
