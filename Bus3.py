from pylon import Bus
from numpy import cos,sin, pi, array, ones, zeros, conj, sum


# General constant
a=cos(2*pi/3)+1j*sin(2*pi/3)
A=array([1.0,a,a*a])# creating three 120-shifted unity voltage.
P=array([[1,1,1],[1,a,a*a],[1,a*a,a]])/3
Pinv=array([[1,1,1],[1,a*a,a],[1,a,a*a]])
ps=cos(pi/6)+1j*sin(pi/6) # this is a 30 degree phase shift
PS=array([[ps,0,0],[0,ps,0],[0,0,ps]])
rootNode="root node"
forkNode="fork node"
endNode="end node"
normalNode="normal node"
satu=ones(3)*A
class Bus3(Bus):
    """Three phase bus model
       I think that every bus should have self.kvbase attribut.
    """
    def __init__(self,voltage=satu,busa=None,busb=None,busc=None,base_kv=1.0,base_mva=1.0,ps=0,loads=[]):
        Bus.__init__(self)
        self.Busa=busa
        self.Busb=busb
        self.Busc=busc
        self.base_kv=base_kv*ones(3)
        self.base_mva=base_mva
        self.typeG=normalNode
        self.loadS=zeros(3,dtype=complex)   #
        self.loads=loads
        self.gens=[]
        self.sourceS=zeros(3) # please assume generation as negative load
        self.shuntS=zeros(3)
        self.E=satu*A# balanced unity voltages. E is line to neutral voltage
        self.connected_from_branch=[] # self is connected to root from branch
        self.connected_to_branch=[] # self is connected to branch
        self.loadI=zeros(3)
        self.sourceI=zeros(3)
        self.childI=zeros(3)
        self.shuntI=zeros(3)
        self.totalI=zeros(3)
        self.Iqpv=zeros(3) # the reactive current injected by a generator in pv-node
    def __repr__(self):
        """ A function that modifies "print this_bus" result
        """
        return self.name

    def set_load(self,S):
        #a,b,c=asarray(S)
        self.loadS=asarray(S)
    def get_load(self):

        return self.loadS

    def calc_loadI(self):
        """A function that calculate current of load directly connected
           to this bus
        """
        L=[load for load in self.loads]
        acc=0
        for load in L:
            load.calculateLineCurrent()
            acc= acc+load.Iterm # this is for wye conected load
        self.loadI=acc
        return self.loadI
    def calc_sourceS(self):
        """A function that calculate power generated to this bus
        """
        g=[gen for gen in self.gens]
        acc=0
        for gen in g:
            acc= acc+gen.S # Negativity is considered in self.calc_sourceI
        self.sourceS=acc
        return self.sourceS
    def calc_sourceI(self):
        """ A function that calculate current from generator connected to this
            bus.
        """
        g=[gen for gen in self.gens]
        acc=0
        for gen in g:
            acc= acc+gen.calculateTerminalCurrent() # Negativity is considered in self.calc_sourceI
        self.sourceI=-acc # its negative because it represents current injected to the bus. We assume load current to be positive.
        return self.sourceI
    def calc_shuntI(self):
        self.shuntI=conj(self.shuntS/self.E)
        return self.shuntI
    def calc_childI(self):
        """ A function that calculate currents flowing from this bus to all
            other child buses.
        """
        self.childI=sum([b.I_from for b in self.connected_to_branch ],axis=0)
        return self.childI
    def calc_totalI(self):
        self.totalI=self.calc_loadI()+self.calc_sourceI()+self.calc_shuntI()+self.calc_childI()#-self.Iqpv
        return self.totalI
    def addConnectedToBranch():
        pass
    def whatIam(self):
        return "Bus3"

    def add_generator(self,g):
	if not (g in self.gens):
	    self.gens.append(g)
	self.calc_sourceS()
	self.calc_sourceI()
