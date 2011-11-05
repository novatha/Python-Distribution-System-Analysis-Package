from pylon import Case
from networkx import DiGraph, shortest_path_length,shortest_path, draw
from pylab import plot,show
from numpy import zeros, imag
from collections import defaultdict

rootNode="root node"
forkNode="fork node"
endNode="end node"
normalNode="normal node"

class Case3(Case):
    """
        This class defines a radial network which will be studied further. It uses directed graph implemented in networkx.
    """
    def __init__(self,name=None, base_mva=1.0, base_kv=20.0, buses=None, branches=None,  generators=None):

        #Case.__init__(self,buses=buses,branches=branches)
        self.buses=buses
        self.base_mva=base_mva
        self.base_kv=base_kv
        self.base_z=base_kv**2/base_mva
        self.branches=branches
        self.generators=generators
        self.G=DiGraph()
        self._buildGraph()
        self.rootbus=[b for b in self.buses if self.G.pred[b] == {}]
        self.forkbuses=[b for b in self.buses if self.G.degree()[b]>2]
        self.normalbuses=[b for b in self.buses if self.G.degree()[b]==2]
        self.endbuses=[b for b in self.buses if self.G.succ[b] == {}]
        self.updateTypeG()
        self.bus_level=shortest_path_length(self.G,source=self.rootbus[0])
        self.updateBusConnBranch()
        self.buses_loadJ=self.updateBusLoadCurrent()
        self.identifyPVbuses()
	self.calculateZpos()
	for b in self.buses:
            b.base_mva=base_mva
        #draw(self.G)
        #pl.show()
    def _buildGraph(self):
        """ add buses as nodes  and branches as edges.
        """
        self.G.add_nodes_from([bus for bus in self.buses])
        self.G.add_edges_from([(b.from_bus,b.to_bus) for b in self.branches])
    def plotGraph(self):
        draw(self.G)
        show()
        return
    def plotVP(self):
        """ Plot voltage profile of distribution system"""
        V=[abs(b.E) for b in self.buses]
        name=[b.name for b in self.buses]
        plot(V)
        show()
        return

    def updateBusConnBranch(self):
	"""
	Updating to which branch a bus is connected.
	"""
        for a in self.branches:
            a.from_bus.connected_to_branch.append(a)
            a.to_bus.connected_from_branch.append(a)
        # as edges are added, there should be a procedure to update bus.connected_to and bus.connected_from
    def updateTypeG(self):
        for b in self.buses:
            if self.rootbus.__contains__(b):
                b.typeG=rootNode
            if self.forkbuses.__contains__(b):
                b.typeG=forkNode
            if self.endbuses.__contains__(b):
                b.typeG=endNode
    def updateBusLoadCurrent(self):
       return [b.calc_loadI() for b in self.buses]

    def getRoot(self):
	"""
	Get the root node of a network which is the substation node.
	"""
        rootNode=[n for n in self.G.node if self.G.pred[n]=={}]
        return rootNode[0]

    def updateBaseKV(self,aNode):
        """
        update base_kv of aNode children. This is done after a voltage regulator changes its tap.
        """
        for b in aNode.connected_to_branch:
            b.updateToBusBaseKV()
            if b.to_bus.connected_to_branch!=[]:
                self.updateBaseKV(b.to_bus)
    def identifyPVbuses(self):
	"""
	identify pv buses in the network. Important when creating the positive sequence matrix needed to update injected current
	"""
	self.pvbuses=[g.bus for g in self.generators if g.type=='PV'] # a list of buses connected to generators controlled as PV bus.
	for b in self.pvbuses:
		b.type='PV' # update the bus.type according to its generator control mode.
    def calculateZpos(self):
	"""
	calculate positive impedance matrix for updating injected current at PV bus  in BFS with DGs

	"""
	self.identifyPVbuses()
	pathToRoot={}
	lines={}
	Z1diag={}
	pvbusindex={}
	d=defaultdict(dict)
	npvb=len(self.pvbuses) # number of pv buses
        Z1=zeros([npvb,npvb],dtype='complex')
	for pvb in self.pvbuses:
		pathToRoot[pvb]=shortest_path(self.G,self.rootbus[0],pvb)
    		lines[pvb]=[b.connected_from_branch[0] for b in pathToRoot[pvb] if b.connected_from_branch!=[]]
        self.pathToRoot=pathToRoot# a list of buses that incident with the path from a pv bus to root bus.
	self.lines=lines # a list of lines that forms the path from a pvbus to the root bus.
	Z_1=zeros([npvb,npvb]) # prepare an all zero matrix with size npvb x npvb
	# This function requires a positif sequence of all lines in self.lines. The positive sequence impedance is available at line.Z012[0]
	# the situation is more intricate if there is a transformer within the path between the PV bus to root bus. If on of the transformer winding
	# is connected in delta, the path for positive sequence is disconnected. The implementation for this kind of situation is postponed.
	# the current implementation assumes that there is NO TRANSFORMER exist in the path. 27th May 2011.
	for pvb in self.pvbuses:
		Z1sum=complex(0,0)
		for l in self.lines[pvb]:
			Z1sum=Z1sum+l.Z012[1][1]
		Z1diag[pvb]=Z1sum
	self.Z1diag=Z1diag # this is the diagonal element of positive impedance matrix, the output of this function.
	ind=0
	for pvb in self.pvbuses:
		pvbusindex[pvb]=ind
		ind=ind+1
	self.pvbusindex=pvbusindex
	for pvb in self.pvbuses:
		for pv in self.pvbuses:
			d[pvb][pv]=set(self.lines[pvb]).intersection(set(self.lines[pv]))
			Z1sum=0
			for l in d[pvb][pv]:
				Z1sum=Z1sum+l.Z012[1][1]
			Z1[pvbusindex[pvb],pvbusindex[pv]]=Z1sum
	self.d=d# this is a two dimensional dictionary. Its diagonal element contains list of lines on the path between a pv bus and root bus.
	        # Its of diagonal contains the common lines between path to root of two pv buses. If we sum the positive sequence impedance
		# of each element, we have already the element of postive sequence impedance matrix we wish to calculate in this function.
	self.Z1=Z1
	self.X1=imag(Z1)
def P_to_S(P,pf):
    theta=arccos(pf)
    Samp=P/pf
    S=Samp*cos(theta)+1j*Samp*sin(theta)
    return S

def polarToRect(magnitude,degree):
    result=magnitude*(cos(degree)+1j*sin(degree))
    return result

