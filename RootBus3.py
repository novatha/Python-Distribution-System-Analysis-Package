from Bus3 import Bus3
from numpy import zeros

class RootBus3(Bus3):
    def __init__(self,busa=None,busb=None,busc=None,base_kv=1):
        Bus3.__init__(self,busa=None,busb=None,busc=None,base_kv=base_kv)
	self.type='PV'
	self.sourceS=[0j,0j,0j]
	self.loadS=zeros(3,dtype=complex)
    def whatIam(self):
        return "RootBus3"
