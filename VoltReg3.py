from pylon import Bus, Generator, Branch, Case
#from networkx import DiGraph, shortest_path_length
from scipy.sparse import hstack, vstack, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, splu
from numpy import cos, array, sin, ones, zeros, pi, arccos, dot, eye, sqrt, conj, sum, rank, append, diag
from numpy.linalg import solve, inv
from Line3 import Line3
from Bus3 import Bus3
from VoltReg1 import VoltReg1

class VoltReg3(Line3):
    """
    
    """
    def __init__(self,from_bus=Bus3(),to_bus=Bus3(),voltregs=[],yt=0,base_mva= 1.0,base_kv=12.0,case_base_mva=1.0,case_base_kv=12.0,E_high= 12.47, E_low=4.16, type_connection='YgYg',type='Step-Down'):
        Line3.__init__(self,from_bus,to_bus)
        #print voltregs
#        self.Z=zeros([3,3])
        self.voltrega=voltregs[0]
        self.voltregb=voltregs[1]
        self.voltregc=voltregs[2]
        self.voltregs=[self.voltrega, self.voltregb, self.voltregc]
        self.from_bus=from_bus
        self.to_bus=to_bus
        self.line_a=eye(3)
        self.line_b=zeros([3,3])
        self.line_c=zeros([3,3])
        self.line_d=zeros([3,3])
        self.line_A=inv(self.line_a)
        self.line_B=dot(inv(self.line_a),self.line_b)
        # attaching voltregs to buses
        ## from_bus
        self.voltrega.from_bus=self.from_bus.Busa
        self.voltregb.from_bus=self.from_bus.Busb
        self.voltregc.from_bus=self.from_bus.Busc
        ## to_bus
        self.voltrega.to_bus=self.to_bus.Busa
        self.voltregb.to_bus=self.to_bus.Busb
        self.voltregc.to_bus=self.to_bus.Busc
        self.a_R=eye(3)
        self.A_R=eye(3)
        self.d_R=zeros([3,3])
        self.Taps=zeros(3)
        self.base_mva=base_mva*ones(3)
        self.base_kv=base_kv*ones(3)
        self.base_A=self.base_mva*1000./(self.base_kv)
        self.childrenBuses=[] # buses influenced by change of taps. They are all nodes connected to self.to_bus. until a transformer is found.
        
    def calculate_a_R(self):
##        for vr in vr3.voltregs:
##            vr.calculate_a_R()
        
        for index in range(3):
            self.voltregs[index].Vs=self.from_bus.E[index]*self.from_bus.base_kv[index]*1000
            self.voltregs[index].I_line=self.to_bus.totalI[index]
            self.voltregs[index].I_from=self.to_bus.totalI[index]*self.base_A[index]
            self.voltregs[index].calculate_a_R()
            self.a_R[index][index]=self.voltregs[index].a_R
            self.Taps[index]=self.voltregs[index].Tap
        self.A_R=inv(self.a_R)
        self.d_R=inv(self.a_R)
        
    def updateIfrom(self):
        self.I_from=self.I_line#dot(self.d_R,self.I_line) # in pu
        for index in range(3):
            self.voltregs[index].I_from=self.I_from[index]*self.base_A[index]

    def updateBackward(self):
        self.updateIto()
        self.updateIline()
        self.updateIfrom()
        
    def updateToBusBaseKV(self):
        self.to_bus.base_kv=dot(self.A_R,self.from_bus.base_kv)
        # this valid for transformer only. Line and regulator should
        # overide this function.    
##        #print "sdfsdfsd===================================================="
##        if self.from_bus.whatIam()!='RootBus3':
##            self.from_bus.E=self.E_to + dot(self.Z,self.I_line)
##            #print self,'E_from0 =',sum(self.E_from)/3
##    def updateForward(self):
##        self.to_bus.E=dot(self.a_R,self.from_bus.E)
    def updateEto(self):
        self.E_to=dot(self.A_R,self.E_from) # in pu
        self.to_bus.E=self.E_to
        self.updateToBusBaseKV()
        
    def updateChildren_base_kv(self):
        for b in self.childrenBuses:
            b.base_kv=self.to_bus.base_kv
    
        
    def updateEfrom(self):
        self.E_from=self.from_bus.E
        for index in range(3):
            self.voltregs[index].Vs=self.E_from[index]
##    def updateIline(self):
##        self.I_line=self.I_to*self.base_i # in A
##    def updateIto(self):
##        self.I_to=self.to_bus.totalI*self.base_i # in A
##        return self.I_to  
