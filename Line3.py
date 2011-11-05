from pylon import Branch
from numpy import dot, cos, pi, sin, array, ones, zeros, conj
from numpy.linalg import inv
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
class Line3(Branch):
    """ Three phase line model
    """
    def __init__(self, from_bus, to_bus, Z=ones([3, 3]),Y=zeros([3,3]),base_kv=1.0):
        Branch.__init__(self,from_bus,to_bus)
        self.Z=Z # This is the series element of the line impedance
        self.Y=Y # This is the shunt element of the line impedance
	self.Z012 = dot(dot(inv(Pinv),Z),Pinv)
        self.base_kv=base_kv
        self.E_from = from_bus.E#zeros(3)
        self.E_to=to_bus.E#zeros(3)
        self.I_from = zeros(3)
        self.I_to=zeros(3)
        #self.Iline=zeros(3)
        self.I_line=self.I_to+dot(self.Y,self.E_to)
        self.line_drop=zeros(3,dtype=complex)
	self.line_loss=zeros(3,dtype=complex)
        #self.A=eye(3, 3) + dot(Z, Y/2)
        #self.B=Z
        #self.C=Y+dot(dot(Y, Z,), Y)
        #self.D=eye(3, 3) + dot(Y/2, Z)
    def getEfrom(self):
        result=self.E_to + dot(self.Z,self.I_line)
        #result=dot(self.A, self.E_to) + dot(self.B, self.I_to)
        return result
    def updateEfrom(self):
        self.E_from=self.from_bus.E
        return self.E_from
    def updateIfrom(self):
        self.I_from=self.I_line +dot(self.Y/2,self.E_from)
    def updateIline(self):
        self.I_line=self.I_to + .5*dot(self.Y,self.E_to)
    def updateIto(self):
        self.I_to=self.to_bus.totalI
        return self.I_to  
    def updateEto(self):
        # the following method is NOT applicable if self.to_bus is a PVbus.
        self.E_to=self.E_from - dot(self.Z,self.I_line)
        if self.to_bus.type!='PV':
            pass # please do something here! 7 June 2011
            # 1. It may be better NOT to update voltages of PV buses here.
            # 2. or simply set it again to scheduled value.
            # 3. then update the reactive power injection instead. See baldwin  p 1187 eq 21
        self.to_bus.E=self.E_to
        return self.E_to
    def updateBackward(self):
        self.updateIto()
        self.updateIline()
        self.updateIfrom()
        #print "sdfsdfsd===================================================="
        if self.from_bus.whatIam()!='RootBus3':
            self.from_bus.E=self.E_to + dot(self.Z,self.I_line)
            #print self,'E_from0 =',sum(self.E_from)/3
    def updateForward(self):
        self.updateEfrom()
        self.updateEto()
    def updateLineDrop(self):
        self.line_drop=dot(self.Z,self.I_line)
        return self.line_drop
    def __repr__(self):
        strs=self.name+'\n'
        return strs
#        str2='E_from = ' +str(self.E_from)+ ' I_from = '+str(self.I_from)+'\n'
#        str3='I_to   = ' +str(self.E_to) + ' I_to   = ' + str(self.I_to)+'\n'    
    def backwardSweep(self):
        pass

#        strs=str1+str2+str3
        return strs
    
    def updateToBusBaseKV(self):
        self.to_bus.base_kv=self.from_bus.base_kv 
        # this valid for line only. Transfomer and regulator should
        # overide this function.
    def updateLineLoss(self):
        Iline=self.I_line
        Zline=self.Z
        deltaV=dot(Zline,Iline)
        self.line_loss=deltaV*conj(Iline)
    
