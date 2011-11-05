from pylon import Bus, Generator, Branch, Case
#from networkx import DiGraph, shortest_path_length
from scipy.sparse import hstack, vstack, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, splu
from numpy import cos, array, sin,ones,zeros,pi, arccos,dot, eye,sqrt,conj, sum, rank, append
from numpy.linalg import solve, inv
from Line3 import Line3
from Bus3 import Bus3

# Constants
a=cos(2*pi/3)+1j*sin(2*pi/3)
A=array([1.0,a,a*a])# creating three 120-shifted unity voltage.
P=array([[1,1,1],[1,a,a*a],[1,a*a,a]])/3
Pinv=array([[1,1,1],[1,a*a,a],[1,a,a*a]])

class Trafo3(Line3):
    """ A model for transformer based on:
    
        Xiao, P.; Yu, D. & Yan, W. A unified three-phase transformer model 
        for distribution load flow calculations Power Systems, 
        IEEE Transactions on, 2006, 21, 153 - 159
        
    """    
    nmode='NEGATIVE MODE'
    pmode='POSITIVE MODE'
    def __init__(self,from_bus=Bus3(),to_bus=Bus3(),yt=0,base_mva= 1.0,base_kv=12.0,case_base_mva=1.0,case_base_kv=12.0,E_high= 12.47, E_low=4.16, type_connection='YgYg',type='Step-Down',ps=0, alpha=1.0,beta=1.0 ):
        Line3.__init__(self,from_bus,to_bus)
        # constants
        self.YI=eye(3)*yt#
        self.YII=array([[2.,-1.,-1.],[-1.,2.,-1.],[-1.,-1.,2.]])*yt/3.
        self.YIII=array([[-1.,1.,0.],[0.,-1.,1.],[1.,0.,-1.]])*yt/sqrt(3.0)
        self.YIV=array([[1.,-1.,0.],[-1.,2.,-1.],[0.,-1.,1.]])*yt/3.
        self.YV=array([[-1.,1.,0.],[0.,-1.,1.],[0.,0.,0.]])*yt/sqrt(3.0)
        self.YVI=array([[1.,0.,0.],[0.,1.,0.],[0.,0.,0.]])*yt
        self.type=type
        self.yt=yt
	self.Yabc=self.YI.copy()# phase admitance matrix
	self.Zabc=inv(self.Yabc) # phase impedance matrix
	self.Z012=dot(dot(P,self.Zabc),Pinv)
        self.type_connection=type_connection
        self.base_mva=base_mva
        self.base_kv=base_kv
        self.E_high=E_high
        self.E_low=E_low
        self.E_from0=0.0 # set initial zero sequence voltage to 0.
        #self.to_bus.base_kv=self.from_bus.base_kv*self.E_low*sqrt(3)/(self.E_high)
        self.ps=ps # how much is secondary side voltage is shifted in radians
        
        # It is important to remember that in the following models, YD connection 
        # will INCREASE phase angle by 30 degrees when we move from Y side to D side.
        # this is a DEFAULT mode named POSITIVE MODE
        
        # When I evaluate the 4 bus feeder test case provided by IEEE, I have found 
        # that It follow reverse convention. In other words, as we move from Y to D
        # side, phase angle are DECREASED by 30 degrees.
        # This is named NEGATIVE MODE
        if self.type_connection=='YgYg': # tested. there are more to model
            self.Ypp=self.YI.copy()
            self.Yps=-self.YI.copy()
            self.Ysp=-self.YI.copy()
            self.Yss=self.YI.copy()
        if self.type_connection=='YgY':
            self.Ypp=self.YII.copy()
            self.Yps=-self.YII.copy()
            self.Ysp=-self.YII.copy()
            self.Yss=self.YII.copy()
        if self.type_connection=='YgD':
            self.Ypp=self.YI.copy()
            self.Yps=self.YIII.copy()
            self.Ysp=self.YIII.transpose().copy()
            self.Yss=self.YII.copy()
        if self.type_connection=='YD':
            self.Ypp=self.YII.copy()
            self.Yps=self.YIII.copy()
            self.Ysp=self.YIII.transpose().copy()
            self.Yss=self.YII.copy()
        if self.type_connection=='DYg':
            self.Ypp=self.YII.copy()
            self.Yps=self.YIII.copy()
            self.Ysp=self.YIII.transpose().copy()
            self.Yss=self.YI.copy()
        if self.type_connection=='DD':
            self.Ypp=self.YII.copy()
            self.Yps=-self.YII.copy()
            self.Ysp=-self.YII.transpose().copy()
            self.Yss=self.YI.copy()
            #print 
        if self.type_connection=='UV':
            self.Ypp=self.YVI.copy()
            self.Yps=self.YV.copy()
            self.Ysp=self.Yps.transpose().copy()
            self.Yss=self.YIV.copy()
        if type=='Step-Up':
            self.reverseMode()
        self.calculateSubmatrices()
    def calculateSubmatrices(self):    
        # additional submatrices to handle singularities    
        if self.type_connection!='YgYg':
            self.Ysp1=self.Ysp.copy()
            self.Ysp2=self.Ysp.copy()
            self.Yss1=self.Yss.copy()
            self.Yss2=self.Yss.copy()
            self.Yps2=self.Yps.copy()
            self.Ypp2=self.Ypp.copy()
            self.Ysp1[2]=ones(3)
            self.Yss1[2]=zeros(3)
            self.Ysp2[2]=zeros(3)
            self.Yps2[2]=ones(3)
            self.Yss2[2]=ones(3)
            self.Ypp2[2]=zeros(3) # used during voltage update  i.e. forward sweep
        #print hstack([self.Ypp,self.Yps])
        #self.YT=vstack([hstack([self.Ypp,self.Yps]),hstack([self.Ysp,self.Yss])])
    def reverseMode(self):
        self.Yps=self.Yps.transpose().copy()
        self.Ysp=self.Ysp.transpose().copy()
        self.calculateSubmatrices()
        
    def updateIto(self):
        """
        Update injected secondary current of transformer a line of to-side of a line.
        """
        self.I_to=self.to_bus.totalI.copy() # the minus sign is important because by definition it is an injected current to the secondary side.
        self.I_tod=self.I_to.copy()
        self.I_tod[2]=0 # I_tod is  used for connection in which Ysp is singular
    def updateIline(self):
        self.I_line =self.I_from.copy()
        return self.I_line
    def updateIfromtemp(self):
        self.updateIntEfrom()
        # remember to update primary voltage before updating this primary current.
        self.I_from_temp= dot(self.Ypp,self.E_temp)+dot(self.Yps,self.E_to)
        if self.type_connection=='YgD':
            # here the idea is that zero component primary current is caused only by
            # primary current. Hence we substract zero componet of secondary voltage.
            # hmm can we actually justify this logic?
            #Ep=self.E_to -sum(self.E_to/3.0)
            self.I_from_temp= dot(self.Ypp,self.E_temp)+dot(self.Yps,self.E_to)
            #print 'updateIfromtemp', self.type_connectio

        self.S_from=self.E_temp*conj(self.I_from)
        self.I_from=self.I_from_temp
    def updateBackward(self):
        """
        Perform the backward sweep process.
        """
        #1. Sum up line segement currents
        self.updateIto()
        #2. calculate intermediate primary voltage
        self.updateIntEfrom()
        #3. Calculate Power at primary side
        # 3.a. Calculate primary current
        self.updateIfromtemp()
        #self.I_from_temp= dot(self.Ypp,self.E_temp)+dot(self.Yps,self.E_to)
        self.calculateSfrom()

    
    def updateForward(self):
        # 1. Calculate primary current injaction
        # update voltage at from_bus
        self.updateEfrom()
        self.E_from0=sum(self.E_from)/3.0
        #self.I_from=conj(self.S_from/self.E_from)

        self.I_fromd=self.I_from.copy()
        self.I_fromd[2]=0
        # 2. Calculate secondary voltage
        # 2.a. if Yps is invertible
        if self.type_connection=='YgYg':
            #b=csr_matrix(self.I_to-dot(self.Ysp,self.from_bus.E))
            b=csr_matrix(self.I_from-dot(self.Ypp,self.from_bus.E))
            #a=csr_matrix(self.Yss)
            a=csr_matrix(self.Yps) # This Yps is singular for other type of connection
            self.E_to=spsolve(a,b)
        # 2.b. if Yps is NOT invertibel i.e. singular 
        if self.type_connection!='YgYg':
            # Calculate secondary voltage positive+negative segquence part
            b=csr_matrix( self.I_fromd -dot(self.Ypp2,self.E_from))
            a=csr_matrix(self.Yps2)
            #self.E_to0=sum(self.to_bus.E)/3.
            self.E_to_temp=spsolve(a,b) # this contains NO zero sequence part
            # Calculate secondary voltage zero sequence part
            self.E_to0=0 # this is a potential error
            # a. when Yss is NOT singular ==> YgYg and DYg
            if self.type_connection=='DYg': #or self.type_connection=='YgYg':
                b=csr_matrix( -self.I_to -dot(self.Ysp,self.E_from))
                a=csr_matrix(self.Yss)
                self.E_to_temp=spsolve(a,b)                
                #self.E_to0=sum(self.E_to)/3.0 #this leads to nonconvergences
                #print '================Update Forwad', self.type_connection
                # in the rest of connection types, zero components of secondary
                # voltage depend on downstream network situation. In situation
                # where there is multigrounded system, it is NOT zero and therefore
                # this assuption is wrong.
                # It is unfortunate that this feature is not yet developed as it
                # is very situational.
                # however, in most of situation this assumption is correct.       
            if self.type_connection=='YgD' or self.type_connection=='YD' :
                self.E_to0=sum(self.E_to)/3.0
##                print 'find me!'
            if self.type_connection=='UV':
##                print 'UV Forward sweep'
                b=csr_matrix( self.I_fromd -dot(self.Ypp2,self.E_from))
                a=csr_matrix(self.Yps2)
                self.E_to_temp=spsolve(a,b)                
                self.E_to0=sum(self.E_to)/3
            # Sum all sequences
            self.E_to=self.E_to_temp+self.E_to0
            #print self.E_to
    def updateIntEfrom(self):
        if self.type_connection=='YgYg':
            b=csr_matrix(-self.I_to - dot(self.Yss, self.E_to))
            a=csr_matrix(self.Ysp)
            self.E_temp=spsolve(a,b)
        if self.type_connection!='YgYg':
            # for other type of connection Ysp is singular. therefore some tricks
            # have to be done. 
            # for definition of Yss1 and Ysp1 see self.__init__()
            b=csr_matrix(-self.I_tod - dot(self.Yss1, self.E_to))
            a=csr_matrix(self.Ysp1)
            
            #self.E_temp=spsolve(a,b) # this contains no zero sequence components
            #print"HHHHHHHHHHHHHHHHHHH"
            # below zero sequence component is added. This is ok for Yg Delta connection
            # what I do not understand at this moment is how about DD and DYg.
##            if self.type_connection=='YgD':
##                self.E_from0=sum(self.E_from)/3.0 
##                #print "****************YgD********************"
            self.E_temp=spsolve(a,b) + self.E_from0
            # this function should not update the from_bus.E because 
            # in current update only current is updated.
        if self.type_connection=='UV':
            A=csr_matrix(self.Ysp[0:2,0:2])
            b=-self.I_tod - dot(self.Yss1, self.E_to)
            B=csr_matrix(b[0:2])
##            print B.shape
            Epm=spsolve(A,B)
            self.E_temp=append(Epm,0) +self.E_from0 # what is the value of E_from[2] here?
        return self.E_temp# This is the intermediate primary voltage.
    def calculateSfrom(self):
        self.S_from=dot(self.E_temp, conj(self.I_from_temp) )
        return self.S_from
    
    
    def updateToBusBaseKV(self):
        self.to_bus.base_kv=self.from_bus.base_kv*self.E_low/self.E_high
        # this valid for transformer only. Line and regulator should
        # overide this function.    
    
    def updateLineLoss(self):
        Ep=self.E_from
        Ip=self.I_from
        Es=self.E_to
        Is=self.I_to
        self.line_loss=Ep*conj(Ip)-Es*conj(Is)
    
    
##    def updateEto(self):
##        #self.I_from=conj(self.S_from/self.E_from)
##        self.I_fromd=conj(self.S_from/self.E_from)
##        self.I_fromd[2]=0.
##        if self.type_connection=='YgYg':
##            #b=csr_matrix(self.I_to-dot(self.Ysp,self.from_bus.E))
##            b=csr_matrix(self.I_from-dot(self.Ypp,self.from_bus.E))
##            #a=csr_matrix(self.Yss)
##            a=csr_matrix(self.Yps)
##            self.E_to=spsolve(a,b)
##        
##            # there are three alternative for calculating zerosequence secondary
##            # voltage
##                # 1. for DYg connection
##        elif self.type_connection=='DYg':
##            b=csr_matrix( -self.I_tod -dot(self.Ysp2,self.E_from))
##            a=csr_matrix(self.Yss2)
##            self.E_to_temp=spsolve(a,b)
##            self.E_to0=sum(self.E_to)/3.0
##            self.E_to=spsolve(a,b) + self.E_to0
##            print "Updating E_to"
##            #b=csr_matrix(self.I_tod-dot(self.Ysp2,self.from_bus.E))
##        else:
##                # 2. The zero sequence voltage is a function of network grounding downstream.
##                #  here it is assume to be zero. This is valid for most of the case.
##            self.E_to0=0
##
##            b=csr_matrix(self.I_fromd-dot(self.Ypp2,self.from_bus.E))
##            #a=csr_matrix(self.Yss2)        
##            a=csr_matrix(self.Yps2)
##            self.E_to=spsolve(a,b) + self.E_to0
##        return self.E_to
