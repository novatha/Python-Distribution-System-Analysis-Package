from pylon import Branch, Bus
from scipy import dot, cos, pi, sin, array, ones, zeros, sqrt, real, imag, arccos, angle, floor, ceil, eye
from numpy import round

class VoltReg1(Branch):
    """
    An implementation of single phase regulator type B as explained in [1]
    Reference
    [1] Kersting, William H.Distribution system modeling and analysis, CRC, 2002
        ISBN 0-8493-0812-7
    """
    def __init__(self,from_bus,to_bus, bandwidth=2.0, N_PT=20.0, CTp=700., CTs=5, Rd=3.0, Xd=9.0, Vset=120, type='B'):
        Branch.__init__(self,from_bus,to_bus)    
        self.from_bus=from_bus # self is attached from bus
        self.to_bus = to_bus # self is attached to bus
        self.bandwidth=bandwidth # my bandwidth in Volt
        self.Vset=Vset # the desired voltage on a 120 V base voltage to be held at the regulating point
        self.Vset_min= self.Vset-self.bandwidth/2.
        self.Vset_max= self.Vset+self.bandwidth/2.
        self.N_PT=N_PT # my potential transformer ratio
        self.CTp=CTp # my primary current rating of the current transformer  in Ampere
        self.CTs=CTs # my secondary current rating of the current transformer  in Ampere
        self.CT=self.CTp/self.CTs
        self.Rd=Rd # Rd --> R' the compensator equvalent setting in Volt
        self.Xd=Xd # Xd --> R' the compensator equvalent setting in Volt
        self.Zd=self.Rd + 1j*self.Xd
        self.Zd_ohm =self.Zd/self.CTs
        self.Tap=0
        self.a_R=1
        self.I_comp=0
        self.I_line=0
        self.I_from=0
        self.Vs=1.0+1j*0
        self.Is=0
        self.VL=1+1j*0
        self.IL=0
        self.V_reg=1+1j*0
        self.V_drop=1+1j*0
        self.type=type
        # Calculate tap position
        # Next please find a_R, eq 7.75 
        # Calculate voltage output
    
    def calculate_a_R(self):
        # a. calculate the actual line current
        #I_line = self.to_bus.totalI
        
        # b. current in compensator
        self.I_comp= self.I_from/self.CT
        #print "I_comp %0.4f / %0.1f "% (abs(self.I_comp),angle(self.I_comp)*180/pi)
##        print "I_from %0.4f / %0.1f / CT %0.0f"% (abs(self.I_from),angle(self.I_from)*180/pi, self.CT)
        # c. input voltage of the compensator is
        
        self.V_reg = self.Vs/self.N_PT
##        print"V_s", self.Vs
##        print "V_reg %0.1f / %0.1f "% (abs(self.V_reg),angle(self.V_reg)*180/pi)
        
        # d. voltage drop in the compensator impedance is
        self.V_drop = self.Zd_ohm * self.I_comp

        # e. the voltage across the voltage relay
        self.V_R = self.V_reg - self.V_drop

        # f. recall that on a 120-V base, one step change on the regulator changes the voltage 0.75 V
        ## please remember that this is valid for type B only. See ch 7 for more explanation
        if abs(self.V_R) < self.Vset_min:
            self.Tap=round((self.Vset_min-abs(self.V_R))/0.75)
            action=['raise',self.Tap]
            self.a_R = 1 - 0.00625*self.Tap
        elif abs(self.V_R) > self.Vset_max:
            self.Tap=round((self.Vset_max-abs(self.V_R))/0.75)
            action=['lower',self.Tap]
            self.a_R = 1 - 0.00625*self.Tap
    def updateCoefficients(self):
        if self.type=='A':
            self.a=1/self.a_R
            self.b=0
            self.c=0
            self.d=self.a_R
        elif self.type=='B':
            self.a=self.a_R
            self.b=0
            self.c=0
            self.d=1/self.a_R
    
    def forwardSweep(self):
        self.updateCoefficients(self)
        self.VL=(self.Vs - self.b*self.IL)/self.a

    def backwardSweep(self):
        self.updateCoefficients(self)
        self.Is=self.c*self.VL + self.d*self.IL
    def updateVset(self,Vset):
        self.Vset=Vset
        self.Vset_min= self.Vset-self.bandwidth/2.
        self.Vset_max= self.Vset+self.bandwidth/2.
# Testing block used during implementation process.

# 1. based on example 7.4 on page 170
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
busa=Bus()
busb=Bus()


## Substation transformer data
kVArated=5000
connection='DYg'
Vprimary=115
Vsecondary=4.16

## The equivalent line impedance from regulator to the load center 
Z=0.3 +1j*0.9


Vs = 4160./sqrt(3) # the rated line-to-ground voltage ot the substation transformer

# potential transformer ratio to provide 120 V to the compensator circuit
N_PT = 2400./120.

# rated current of the substatio transformer
Irated=kVArated/(sqrt(3)*4444.16)

# select a primary Current Transformer rating as 700 A
CTp=700
# if the current at the compensator is chosen to be 5 A, the CT ratio is
CTs=5
CT=CTp/CTs

# R' and X' setting of the compensator in V
Zd=Z*CTp/N_PT # in V. This NOT an impedance but compensator setting
Zd_ohm =Zd/CTs

# 2. based on example 7.5 on page 172
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load data
kVAload = 2500 # kVA
Vload=4.16 # kV
pf= 0.9 # lagging pf
load_angle=-arccos(pf)

# regulator setting
Zd= 10.5 + 1j* 31.5 # compensator setting in V 
Rd= real(Zd)
Xd=imag(Zd)
Vset = 120. # compensator desired voltage in V
bandwidth = 2. # in V
Vset_min= Vset-bandwidth/2.
Vset_max= Vset+bandwidth/2.
####
# Determine the tap position of regulator tha will hold the load center voltage at the desired voltage and within bandwidth
###

# a. calculate the actual line current

I_line = (cos(load_angle)+1j*sin(load_angle))*kVAload/(sqrt(3)*Vload)

# b. current in compensator
I_comp= I_line/CT

# c. input voltage of the compensator is
V_reg = Vs/N_PT

# d. voltage drop in the compensator impedance is
V_drop = Zd_ohm * I_comp

# e. the voltage across the voltage relay
V_R = V_reg - V_drop

# f. recall that on a 120-V base, one step change on the regulator changes the voltage 0.75 V

if V_R < Vset_min:
    Tap=floor((Vset_min-abs(V_R))/0.75)
    action=['raise',Tap]
    a_R = 1 - 0.00625*Tap
elif V_R > Vset_max:
    Tap=ceil((Vset_max-abs(V_R))/0.75)
    action=['raise',Tap]
    a_R = 1 - 0.00625*Tap
### g. regulator model coefficients
##
##a= a_R
##b=0
##c=0
##d=1/a_R
##
### 3. Example 7.6 on page 173
###~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Using the results of Examples 7.5, calculate the actual voltage at the load center 
### assuming the  2500 kVA at 4.16 kV is measured at the substation transformer 
### low-voltage terminals
##
### a. the actual voltage and current at the load side terminal of regulator are
##VL=Vs/a
##Is=I_line
##IL=Is/d
##
### b. the actual line-to-ground voltage at the load center is
##
##VLC = VL - Z*IL
##
### c. the line-to-ground voltage at the load centes on a 120-V base
##
##VLC_120 = VLC/N_PT
##
##vr=VoltReg1(busa,busb,Rd=Rd, Xd=Xd)
##vr.Vs=Vs
##vr.I_line=I_line
##vr.calculate_a_R()