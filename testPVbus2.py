from bfs_pfdg import bfs_pfdg,bfs_pf
import cProfile
from numpy import array, angle, pi, sqrt, dot, conj, ones, real
from pylon import Bus
from Bus3 import Bus3
from RootBus3 import RootBus3
from Line3 import Line3
from Trafo3 import Trafo3
from Case3 import Case3
from Load3 import Load3
from Generator3 import Generator3
from tool import P_to_S, polarToRect
from VoltReg3 import VoltReg3
from VoltReg3 import VoltReg3, VoltReg1
from tool import *
import inspect
thisfilename=inspect.getfile( inspect.currentframe() )
print '\n\nExecuting:', thisfilename


a=cos(2*pi/3)+1j*sin(2*pi/3)
A=array([1.0,a,a*a])/3.0# creating three 120-shifted unity voltage
zeroload=Load3() # default load

MVAbase=6.0
kVbase=12.47/sqrt(3) # line-to-ground kV
Zbase=kVbase**2/MVAbase#in the present of transformer, zbase as well as kVbase on both side of transformer are differents.

# Buses
busa=Bus()
bus1=RootBus3(busa,busa,busa,base_kv=kVbase)
bus2=Bus3(busa,busa,busa,base_kv=kVbase)
bus3=Bus3(busa,busa,busa,base_kv=kVbase)
bus4=Bus3(busa,busa,busa,base_kv=kVbase)
busvr=Bus3(busa,busa,busa,base_kv=kVbase)
bus5=Bus3(busa,busa,busa,base_kv=kVbase)

# load at bus4

##Sa=.75/MVAbase
##pfa=0.85
##Sb=1./MVAbase
##pfb=0.9
##Sc=1.25/MVAbase
##pfc=0.95
Sa=.175/MVAbase
pfa=0.8
Sb=.195/MVAbase
pfb=0.9
Sc=.107/MVAbase
pfc=0.75
Sa=P_to_S(Sa*pfa,pfa)
Sb=P_to_S(Sb*pfb,pfb)
Sc=P_to_S(Sc*pfc,pfc)
S3=array([Sa,Sb,Sc])
load1=Load3(Sa,Sb,Sc,connection='wye',bus=bus4)
bus1.loads=[zeroload]
bus2.loads=[load1]
bus3.loads=[load1]
bus4.loads=[load1]
bus123=[bus1,bus2,busvr,bus3,bus4,bus5]
busvr.loads=[zeroload]
gen1=Generator3(bus1)
gen2=Generator3(bus4,type='PV')
gen3=Generator3(bus5,type='PV')
gen2.updateS(ones(3)*(.01+0.00j))
gen3.updateS(ones(3)*(.0+0.01j))
gen3.E=gen3.E
# Load at bus5
bus5.loads=[load1]

# line data
# from substation to transfomer
zd=array([[0.1414 + 1j*0.5353, 0.0361 + 1j*0.3225, 0.0361 + 1j*0.2752],
         [0.0361 + 1j*0.3225, 0.1414 + 1j*0.5353, 0.0361 + 1j*0.2955],
         [0.0361 + 1j*0.2752, 0.0361 + 1j*0.2955, 0.1414 + 1j*0.5353] ]) # Ohm
#zd=eye(3,3)
Z1=zd/(Zbase)
branch1=Line3(bus1,bus2,Z1)

#transformer data
t_connection='Step-Down'
t_MVA=6.0
t_kVLL_high=12.47
t_zbase= t_kVLL_high**2/t_MVA
t_kVLL_low=2.4*sqrt(3)
t_R=0.01 # pu
t_X=0.06 # pu
t_a=(MVAbase/t_MVA)*(t_kVLL_high/kVbase)**2 # adjustment factor
t_z=(t_R + 1j* t_X)*t_a
t_y=1/t_z
trafo1=Trafo3(bus2,busvr,yt=t_y,base_mva=t_MVA,type_connection='DYg',E_high=t_kVLL_high,E_low=t_kVLL_low)

kVbase3=t_kVLL_low*kVbase/t_kVLL_high # Line to ground

# updating base kv
bus3.base_kv=kVbase3*ones(3)
bus4.base_kv=kVbase3*ones(3)

# from transformer to load line
zd1=array([[0.1907 + 1j*0.5035, 0.0607 + 1j*0.2302, 0.0598 + 1j*0.1751],
         [0.0607 + 1j*0.2302, 0.1939 + 1j*0.4885, 0.0614 + 1j*0.1931],
         [0.0598 + 1j*0.1751, 0.0614 + 1j*0.1931, 0.1921 + 1j*0.4970]])# Ohm
zbase3=kVbase3**2/MVAbase

# voltage regulator
Rd=7.3
Xd=14.2

vr=VoltReg1(Bus(),Bus(),Rd=Rd, Xd=Xd)
vr.calculate_a_R()
vrs=[vr,vr,vr]
vr3=VoltReg3(voltregs=vrs,base_mva=MVAbase,base_kv=kVbase3 )
busvr.base_kv=kVbase3*ones(3)
vr3.from_bus=busvr
vr3.to_bus=bus3
vr3.childrenBuses=[bus4]
N_PT=20
Vset=121
Zd= Rd +1j*Xd
for vr in vr3.voltregs:
    vr.Zd=Zd
    vr.updateVset(Vset)
    vr.Zd_ohm=vr.Zd/vr.CTs
    vr.N_PT=N_PT
    vr.CTp=1000
    vr.CT=vr.CTp/vr.CTs

Z2=zd1/zbase3
branch3=Line3(bus3,bus4,Z2,base_kv=kVbase3)
# line from bus3 to bus5
branch4=Line3(bus3,bus5,Z2,base_kv=kVbase3)
# branchess
branch12=[branch1,trafo1,vr3,branch3,branch4]
kasus=Case3(base_mva=MVAbase, base_kv=kVbase, buses=bus123, branches=branch12,generators=[gen1,gen2,gen3])

bs=bfs_pfdg(kasus)
bs.tolerance=1e-6
[g.S for g in bs.case.generators]
bs.solve_dg()
[g.S for g in bs.case.generators]
#vr3.calculate_a_R()
#rn=kasus.getRoot()
#bs.solve()
#bs.solve_dg()
bs.printResult()
print

##if bs.checkConvergence():
##    print "\nThe power flow calculation was converged after %.0f iterations and %.3f seconds." % (bs.itercount,bs.duration)
##    print 'The tolerance was set to ',bs.tolerance,'.'
##    print "transformator connection type was", trafo1.type_connection,'.'
##
##    E=[b.E*b.base_kv for b in kasus.buses]
##    E2=[bus2]
##    E2=[b.E*b.base_kv for b in E2]
##    print 'Voltages are as follows:\n'
##    printVoltagell(E2)
##    E34=[bus3, bus4]
##    E34=[b.E*b.base_kv for b in E34]
##    printVoltage(E34)
##    printCurrent(kasus)
##    print '\nEnd voltage at 120 v base' , abs(bus4.E*bus4.base_kv*1000/20)
##
