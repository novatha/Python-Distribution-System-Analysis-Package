from numpy import arccos, cos, sin, angle, pi, array, dot
from Trafo3 import Trafo3

def P_to_S(P,pf):
    """
    Calculate complex power given real power and lagging power factor
    """
    theta=arccos(pf)
    Samp=P/pf
    S=Samp*cos(theta)+1j*Samp*sin(theta)
    return S
def printPolar(E,kabs=1.0):
    """
    print polar form of a complex number E. angle is in degree
    """
    for bbb in E:
        
        print "%.0f  %.1f " % (abs(bbb*kabs),angle(bbb)*180/pi)



def polarToRect(magnitude,degree):
    """
    convert polar complex number to its rectangular form
    """
    result=magnitude*(cos(degree)+1j*sin(degree))
    return result


def printCurrent(kasus):
    """
	Print currents in each line.
    """
    print "Currents (magnitude in V/ angle in degree)"
    BR=[br for br in kasus.branches if not type(br)==Trafo3 ]
    for br in BR:
        print
        base_i=(kasus.base_mva)/br.from_bus.base_kv
        for index in range(3):
            brb=br.I_line[index]
	    print "%4.1f   %.1f " % (abs(brb)*base_i[index]*1000,angle(brb)*180/pi)

def printResult(E,kasus):
    printVoltage(E)
    printCurrent(kasus)
def printVoltage(E):
    """
    Print Line-to-neutral Voltage

    """
    
    for bb in E:
        printPolar(bb,1000)
        print
        print
def printVoltagell(E):
    """
	Print line-to-line voltage

    """
    D=array([[1.0,-1.0,0.0],[0.0,1.0,-1.0],[-1.0,0.0,1.0]])
##    E=[b.E*b.base_kv for b in kasus.buses]
    for bb in E:
        bb=dot(D,bb)
        for bbb in bb:
            
            print "%.1f   %.1f " % (abs(bbb*1000),angle(bbb)*180/pi)
        print
        print

def printVoltageln(E):
    """
	Print line-to-neutral voltages

    """
    W=array([
    [2.0,1.0,0.0],[0.0,2.0,1.0],[1.0,0.0,2.0]
    ])/3.0
    D=array([[1.0,-1.0,0.0],[0.0,1.0,-1.0],[-1.0,0.0,1.0]])
    for bb in E:
        print bb
        bbll=dot(D,bb)
        bbln=dot(W,bbll)
        print
        for bbb in bbln:
            
            print "%.0f   %.1f " % (abs(bbb*1000),angle(bbb)*180/pi)
        print
    I=[br.I_line for br in kasus.branches]


def printVoltageng():
    """
    Print Line-to-Ground Voltage

    """
    W=array([
    [2.0,1.0,0.0],[0.0,2.0,1.0],[1.0,0.0,2.0]
    ])/3.0
    D=array([[1.0,-1.0,0.0],[0.0,1.0,-1.0],[-1.0,0.0,1.0]])
    Elg=bus2.E    
    Ell=dot(D,Elg)
    Eln=dot(W,Ell)
    Eng=Elg-Eln
    for bbb in Eng:
        print "%.0f   %.1f " % (abs(bbb*1000),angle(bbb)*180/pi)
