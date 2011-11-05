from bfs_pf2 import bfs_pf

from numpy import array, sign, pi, exp, angle, dot, diag, real, imag, cos, sin, conj,zeros,ones,abs
from numpy.linalg import solve, inv


a=cos(2*pi/3)+1j*sin(2*pi/3)
A=array([1.0,a,a*a])/3.0# creating three 120-shifted unity voltage. Used to calculate positive sequence component

class bfs_pfdg(bfs_pf):
        """
        Implements a backward forward sweep method with Distributed generation modeled as PV bus. [1]
        [1] Thomas L. Baldwin, S. A. L. Distribution Load Flow Methods for Shipboard Power Systems IEEE Transactions on Industry Applications, 2004, 40, 1183-1190
        [2] Khushalani, S.; Solanki, J. & Schulz, N. Development of Three-Phase Unbalanced Power Flow Using PV and PQ Models for Distributed Generation and Study of the Impact of DG Models Power Systems, IEEE Transactions on, 2007, 22, 1019 -1025
        """
        def __init__(self,case):
            bfs_pf.__init__(self,case)
            self.iter_max_dg=20
            self.deltaIabc=[]
            self.deltaI=[]
            self.deltaV=[]
            self.deltaQpv=[]
            self.deltaIqmin=[]
            self.deltaIqmax=[]
            self.Qgnew=[]
            self.case.identifyPVbuses()
            self.setPV=set(self.case.pvbuses) # a set of PV buses
            self.setPVhitlimit=set()
            self.Vmismatch=1.0

        def solve_dg(self):
            repeatmainloop=True
            self.n_iter=0
            npv=len(self.case.pvbuses)
            self.deltaV=zeros(npv)
            while repeatmainloop and self.n_iter<self.iter_max_dg:
                self.n_iter=self.n_iter+1
#                print 'iteration ke',self.n_iter
                #1. run LF with no pvbuses.
                self.solve()
                if not self.checkConvergence():
                    print 'Initial power flow was not succesfull!. Stopping.'
                    repeatmainloop=False
                    break
                if self.checkConvergence():
                    self.case.identifyPVbuses()
                    self.case.calculateZpos()# if there is no PV bus, this impedance matrix should be zero.
                    npv=len(self.case.pvbuses)
    #              print 'npv',npv
                    self.k=0 # number of converged pv buses.
                    self.deltaV=array([abs(b.gens[0].Esch[0]) -abs(dot(A, b.E)) for b in self.case.pvbuses]) # positive sequence voltage mismatch (npv x 1)
                    dVold=self.deltaV
#                    print 'deltaV',self.deltaV
                    self.deltaI=1j*solve(abs(self.case.Z1),self.deltaV)
#                    print "Z1", self.case.Z1
#                    print 'deltaI',self.deltaI
                    self.calculateDeltaIeachPhase()
                    self.evaluateqlimit()
                    if self.setPV==self.setPV.intersection(self.setPVhitlimit):
                        print 'All PV buses have been converted into PQ buses'
                        for g in self.case.generators:
                            g.calculateTerminalCurrent()
                        self.solve()
                        repeatmainloop = False
                        break
                    if abs(self.deltaV.all()) < self.tolerance:
                        repeat=False
                        print 'All PV buses are converged.'
                        for g in self.case.generators:
                            g.calculateTerminalCurrent()
                        self.solve()
                        repeatmainloop=False
                        break
                    for g in self.case.generators:
                        g.calculateTerminalCurrent()
                    self.setPVhitlimit=set()

        def calculateDeltaIeachPhase(self):
            #pv=len(self.deltaI)
            angleV=[angle(b.E) for b in self.case.pvbuses]
            signV=[sign(dV)*ones(3) for dV in self.deltaV]
            h=[exp(1j*(x*pi/2 + y)) for (x,y) in zip(signV,angleV)]
            self.deltaIabc=[x*y for (x,y) in zip(self.deltaI,h)]
#            print (self.deltaIabc)
            # then assign these values to pvbuses
            for indeks in range(len(self.case.pvbuses)):
                    self.case.pvbuses[indeks].Iqpv= self.case.pvbuses[indeks].Iqpv+(self.deltaIabc[indeks])
#                    print 'deltaIabce[indeks]]',(self.deltaIabc)
#                    print 'Iqpv',self.case.pvbuses[indeks].Iqpv
                    self.case.pvbuses[indeks].calc_totalI()
        def evaluateqlimit(self):
            """
            Calculate min and max limit of reactive injected current of a generator modeled as PV
            """
            gens=[b.gens[0] for b in self.case.pvbuses]
            for g in gens:
                Qgold=imag(g.S)
                deltaQgnew=imag(g.bus.E*conj(g.bus.Iqpv))
                Qgnew=sum(Qgold+deltaQgnew)
#                print 'deltaQgnew',deltaQgnew
#                print 'Iqpv',g.bus.Iqpv
                if g.Qmin <= Qgnew <= g.Qmax:
                      #  print 'Good! Qg of',g.name,' is within limits'
                        g.S=g.S+1j*deltaQgnew
                        continue
                if Qgnew < g.Qmin:
                        print 'Minimum Q violation at bus', g.bus
                        g.type='PQ'
                        g.S=real(g.S)+1j*g.Qmin
                        self.setPVhitlimit= self.setPVhitlimit.union(set([g]))
                        g.S=real(g.S)+1j*g.Qmin*ones(3)/3.0
                if Qgnew > g.Qmax:
                        g.type='PQ'
                        g.S=real(g.S)+1j*g.Qmax
                        self.setPVhitlimit= self.setPVhitlimit.union(set([g]))
                        g.S=real(g.S)+1j*g.Qmax*ones(3)/3.0
