from numpy import array, r_, Inf, hstack, vstack, abs, max, real, imag
from numpy.linalg import norm
from iterator import Traverser
from time import time
from tool import *
from cetak import *

class bfs_pf(Traverser):
    """Implement backward and forward sweep method
    """
    def __init__(self,case):
        Traverser.__init__(self,case)
        self.max_iter=100
        self.tolerance=1e-4
        self.verbose=False
        self.oldE=array([b.E for b in self.case.buses])
        self.presentE=array([b.E for b in self.case.buses])
        self.itercount=0
        self.duration=0
	self.myfilename=__file__
    def updateCurrent(self):
        """ Performing backward sweep level by level.
        """
        max_level=max(self.case.bus_level.values())

        # start with the maximum level. decrease step bu step.
        k=max_level # the level counter
        while k>0 :
            buses=[bus for bus in self.getBusesAtLevel(k)]
            # print 'k= ',k,'\nbuses =',buses
            for bus in buses:
                # update load current at each bus at level k.
                bus.calc_loadI()
                #print 'loadI = ',bus.loadI
                # update shunt current
                bus.calc_shuntI()
                # update current to branches connected to this bus from higher level.
                bus.calc_childI()
                # update branch curret at each branch ending at level k.
                totalI=bus.calc_totalI()
                for branch in bus.connected_from_branch:
                    branch.I_to=totalI
                    branch.updateBackward()
            # go to the k-1 level
            # print
            k=k-1
    def updateVoltage(self):
        """ Performing forward sweep level by level. starting from level 1
        """

        # calculate max level
        max_level=max(self.case.bus_level.values())
        # start from level 1
        k=1
        while k<=max_level:
            # pick a bus
            for bus in self.getBusesAtLevel(k):
                # pick a branch having the bus as its from_bus
                for branch in bus.connected_from_branch:
                    # update voltages at to_bus
                    branch.updateForward() # line and transformer may have different implementation here
                    # set the bus voltage to branch.E_to
                    branch.to_bus.E=branch.E_to
            # after all buses are updated, go to next level
            k=k+1
    def oneIteration(self):
        self.updateCurrent()
        self.updateVoltage()
        self.presentE=array([b.E for b in self.case.buses])
        #print max(abs(self.oldE-self.presentE))
        #print abs(self.presentE[1])

    def solve(self):
        awal=time()
        repeat= True
        k=0
        while repeat and k<self.max_iter:
            self.oldE=array([b.E for b in self.case.buses])
            self.oneIteration()
##            self.updateCurrent()
##            self.updateVoltage()
            converged=self.checkConvergence()
            repeat=not converged
            #print converged
            k=k+1
        if converged:
            pass
            #print 'Conververged after ',k,'iterations.'
        if not converged:
            pass
            #print 'Sorry! I gave up after ',self.max_iter,' iterations.'
        self.itercount=k
        akhir=time()
        self.duration=akhir-awal # in seconds
	#if self.checkConvergence():
	for br in self.case.branches:
            br.updateLineLoss() # calculate line losses only when iteration is converged.
    def checkConvergence(self):
        mis=self.oldE-self.presentE
        F = r_[abs(mis)]
        normF = max(abs(F))
        if normF < self.tolerance:
            converged = True
        else:
            converged = False
            if self.verbose:
                logger.info("Difference: %.3f" % (normF - self.tolerance))
        return converged

    def printE(self):
        for bus in self.case.buses:
            print bus, bus.E
    def profileMe(self):
        cProfile.runcall('self.solve()', 'bssolvetxt')
    def update_loadS(self):
	for b in self.case.buses:
	    b.loadS=sum([l.S for l in b.loads])

    def printResult(self):
        if self.checkConvergence():
            self.update_loadS()
            kVAbase=self.case.base_mva*1000
	    print "===================================================================="	
	    print " pypdsa: Python Package for Electrical Distribution System Analysis"
#	    print " by Novalio Daratha"
            print "===================================================================="
#            print "\nThe power flow calculation was converged after %.0f iterations and %.3f seconds." % (self.itercount,self.duration)
#            print 'The tolerance was set to ',self.tolerance,'.\n'
            print '\nMVAbase = %.3f'%self.case.base_mva
            print 'kVbase  = %.3f' %self.case.base_kv
            print '\n========'
            print 'BUS DATA'
            print '========'
            table=[[b.name,abs(b.E[0]),abs(b.E[1]),abs(b.E[2]),angle(b.E[0])*180/pi,angle(b.E[1])*180/pi,angle(b.E[2])*180/pi,real(b.sourceS[0]),real(b.sourceS[1]),\
                    real(b.sourceS[2]),imag(b.sourceS[0]),imag(b.sourceS[1]),imag(b.sourceS[2]), real(b.loadS[0]),real(b.loadS[1]),real(b.loadS[2]),\
                    imag(b.loadS[0]),imag(b.loadS[1]),imag(b.loadS[2])] for b in self.case.buses]
            title_row='Bus name\t\tV(pu)\t\t\tAngle(deg) \t\t\tPg(pu)  \tQg(pu) \t\t\tPL(pu) \t\t\tQL(pu)\n                a        b      c        a        b      c         a        b      c      a        b      c      a        b      c        a        b      c'

            table2=[[br.name, br.from_bus, br.to_bus,real(br.I_line[0]),real(br.I_line[1]), real(br.I_line[2]),imag(br.I_line[0]),imag(br.I_line[1]),imag(br.I_line[2]),kVAbase*real(br.line_loss[0]),kVAbase*real(br.line_loss[1]),kVAbase*real(br.line_loss[2]),kVAbase*imag(br.line_loss[0]),kVAbase*imag(br.line_loss[1]),kVAbase*imag(br.line_loss[2])] for br in self.case.branches]
            import sys
            out=sys.stdout
            print title_row
            pprint_table(out,table)

            print '\n========='
            print 'LINE DATA'
            print '========='
	    for br in self.case.branches:
		br.updateLineLoss()

            title_row_line='Name\t\t From Bus  To Bus\t Real Current (pu)\tReactive Current (pu) \t    Active Losses (kw)\t   Reactive Losses (kvar)'
            phase_row_line='\t\t\t\t       a        b        c       a        b        c         a        b        c         a        b        c'
            print title_row_line
            print phase_row_line
            pprint_table(out,table2)
