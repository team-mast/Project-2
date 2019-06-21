#
# This script will v correct the orbit without any knowledge 
# of the lattice. It will use BPM signals to flatten the orbit
# and the dipole correctors as independent parameters.
#

import sys
import math
import types
import time

from jarray import *
from java.lang import *
from java.util import *
from java.io import *
from java.util import ArrayList

from xal.smf.data import XMLDataManager
from xal.smf import AcceleratorSeqCombo
from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

from xal.smf.impl import BPM, VDipoleCorr

#----------------------------------------
# Calsses of OpenXAL Solver Package 
#----------------------------------------
from xal.extension.solver import Scorer
from xal.extension.solver import Trial
from xal.extension.solver import Variable
from xal.extension.solver import Stopper
from xal.extension.solver import SolveStopperFactory
from xal.extension.solver import ProblemFactory
from xal.extension.solver import Solver
from xal.extension.solver import Problem
from xal.extension.solver.algorithm import SimplexSearchAlgorithm
from xal.extension.solver.algorithm import RandomShrinkSearch
from xal.extension.solver.hint import Hint
from xal.extension.solver.hint import InitialDelta

#===============================================================
#              MAIN PROGRAM
#===============================================================
# read the accelerator & make the sequence
accl = XMLDataManager.loadDefaultAccelerator()

#====== Let's construct accelerator ===========
ccl1 = accl.getSequence("CCL1")
ccl2 = accl.getSequence("CCL2")
ccl3 = accl.getSequence("CCL3")
ccl4 = accl.getSequence("CCL4")
DTL6 = accl.getSequence("DTL6")


#+++++++++++++++++++++++++++++++++++++++++++++
#accSeq = AcceleratorSeqCombo("SEQUENCE", [ccl1,ccl2,ccl3,ccl4])

#DTL1 = accl.getSequence("DTL1")
#DTL2 = accl.getSequence("DTL2")
#DTL3 = accl.getSequence("DTL3")
#DTL4 = accl.getSequence("DTL4")
#DTL5 = accl.getSequence("DTL5")
#DTL6 = accl.getSequence("DTL6")

#+++++++++++++++++++++++++++++++++++++++++++++
accSeq = AcceleratorSeqCombo("SEQUENCE", [DTL6,ccl1,ccl2,ccl3,ccl4])
accSeq2=AcceleratorSeqCombo("SEQUENCE", [ccl1,ccl2,ccl3,ccl4])

BPMPos = []	

bpms = accSeq.getAllNodesOfType("BPM")
for bpm in bpms:
	print "debug bpm=",bpm.getId()," pos[m]=",accSeq.getPosition(bpm)
	BPMPos.append(accSeq.getPosition(bpm))
        print "BPM Position",BPMPos
print "=============================================="

bpm_y_min=[]
bpm_y_max=[]
dcvs = accSeq2.getAllNodesOfType("DCV")
for dcv in dcvs:
	print "debug dcv=",dcv.getId(),"  pos[m]=",accSeq2.getPosition(dcv)
	#---- initial set all dcv to zero as a test
        dcv.setField(-.012)
        time.sleep(1.)
        for bpm in bpms:
            x=bpm.getYAvg()
            print x
            bpm_y_min.append(x)
        dcv.setField(.012)
        time.sleep(1.)
        for bpm in bpms:
            y=bpm.getYAvg()
            print y
            bpm_y_max.append(y)

	
print "=============================================="
print "=========="

#-------------------------------------------------------------------
#


scenario=Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_LIVE)

particle_tracker=AlgorithmFactory.createParticleTracker(accSeq)
particle_tracker.setRfGapPhaseCalculation(True)
probe=ProbeFactory.createParticleProbe(accSeq,particle_tracker)

scenario.resync()
scenario.setProbe(probe)
scenario.run()

traj=scenario.getProbe().getTrajectory()

kick_change=[]
response_Mtrx=[]
for dcv in range(len(dcvs)):
        ##set field stuff here
    field_min =-.12
    field_max=.12
    state_init=traj.initialState()
    momentum=state_init.getMomentum()/1.0e9 #GeV/c
    Leff_0=dcvs[0].getEffLength()
    kick_change.append(Leff_0*(field_max-field_min)/3.33564*momentum)

print "kick change"
print kick_change
print "BPM MAX"
print bpm_y_max
print "BPM MIN"
print bpm_y_min
Res_mtrx=[]

for dcv in range(len(dcvs)):
    Res_mtrx.append((bpm_y_max[dcv]-bpm_y_min[dcv])/kick_change[dcv])
print Res_mtrx

print "AT THE ENNNNNND"
    
        

#-------------------------------------------------------------------
# Start of orbit optimization ======================================
#-------------------------------------------------------------------

#-------------------------------------
# Scorer Interface implementation
#-------------------------------------
class OrbitScorer(Scorer):
        """ 
        Calculate the avreage deviation of the orbit (at BPM positions)
        from the center by using data from trial_point variables
        variables is Java's ArrayList()
        """
        def __init__(self,bpms,dcvs,Res_mtrx):
                self.bpms = bpms
                self.dcvs = dcvs
                self.variables = variables
                #--------------------------------
                self.diff2_min = Double.MAX_VALUE
                self.count = 0
        def score(self,trial,variables):
                self.count += 1
                #---------------------------------
                diff2_dcv = 0.
                limit = 0.98
                penalty_rate = 1000000.0
                bpm_y=[]
                Sum=0.
                bpm_y_new=[]
                for bpm in bpms:
                        y=bpm.getYAvg()
                        bpm_y.append(y)
                for dcv_ind in range(self.variables.size()):
                        var = self.variables.get(dcv_ind)
                        kick=trial.getTrialPoint().getValue(var)
                        min_kick = self.variables.get(dcv_ind).getLowerLimit()*limit
                        max_kick = self.variables.get(dcv_ind).getUpperLimit()*limit
                        if(kick < min_kick):
                                diff2_dcv = penalty_rate*(kick - min_kick)**2
                        if(kick > max_kick):
                                diff2_dcv = penalty_rate*(kick - max_kick)**2
                        Sum+=(Res_mtrx[dcv_ind]*kick)
                bpm_y_new.append(bpm_y+Sum)
                diff2_dcv /= self.variables.size()
                #-----------------------------------------------
                diff2 = 0.
                for y in bpm_y_new:
                        diff2 += y**2
                diff2 /= len(bpm_y_new)
                #--------------------------------------
                diff2 += diff2_dcv
                if(self.diff2_min > diff2):
                        self.diff2_min = diff2
                        print "debug solver count = ",self.count,"  sqrt(diff2)= %12.5g "%math.sqrt(diff2)
                return diff2
def getDCH_Field_Arr_for_Trial(self,trial):
                #------ return dch field array for the trial point 
                field_arr = []
                for dcv_ind in range(self.variables.size()):
                        var = self.variables.get(dcv_ind)
                        field =  trial.getTrialPoint().getValue(var)
                        field_arr.append(field)
                return field_arr

#---- Initial step in parameters. During optimization
#---- these steps will be reduced inside the optimizer. 
delta_hint = InitialDelta()

#---- optimizing variabes
variables = ArrayList()

field_max =  0.012
field_min = -0.012

field_step = (field_max - field_min)/30

for dcv_ind in range(len(dcvs)):
        dcv = dcvs[dcv_ind]
        field = dcv.getField()
        var = Variable(dcv.getId(),field, field_min, field_max)
        variables.add(var)
        delta_hint.addInitialDelta(var,field_step)

scorer = OrbitScorer(bpms,dcvs,Res_mtrx)

n_iterations = 200
maxSolutionStopper = SolveStopperFactory.maxEvaluationsStopper(n_iterations)
#solver = Solver(SimplexSearchAlgorithm(),maxSolutionStopper)
solver = Solver(RandomShrinkSearch(),maxSolutionStopper)
problem = ProblemFactory.getInverseSquareMinimizerProblem(variables,scorer,0.00000001)
problem.addHint(delta_hint)
solver.solve(problem)

#------- get optimization results
trial = solver.getScoreBoard().getBestSolution()
dch_field_arr = scorer.getDCH_Field_Arr_for_Trial(trial)


print "========== put new dcv fields into VA ======"
#---- send all results to EPICS
for dcv_ind in range(len(dcvs)):
	dcv = dcvs[dcv_ind]
	field = dcv_field_arr [dcv_ind]
	dcv.setField(field)
	print "dcv=",dcv.getId()," field= %+8.6f "%field	
print "Done." 



sys.exit()
