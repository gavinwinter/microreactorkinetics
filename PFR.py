# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:53:27 2020

@author: Gavin
"""

import numpy as np
from scipy.integrate import odeint

def pfr(in_file1,c0,Vol,v0):
    
    # Read in relevant values from rate constants input file
    with open(in_file1,"r") as f:
        lines=np.array(f.readlines())
    
    for i in np.arange(len(lines)):
        if lines[i].startswith("k1 ="):
            k1 = float(lines[i].strip("k1 ="))
        elif lines[i].startswith("k1r ="):
            k1r = float(lines[i].strip("k1r ="))
        elif lines[i].startswith("k2 ="):
            k2 = float(lines[i].strip("k2 ="))
        elif lines[i].startswith("k2r ="):
            k2r = float(lines[i].strip("k2r ="))
        elif lines[i].startswith("k3 ="):
            k3 = float(lines[i].strip("k3 =")) 
        elif lines[i].startswith("k3r ="):
            k3r = float(lines[i].strip("k3r ="))
        elif lines[i].startswith("kLTG ="):
            kLTG = float(lines[i].strip("kLTG ="))
        elif lines[i].startswith("kLDG ="):
            kLDG = float(lines[i].strip("kLDG ="))
        elif lines[i].startswith("kLMG ="):
            kLMG = float(lines[i].strip("kLMG ="))
        elif lines[i].startswith("aTG ="):
            aTG = float(lines[i].strip("aTG ="))
        elif lines[i].startswith("aDG ="):
            aDG = float(lines[i].strip("aDG ="))
        elif lines[i].startswith("aMG ="):
            aMG = float(lines[i].strip("aMG ="))
    
    # Define function for reaction
    def rxn(F, V, v0): 
        # E: EtOH phase
        cTGE = F[0]/v0
        cDGE = F[1]/v0
        cMGE = F[2]/v0
        cCAT = F[3]/v0
        cEtOH = F[4]/v0
        cEE = F[5]/v0
        cG = F[6]/v0
        # O: Oil phase
        cTGO = F[7]/v0
        cDGO = F[8]/v0
        cMGO = F[9]/v0
        
        # Elementary reactions
        r1 = (k1*cTGE*cEtOH*cCAT)-(k1r*cDGE*cEE*cCAT)
        r2 = (k2*cDGE*cEtOH*cCAT)-(k2r*cMGE*cEE*cCAT)
        r3 = (k3*cMGE*cEtOH*cCAT)-(k3r*cG*cEE*cCAT)
        # Mass transport btwn EtOH and Oil phase
        rmTG = kLTG*aTG*(cTGO-cTGE)
        rmDG = kLDG*aDG*(cDGO-cDGE)
        rmMG = kLMG*aMG*(cMGO-cMGE)
        
        # dFdV Molar flow rates of respective reactants (in mol/s)
        dTGEdV = -r1+rmTG
        dDGEdV = r1-r2+rmDG
        dMGEdV = r2-r3+rmMG
        dCATdV = 0
        dEtOHdV = -r1-r2-r3
        dEEdV = r1+r2+r3
        dGdV = r3
        dTGOdV = -rmTG
        dDGOdV = -rmDG
        dMGOdV = -rmMG
        
        return[dTGEdV, dDGEdV, dMGEdV, dCATdV, dEtOHdV, dEEdV, dGdV, dTGOdV, dDGOdV, dMGOdV]
    
    # Indicate increment for volume and total reactor volume
    V = np.linspace(0,Vol,1000)    # L
    
    F0 = [i*v0 for i in c0]  # mol/s
    F = odeint(rxn, F0, V, args=(v0,), hmax=0.1)
    # Convert from molar flow rate (mol/s) to concentration
    conc = F[:,:]/v0

    return conc