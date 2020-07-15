# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:53:27 2020

@author: Gavin
"""

import numpy as np
from scipy.optimize import fsolve

def cstr(in_file1,c0,Vol,v0):
    
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
            
    
    # Indicate increment for volume and total reactor volume
    V = np.linspace(0,Vol,1000)    # L

    # Define function for reaction
    def rxn(c, c0, V): 
        cTGE, cDGE, cMGE, cCAT, cEtOH, cEE, cG, cTGO, cDGO, cMGO = c
        
        # Define residence time
        tau = V/v0
        
        # E: EtOH phase
        cTGEi = c0[0]
        cDGEi = c0[1]
        cMGEi = c0[2]
        cCATi = c0[3]
        cEtOHi = c0[4]
        cEEi = c0[5]
        cGi = c0[6]
        # O: Oil phase
        cTGOi = c0[7]
        cDGOi = c0[8]
        cMGOi = c0[9]
        
        # Elementary reactions
        r1 = (k1*cTGE*cEtOH*cCAT)-(k1r*cDGE*cEE*cCAT)
        r2 = (k2*cDGE*cEtOH*cCAT)-(k2r*cMGE*cEE*cCAT)
        r3 = (k3*cMGE*cEtOH*cCAT)-(k3r*cG*cEE*cCAT)
        # Mass transport btwn EtOH and Oil phase
        rmTG = kLTG*aTG*(cTGO-cTGE)
        rmDG = kLDG*aDG*(cDGO-cDGE)
        rmMG = kLMG*aMG*(cMGO-cMGE)
        
        # Concentration of respective species
        dTGEdt = -r1+rmTG
        dDGEdt = r1-r2+rmDG
        dMGEdt = r2-r3+rmMG
        dCATdt = 0
        dEtOHdt = -r1-r2-r3
        dEEdt = r1+r2+r3
        dGdt = r3
        dTGOdt = -rmTG
        dDGOdt = -rmDG
        dMGOdt = -rmMG
        
        # Mole balances in CSTR for each species
        f1 = cTGEi-cTGE+dTGEdt*tau
        f2 = cDGEi-cDGE+dDGEdt*tau
        f3 = cMGEi-cMGE+dMGEdt*tau
        f4 = cCATi-cCAT+dCATdt*tau
        f5 = cEtOHi-cEtOH+dEtOHdt*tau
        f6 = cEEi-cEE+dEEdt*tau
        f7 = cGi-cG+dGdt*tau
        f8 = cTGOi-cTGO+dTGOdt*tau
        f9 = cDGOi-cDGO+dDGOdt*tau
        f10 = cMGOi-cMGO+dMGOdt*tau
        
        return[f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]
    
    # Run calculation of conversion over course of many iterations
    lines = []
    for i in range(len(V)):        
        if i == 0:
            lines.append(fsolve(rxn,[0,0,0,0,0,0,0,0,0,0], args=(c0,V[i])))
        else:
            lines.append(fsolve(rxn,lines[-1], args=(c0,V[i])))
    conc = np.vstack(lines)
    
    return conc