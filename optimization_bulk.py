# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 09:36:36 2020

@author: Gavin
"""

# Import modules
import numpy as np
import CSTR
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pandas as pd

saveDir = "basecatbulk"

# Total flowrate in L/hr
v0 = 68100

# Reactor volume in L
Vol = 200000

V_lowerlim = 0
V_upperlim = Vol

V = np.linspace(V_lowerlim,V_upperlim,1000)

reac_constants = ["bulkreactor2_constants_basecat","bulkreactor3_constants_basecat","bulkreactor1_constants_basecat","bulkreactor4_constants_basecat"]
reac_conditions = ["Sunflower oil / MeOH / 1wt% KOH / 65°C / 600 rpm [35]","Sunflower oil / MeOH / 1wt% NaOH / 60°C / 400 rpm [36]","Sunflower oil / MeOH / 0.2wt% NaOH / 50°C / 300 rpm [37]","Sunflower oil / MeOH / 1wt% KOH / 50°C / 500 rpm [38]"]

conv = np.empty((len(reac_constants),len(V)))
perc_yield = np.empty((len(reac_constants),len(V)))
DGrelconc = np.empty((len(reac_constants),len(V)))
MGrelconc = np.empty((len(reac_constants),len(V)))
EErelconc = np.empty((len(reac_constants),len(V)))

for m in range(len(reac_constants)):
    
        # Input file containing rate constants
        in_file1 = reac_constants[m]
        # Input file containing initial concentration
        in_file2 = "feed1"
        # Input file containing reactor network scheme w/ reactor volumes
        in_file3 = "CSTR"
        
                
        # Read in relevant values from initial concentration input file
        with open(in_file2,"r") as f:
            lines=np.array(f.readlines())        
                
        c0 = []    
        for i in np.arange(len(lines)):
            if lines[i].startswith("cTGE ="):
                c0.append(float(lines[i].strip("cTGE =")))
            elif lines[i].startswith("cDGE ="):
                c0.append(float(lines[i].strip("cDGE =")))
            elif lines[i].startswith("cMGE ="):
                c0.append(float(lines[i].strip("cMGE =")))
            elif lines[i].startswith("cCAT ="):
                c0.append(float(lines[i].strip("cCAT =")))
            elif lines[i].startswith("cEtOH ="):
                c0.append(float(lines[i].strip("cEtOH =")))
            elif lines[i].startswith("cEE ="):
                c0.append(float(lines[i].strip("cEE =")))
            elif lines[i].startswith("cG ="):
                c0.append(float(lines[i].strip("cG =")))
            elif lines[i].startswith("cTGO ="):
                c0.append(float(lines[i].strip("cTGO =")))
            elif lines[i].startswith("cDGO ="):
                c0.append(float(lines[i].strip("cDGO =")))
            elif lines[i].startswith("cMGO ="):
                c0.append(float(lines[i].strip("cMGO =")))
        
        conc = CSTR.cstr(in_file1,c0,Vol,(v0/3600))
        conc = np.vstack(conc)

        """
        if m==1:
            print(reac_constants[m])
            labels = ["cTGE","cDGE","cMGE","cCAT","cEtOH","cEE","cG","cTGO","cDGO","cMGO"] 
            df = pd.DataFrame(conc,columns=labels)
            df.to_csv("{}/conc.csv".format(saveDir))
        """
            
        for n in range(len(V)):
            # Calculate product stream results
            conv[m,n] = 100*(c0[7]-conc[n][7])/(c0[7]) 
            perc_yield[m,n] = 100*conc[n][5]/(c0[7]*3)
            totconc=0
            for i in range(0,len(conc[n])):
                totconc += conc[n][i]
            DGrelconc[m,n] = 100*((conc[n][1]+conc[n][7])/totconc)
            MGrelconc[m,n] = 100*((conc[n][2]+conc[n][9])/totconc)
            EErelconc[m,n] = 100*(conc[n][5]/(totconc-conc[n][3]-conc[n][4]-conc[n][6]))
                
def plotter(varName,axisName,title,manual_lowerlim=None,manual_upperlim=None):
    plt.figure(figsize=(10,9))
    for m in range(len(reac_constants)):
        plt.plot(V, varName[m], label=reac_conditions[m])
    plt.xlim(V_lowerlim,V_upperlim)
    if manual_lowerlim and manual_upperlim:
        plt.ylim(manual_lowerlim,manual_upperlim)
    elif manual_lowerlim:
        plt.ylim(manual_lowerlim,np.amax(varName))
    elif manual_upperlim:
        plt.ylim(np.amin(varName),manual_upperlim)
    else:
        plt.ylim(np.amin(varName),np.amax(varName))
        
    plt.ylabel(axisName,fontsize=17)
    plt.xlabel('Reactor volume (L)',fontsize=17)
    #plt.rcParams.update({'font.size': 12})
    # plt.title('Conventional CSTR for Transesterification with Homogeneous Base Catalysis')    
    plt.legend(loc='lower right',fontsize=15)
    plt.savefig("{}/{}.png".format(saveDir,title))
    return


plotter(conv, '% Coversion of Triglyceride','conversion',manual_upperlim=100)
plotter(perc_yield, '% Yield of FAME','percentyield')
plotter(DGrelconc, 'Diglyceride Relative Concentration (mol%)','DGrelconc')
plotter(MGrelconc, 'Monoglyceride Relative Concentration (mol%)','MGrelconc')
plotter(EErelconc, 'FAME Purity in Separated Oil Phase (mol%)','EErelconc')
plt.show()