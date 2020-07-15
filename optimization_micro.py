# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 09:36:36 2020

@author: Gavin
"""

# Import modules
import numpy as np
import PFR
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pandas as pd

saveDir = "graphs"

mesh_size = 100

# Total channel (reactor) length in m
l_lowerlim = 0.25   # 50
l_upperlim = 5      # 200

# Total flowrate in mL/hr
tot_flowrate = 68100000

# Reactor and channel dimensions in m
# Total channel (reactor) length is specified in the function itself, which is dependent on the specified channel (reactor) length
switchback_len = 50e-3
channel_width = 1200e-6
channel_depth = 1200e-6
addl_len = 10e-3
chip_width = switchback_len+1e-2
chip_thick = channel_depth+2400e-6

# Constants required for determining pumping power
# dynamic viscosity
mu = 0.01985   # Ns/m^2 (from soybean oil @ 60C, from interpolation)
# density
rho = 0.917   # kg/m^3 (from soybean oil)
# local loss (pressure drop) coefficient
KL = 0.2   # for 180 deg bend, flanged
# flow coefficient
Cv = 56.9   # for square channel with laminar flow
# discharge coefficient
Cd = 0.5959
# inlet to microreactors (full pipe diameter in m)
Di = 0.1016
# pump efficiency
eff = 0.75   # approx value from sample pump curve in unit ops
# gravity
g = 9.81   # m/s^2

 

# Number of individual microreactor units
units_lowerlim = 200000     # 200000
units_upperlim = 2000000    # 2000000

# Reactor volume in μL
V_lowerlim = l_lowerlim*(channel_depth*channel_width)/(1e-9)
V_upperlim = l_upperlim*(channel_depth*channel_width)/(1e-9)

# Flow rate in mL/hr
v0_lowerlim = (tot_flowrate)/(units_lowerlim)
v0_upperlim = (tot_flowrate)/(units_upperlim)


V = np.linspace(V_lowerlim,V_upperlim,mesh_size)
v0 = np.linspace(v0_lowerlim,v0_upperlim,mesh_size)
V, v0 = np.meshgrid(V, v0)


conv = np.empty((len(V),len(v0)))
perc_yield = np.empty((len(V),len(v0)))
DGrelconc = np.empty((len(V),len(v0)))
MGrelconc = np.empty((len(V),len(v0)))
EErelconc = np.empty((len(V),len(v0)))
for m in range(len(V)):
    for n in range(len(v0)):
    
        # Input file containing rate constants
        in_file1 = "microreactor_constants_basecat"
        # Input file containing initial concentration
        in_file2 = "feed1"
        # Input file containing reactor network scheme w/ reactor volumes
        in_file3 = "PFR"
        
                
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
        
        conc = PFR.pfr(in_file1,c0,(V[m,n]*1e-6),(v0[m,n]*(1e-3)/3600))
        conc = np.vstack(conc)

        # Calculate product stream results
        conv[m,n] = 100*(c0[7]-conc[:,7][-1])/(c0[7]) 
        perc_yield[m,n] = 100*conc[:,5][-1]/(c0[7]*3)
        totconc=0
        for i in range(0,len(conc[-1])):
            totconc += conc[-1][i]
        DGrelconc[m,n] = 100*((conc[:,1][-1]+conc[:,7][-1])/totconc)
        MGrelconc[m,n] = 100*((conc[:,2][-1]+conc[:,9][-1])/totconc)
        EErelconc[m,n] = 100*(conc[:,5][-1]/(totconc-conc[:,3][-1]-conc[:,4][-1]-conc[:,6][-1]))
        
        # Calculate footprint and pumping power
        units = tot_flowrate/v0
        channel_len = V*(1e-9)/(channel_depth*channel_width)
        num_switchbacks = (channel_len/switchback_len)
        chip_len = (num_switchbacks)*(2*channel_width)+addl_len
        
        # Calculate orifice characteristics
        beta = ((channel_depth+channel_width)/2)/(Di)
        
        footprint = units*chip_len*chip_width*chip_thick
        
        delta_p_major = (Cv*mu*channel_len*v0*(1e-3)*(1e-3)/3600)/(2*(channel_depth*channel_width)**2)
        delta_p_minor = (KL*num_switchbacks*rho*(v0*(1e-3)*(1e-3)/3600)**2)/(2*(channel_depth*channel_width)**2)
        #delta_p_orif = (rho/2)*(1-(beta**4))*((v0*(1e-3)*(1e-3)/3600)/(Cd*channel_depth*channel_width))**2
        delta_p = delta_p_major + delta_p_minor
        pump_power = units*delta_p*((v0*(1e-3)*(1e-3)/3600)/(eff*g))
        
        # Make dataframe with summary of results
        labels=["# Microreactor Units","Channel Length (m)",
                "Total Flow Rate (mL/hr)","Flowrate thru each Channel (mL/hr)",
                "Channel Depth (um)", "Channel Width (um)","% Coversion TG",
                "% Yield FAME","DG Rel Conc (mol%)","MG Rel Conc (mol%)",
                "FAME Purity (mol%)","Footprint (m^3)","Pumping Power (W)",
                "Pressure drop (Pa)"]
        data = {"# Microreactor Units": units[m,n],
                "Channel Length (m)": channel_len[m,n], 
                "Total Flow Rate (mL/hr)": tot_flowrate, 
                "Flowrate thru each Channel (mL/hr)": v0[m,n], 
                "Channel Depth (um)": channel_depth, 
                "Channel Width (um)": channel_width, 
                "% Coversion TG": conv[m,n], 
                "% Yield FAME": perc_yield[m,n], 
                "DG Rel Conc (mol%)": DGrelconc[m,n], 
                "MG Rel Conc (mol%)": MGrelconc[m,n], 
                "FAME Purity (mol%)": EErelconc[m,n], 
                "Footprint (m^3)": footprint[m,n], 
                "Pressure drop (Pa)": delta_p[m,n], 
                "Pumping Power (W)": pump_power[m,n]}
        labels2 = ["cTGE","cDGE","cMGE","cCAT","cEtOH","cEE","cG","cTGO","cDGO","cMGO"] 
        data2 = {"cTGE": conc[:,0][-1],
                 "cDGE": conc[:,1][-1],
                 "cMGE": conc[:,2][-1],
                 "cCAT": conc[:,3][-1],
                 "cEtOH": conc[:,4][-1],
                 "cEE": conc[:,5][-1],
                 "cG": conc[:,6][-1],
                 "cTGO": conc[:,7][-1],
                 "cDGO": conc[:,8][-1],
                 "cMGO": conc[:,9][-1]}      
                
        if m==0 and n==0:
            existing_df = pd.DataFrame(columns=labels)
            df = existing_df.append(data,ignore_index=True)
            existing_df2 = pd.DataFrame(columns=labels2)
            df2 = existing_df2.append(data2,ignore_index=True)
        else:
            df = df.append(data,ignore_index=True)
            df2 = df2.append(data2,ignore_index=True)

def plotter(varName,axisName,title,manual_lowerlim=None,manual_upperlim=None,color=None,rotation=None):
    # Plot the surface
    fig = plt.figure(figsize=(9,7))
    ax = fig.gca(projection='3d') 
    if color:
        surf = ax.plot_surface(channel_len, (units*1e-3), varName, cmap=color)
    else:
        surf = ax.plot_surface(channel_len, (units*1e-3), varName, cmap=cm.viridis)
    ax.set_xlim(l_lowerlim,l_upperlim)
    ax.set_ylim((units_lowerlim*1e-3),(units_upperlim*1e-3))
    if manual_lowerlim and manual_upperlim:
        ax.set_zlim(manual_lowerlim,manual_upperlim)
    elif manual_lowerlim:
        ax.set_zlim(manual_lowerlim,np.amax(varName))
    elif manual_upperlim:
        ax.set_zlim(np.amin(varName),manual_upperlim)
    else:
        ax.set_zlim(np.amin(varName),np.amax(varName))
        
    ax.set_zlabel(axisName)
    ax.set_ylabel('$10^3$ microreactor units')
    ax.set_xlabel('Reactor channel length (m)')
    ax.zaxis.label.set_size(13)
    ax.yaxis.label.set_size(13)
    ax.xaxis.label.set_size(13)
    
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    if rotation:
        ax.view_init(20,rotation)
    else:
        ax.view_init(20,225)
    plt.savefig("{}/{}.png".format(saveDir,title))
    return


plotter(conv, '% Coversion of Triglyceride','conversion',manual_upperlim=100)
plotter(perc_yield, '% Yield of FAME','percentyield')
plotter(DGrelconc, 'Diglyceride Relative Concentration (mol%)','DGrelconc',rotation=135)
plotter(MGrelconc, 'Monoglyceride Relative Concentration (mol%)','MGrelconc',rotation=135)
plotter(EErelconc, 'FAME Purity in Separated Oil Phase (mol%)','EErelconc')
plotter(footprint, 'Footprint ($m^3$)','footprint',color='inferno')
plotter(pump_power, 'Pumping Power (W)','pumppower',color='inferno',rotation=135)
df.to_csv("{}/summary.csv".format(saveDir))
df2.to_csv("{}/conc.csv".format(saveDir))
plt.show()