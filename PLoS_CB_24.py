# PLoS_CB_24.py

# This Python file describes a network of bursting and spiking cells that
# are electrically coupled, as described in the paper "Conversion of
# Spikers to Bursters in Pituitary Cell Networks: Is it Better to 
# Disperse for Maximum Exposure or Circle the Wagons?", Mehran Fazli
# and Richard Bertram, PLoS Computaitonal Biology, 20:e1011811, 2024.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
from sklearn.metrics import pairwise_distances
import seaborn as sns
from random import random
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})
plt.rcParams.update({'font.size': 20})
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)



def ode_solver(IC,Nt,ggj,gbk):
    vn=-5
    kc=0.12
    ff=0.005
    vca=60
    vk=-75
    vl = -50.0
    gk=2.5
    cm=5
#    gbk=1
    gca=2.1
    gsk=2
    vm=-20
    vb=-5
    sn=10
    sm=12
    sbk=2
    taun=30
    taubk=5
    ks=0.4
    alpha=0.0015
    gl=0.2
#    ggj=0.05
    dt=0.5
    
    #//////////////////////////////////
    # LOOP OF INITIAL CONDITION



    vtim=np.zeros((num_cell,Nt))
    ntim=np.zeros((num_cell,Nt))
    ctim=np.zeros((num_cell,Nt))
    btim=np.zeros((num_cell,Nt))

    vv=np.zeros(num_cell)
    nn=np.zeros(num_cell)
    cc=np.zeros(num_cell)
    bb=np.zeros(num_cell)
    T_n=np.zeros(num_cell)
    T_ij=np.zeros((num_cell,num_cell))
    C_ij=np.zeros((num_cell,num_cell))
    #/////////////// initialize the variable matrices

    for i in range(0,num_cell):
        vv[i]= IC[0,i]
        nn[i]= IC[1,i]
        cc[i]= IC[2,i]
        bb[i]= IC[3,i]


    for t in range(0,Nt): #/////////////////////////////////////////////////SOLVER TIME-STEP ITERATION

        ninf=0
        bkinf=0
        minf=0
        cinf=0
        ica=0
        isk=0
        ibk=0
        ikdr=0
        il=0

        k1v=0
        k2v=0
        k3v=0
        k4v=0
        k1n=0
        k2n=0
        k3n=0
        k4n=0
        k1c=0
        k2c=0
        k3c=0
        k4c=0
        k1b=0
        k2b=0
        k3b=0
        k4b=0

        vv_old=np.zeros(num_cell)
        for i in range(0,num_cell):
            vv_old[i]=vv[i]


        for i in range(0,num_cell): #/////////////////////// LOOP IN CELLS

            #//////////////////////////////////////////////////////////////////////////// FIRST STEP RK
            ninf=1/(1+np.exp((vn-vv[i])/sn))
            bkinf=1/(1+np.exp((vb-vv[i])/sbk))
            minf=1/(1+np.exp((vm-vv[i])/sm))
            cinf=((cc[i])**2)/(((cc[i])**2)+ks*ks)

            ica=gca*minf*(vv[i]-vca)
            isk=gsk*cinf*(vv[i]-vk)
            ibk=gbk[i]*bb[i]*(vv[i]-vk)
            ikdr=gk*nn[i]*(vv[i]-vk)
            il = gl*(vv[i]-vl)

            igj=0;
            for h in range(0,num_cell):
                igj=igj+SN_sim[i][h]*ggj*(vv[i]-vv_old[h])

            k1v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k1n = dt*((ninf-nn[i])/taun)
            k1c = dt*(-ff*(alpha*ica+kc*cc[i]))
            k1b = dt*(-(bb[i]-bkinf)/taubk)
            #//////////////////////////////////////////////////////////////////////////// SECOND STEP RK
            ninf=1/(1+np.exp((vn-(vv[i] + 0.5*k1v))/sn))
            bkinf=1/(1+np.exp((vb-(vv[i] + 0.5*k1v))/sbk))
            minf=1/(1+np.exp((vm-(vv[i] + 0.5*k1v))/sm))
            cinf=((cc[i] + 0.5*k1c)**2)/(((cc[i] + 0.5*k1c)**2)+ks*ks)

            ica=gca*minf*((vv[i] + 0.5*k1v)-vca)
            isk=gsk*cinf*((vv[i] + 0.5*k1v)-vk)
            ibk=gbk[i]*(bb[i] + 0.5*k1b)*((vv[i] + 0.5*k1v)-vk)
            ikdr=gk*(nn[i] + 0.5*k1n)*((vv[i] + 0.5*k1v)-vk)
            il = gl*((vv[i] + 0.5*k1v)-vl)
            igj=0
            for h in range(0,num_cell):
                igj=igj+SN_sim[i][h]*ggj*((vv[i] + 0.5*k1v)-vv_old[h])

            k2v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k2n = dt*((ninf-(nn[i] + 0.5*k1n))/taun)
            k2c = dt*(-ff*(alpha*ica+kc*(cc[i] + 0.5*k1c)))
            k2b = dt*(-((bb[i] + 0.5*k1b)-bkinf)/taubk)
            #//////////////////////////////////////////////////////////////////////////////// THIRD STEP RK
            ninf=1/(1+np.exp((vn-(vv[i] + 0.5*k2v))/sn))
            bkinf=1/(1+np.exp((vb-(vv[i] + 0.5*k2v))/sbk))
            minf=1/(1+np.exp((vm-(vv[i] + 0.5*k2v))/sm))
            cinf=((cc[i] + 0.5*k2c)**2)/(((cc[i] + 0.5*k2c)**2)+ks*ks)

            ica=gca*minf*((vv[i] + 0.5*k2v)-vca)
            isk=gsk*cinf*((vv[i] + 0.5*k2v)-vk)
            ibk=gbk[i]*(bb[i] + 0.5*k2b)*((vv[i] + 0.5*k2v)-vk)
            ikdr=gk*(nn[i] + 0.5*k2n)*((vv[i] + 0.5*k2v)-vk)
            il = gl*((vv[i] + 0.5*k2v)-vl);
            igj=0;
            for h in range(0,num_cell):
                igj=igj+SN_sim[i][h]*ggj*((vv[i] + 0.5*k2v)-vv_old[h])

            k3v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k3n = dt*((ninf-(nn[i] + 0.5*k2n))/taun)
            k3c = dt*(-ff*(alpha*ica+kc*(cc[i] + 0.5*k2c)))
            k3b = dt*(-((bb[i] + 0.5*k2b)-bkinf)/taubk)
            #////////////////////////////////////////////////////////////////////////////////// FOURTH STEP RK 
            ninf=1/(1+np.exp((vn-(vv[i] + 0.5*k3v))/sn))
            bkinf=1/(1+np.exp((vb-(vv[i] + 0.5*k3v))/sbk))
            minf=1/(1+np.exp((vm-(vv[i] + 0.5*k3v))/sm))
            cinf=((cc[i] + 0.5*k3c)**2)/(((cc[i] + 0.5*k3c)**2)+ks*ks)

            ica=gca*minf*((vv[i] + 0.5*k3v)-vca)
            isk=gsk*cinf*((vv[i] + 0.5*k3v)-vk)
            ibk=gbk[i]*(bb[i] + 0.5*k3b)*((vv[i] + 0.5*k3v)-vk)
            ikdr=gk*(nn[i] + 0.5*k3n)*((vv[i] + 0.5*k3v)-vk)
            il = gl*((vv[i] + 0.5*k3v)-vl)
            igj=0;
            for h in range(0,num_cell):
                igj=igj+SN_sim[i][h]*ggj*((vv[i] + 0.5*k3v)-vv_old[h])

            k4v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k4n = dt*((ninf-(nn[i] + 0.5*k3n))/taun)
            k4c = dt*(-ff*(alpha*ica+kc*(cc[i] + 0.5*k3c)))
            k4b = dt*(-((bb[i] + 0.5*k3b)-bkinf)/taubk)
            #////////////////////////////////////////////////////////////////////////////////// FINAL STEP RK

            vv[i] = vv[i] + (1.0/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
            nn[i] = nn[i] + (1.0/6.0)*(k1n + 2*k2n + 2*k3n + k4n)
            cc[i] = cc[i] + (1.0/6.0)*(k1c + 2*k2c + 2*k3c + k4c)
            bb[i] = bb[i] + (1.0/6.0)*(k1b + 2*k2b + 2*k3b + k4b)
            vtim[i,t]=vv[i]
            ntim[i,t]=nn[i]
            ctim[i,t]=cc[i]
            btim[i,t]=bb[i]
            
    return vtim, ntim, ctim, btim
##############################################################################
##############################################################################
##############################################################################
##############################################################################




# adjacency matrix
SN_sim=[[0, 1], 
        [1, 0]];

# the network can be expanded for example:
'''
SN_sim=[[0, 1, 0], 
        [1, 0, 1], 
        [0, 1, 0]];
'''

G_sim_sc = nx.Graph()
num_cell=len(SN_sim)
for i in range(num_cell):
    for j in range(num_cell):
        if SN_sim[i][j] == 1: 
              G_sim_sc.add_edge(i,j)
color_opt=['tab:red','tab:blue']
cell_type_opt=['burster', 'spiker']

cols=[color_opt[0],color_opt[1]]# color_opt[0] for burster and color_opt[0] for spiker, need update with network changes 
cell_type=[cell_type_opt[0],cell_type_opt[1]]# cell_type_opt[0] for burster and cell_type_opt[0] for spiker, need update with network changes 

plt.figure(figsize=(4,4))

nx.draw_planar(G_sim_sc, with_labels = False, font_size=30, font_color='w', node_color=cols, node_size=2000, edge_color="gray", width=5)
plt.axis('off')
plt.show()


gbk=np.array([1,0]) # 1 for burster and 0 for spiker, need update with network changes 
ggj=0.05 # medium coupling strength
# ggj=0 for disconnected cells

dt=0.5
sec=100



Ntim=int((sec*1000*(1/dt)))
IC=np.zeros((4,num_cell))
    
for i in range(0,num_cell):
    IC[0,i]= -60
    IC[1,i]= 0.1
    IC[2,i]= 0.1
    IC[3,i]= 0.1

vtraj, ntraj, ctraj, btraj=ode_solver(IC,Ntim,ggj,gbk)



######################################################################


fig=plt.figure(figsize=(18,12))
    
ax = fig.add_subplot(3,1,1)
xx=np.arange(0, 2.5, 0.0005)
for i in range(num_cell):
    plt.plot(xx,vtraj[i,Ntim-5000:Ntim], c=cols[i], lw=4, label=cell_type[i])   
#plt.xlabel('time (sec)')
plt.ylabel('$V$ (mV)', fontsize=25)
plt.ylim([-70, 20])
plt.legend(loc='upper left', fontsize = 15)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.setp(ax.get_xticklabels(), fontsize=25)
plt.setp(ax.get_yticklabels(), fontsize=25)

ax = fig.add_subplot(3,1,2)
for i in range(num_cell):
    plt.plot(xx,ctraj[i,Ntim-5000:Ntim], c=cols[i], lw=4, label=cell_type[i])  
#plt.xlabel('time (sec)')
plt.ylim([0.25, 0.37])
plt.ylabel('$c$ ($\mu$ M)', fontsize=25)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.setp(ax.get_xticklabels(), fontsize=25)
plt.setp(ax.get_yticklabels(), fontsize=25)

ax = fig.add_subplot(3,1,3)
alp=5
half=0.6
for i in range(num_cell):
    plt.plot(xx,1/(1+np.exp(-alp*((ctraj[i,Ntim-5000:Ntim]-0.27)/(0.352-0.27)-half))), c=cols[i], lw=4, label=cell_type[i])  

plt.xlabel('time (sec)', fontsize=25)
plt.ylim([0, 1])
plt.yticks([0.2, 0.5, 0.8])
plt.ylabel('$s$', fontsize=25)# secretion
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.setp(ax.get_xticklabels(), fontsize=25)
plt.setp(ax.get_yticklabels(), fontsize=25)

plt.show()
