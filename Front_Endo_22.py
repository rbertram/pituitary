# Front_Endo_22.py

# This Python file describes a collection of electrically coupled 
# pseudo-plateau bursting cells as described in the paper
# "Network Properties of Electrically Coupled Bursting Pituitary
# Cells", Mehran Fazli and Richard Bertram, Frontiers in Endocrinology,
# 13:936160, 2022.
 
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from random import random


def ode_solver(IC,Nt,ggj):
    # parameter values
    vn=-5
    kc=0.12
    ff=0.005
    vca=60
    vk=-75
    vl = -50.0
    gk=2.5
    cm=5
    gbk=1
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
#    ggj=0.002
    dt=0.5

    #//////////////////////////////////
    # LOOP OVER INITIAL CONDITIONS



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
            ibk=gbk*bb[i]*(vv[i]-vk)
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
            ibk=gbk*(bb[i] + 0.5*k1b)*((vv[i] + 0.5*k1v)-vk)
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
            ibk=gbk*(bb[i] + 0.5*k2b)*((vv[i] + 0.5*k2v)-vk)
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
            ibk=gbk*(bb[i] + 0.5*k3b)*((vv[i] + 0.5*k3v)-vk)
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
            if t>=Nt-20000:
                for i in range(num_cell-1):
                    for h in range(i+1,num_cell):
                        if vv[i]>=-40 and vv[h]>=-40:
                            T_ij[i,h]=T_ij[i,h]+1

                        if i==0 and h==1:
                            if vv[i]>=-40:
                                T_n[i]=T_n[i]+1

                            if vv[h]>=-40:
                                T_n[h]=T_n[h]+1;
                        if i==0 and h>1:
                            if vv[h]>=-40:
                                T_n[h]=T_n[h]+1

    for i in range(num_cell-1):
        for h in range(i+1,num_cell):
            C_ij[i,h]=T_ij[i,h]/((T_n[i]*T_n[h])**0.5)
            C_ij[h,i]=T_ij[i,h]/((T_n[i]*T_n[h])**0.5)
    return C_ij, vtim, ntim, ctim, btim
# a two node network 
SN_sim=[[0, 1], 
        [1, 0]];
G_sim_sc = nx.Graph()
num_cell=len(SN_sim)
for i in range(num_cell):
    for j in range(num_cell):
        if SN_sim[i,j] == 1: 
              G_sim_sc.add_edge(i,j)
ggj=0.002 # gap junction conductance
dt=0.5
sec=100
num_set=100 # number of initial conditions
Ntim=int((sec*1000*(1/dt)))
sim=np.zeros((num_set,num_cell,num_cell))
ICall=np.zeros((num_set,4,num_cell))
vall=np.zeros((num_set,num_cell,Ntim))
nall=np.zeros((num_set,num_cell,Ntim))
call=np.zeros((num_set,num_cell,Ntim))
ball=np.zeros((num_set,num_cell,Ntim))

for setn in range(1,num_set+1):
    for i in range(0,num_cell):
        ICall[setn-1,0,i]= -90*(random())+20
        ICall[setn-1,1,i]= random()/2
        ICall[setn-1,2,i]= random()/2
        ICall[setn-1,3,i]= random()
    AA=ICall[setn-1,:,:]
    sim[setn-1,:,:], vall[setn-1,:,:], nall[setn-1,:,:], call[setn-1,:,:], ball[setn-1,:,:]=ode_solver(AA,Ntim,ggj)
