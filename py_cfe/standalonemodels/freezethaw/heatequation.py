#!/usr/bin/env python
# coding: utf-8

# ## Crank Nicolson discretization of Heat Equation

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import inv

# depth
N = -1

#HFUS = 2.96E6 #for saturated sandy soil
HFUS = 0.3336E06 # latent heat of fusion (from noahmp.F)
DF = np.zeros(1) 
ST = np.zeros(1)
Z = np.zeros(1)

def set_Z(Z_):
    global Z, N
    Z = Z_
    N = len(Z)
    
def set_initial_conditions(initT):
    global ST
    ST = initT


def get_initial_conditions():
    return ST

def set_temperature(T):
    global ST
    ST = T
    
def set_thermal_conductivity(TCond):
    global DF
    DF = TCond
# *************************************************************
#boundary conditions
TBOT = 273.15
TPST = 273.15
OPT_TBOT = 1 # 1 = zero flux, 2 = pescribed temperature

# *************************************************************
def set_boundary_conditions(Ttop, Tbot=-999):
    global TBOT, TPST, OPT_TBOT
    TBOT = Tbot
    TPST = Ttop
    if TBOT == -999:
        OPT_TBOT = 1
    
# Ground heat flux at the top surface
def ground_flux(st):
    
    GHF = - DF[0] * (st - TPST)/ (Z[0])
    return GHF


def AssembleM(N, DT, ST):
    FLUX = np.zeros(N)
    AI = np.zeros(N)
    BI = np.zeros(N)
    CI = np.zeros(N)
    RHS = np.zeros(N)
    LAMBD = np.zeros(N)
    
    # compute matrix coefficient using Crank-Nicolson discretization scheme
    for i in range(0,N):
        if i == 0:
            #h = (Z[i+1] - Z[i])
            h = Z[i]
            DTDZ  = (ST[i+1] - ST[i])/ h
            LAMBD[i] = DT/(4* h * HFUS)
            GHS = ground_flux(ST[i])
            FLUX[i] = LAMBD[i] * (DF[i] * DTDZ + GHS)
            
        elif (i < N -1):
            h = Z[i+1] - Z[i]
            h1 = Z[i] - Z[i-1]
            LAMBD[i] = DT/(4*h* HFUS)
            a_ = - LAMBD[i] * DF[i] /h1
            c_ = - LAMBD[i] * DF[i+1] /h
            b_ = 1 + a_ + c_
            FLUX[i] = -a_ * ST[i-1] + b_ * ST[i] - c_ * ST[i+1]
        elif (i == N-1):
            h = (Z[i] - Z[i-1])
            LAMBD[i] = DT/(4* h* HFUS)
            if OPT_TBOT == 1:
                BOTFLX = 0.
            elif OPT_TBOT == 2:
                DTDZ1 = (ST[i] - TBOT) / ( Z[i] - Z[i-1])
                BOTFLX  = - DF[i] * DTDZ1
            DTDZ = (ST[i] - ST[i-1] )/ h
            FLUX[i]  = LAMBD[i] * (-DF[i]*DTDZ  + BOTFLX ) 
            
    # put coefficients in the corresponding vectors A,B,C, RHS
    for i in range(0,N):
        if i == 0:
            AI[i] = 0
            CI[i] = - LAMBD[i] *DF[i]/Z[i]
            BI[i] = 1 - CI[i]
        elif (i < N-1):
            AI[i] = - LAMBD[i] * DF[i]/(Z[i] - Z[i-1])
            CI[i] = - LAMBD[i] * DF[i+1]/(Z[i+1] - Z[i])
            BI[i] = 1 - AI[i] - CI[i]
        elif (i == N-1):
            AI[i] = - LAMBD[i] * DF[i]/(Z[i] - Z[i-1])
            CI[i] = 0
            BI[i] = 1 - AI[i]
        RHS[i] = FLUX[i]
        
    # add the previous timestep ST to the RHS at the boundaries
    for i in range(N):
        if i ==0:
            RHS[i] = ST[i] + RHS[i]
        elif i == N-1:
            RHS[i] = ST[i] + RHS[i]
        else:
            RHS[i] = RHS[i]
    
    # Set up the enrite NxN matrix 
    temp_list = []
    for c in range(N):
        col = [0.,]*N
        temp_list.append(list(col))

    M = np.matrix(temp_list)

    for i in range(N):
        if i==0:
            M[i,i] = BI[i]
            M[i,i+1] = CI[i]
        elif (i < N-1):
            M[i,i] = BI[i]
            M[i,i+1] = CI[i]
            M[i,i-1] = AI[i]
        else:
            M[i,i] = BI[i]
            M[i,i-1] = AI[i]
    
    return M, RHS

def UpdateT(DT, ST):
    M, RHS = AssembleM(N, DT, ST)

    RHS_T = np.reshape(RHS,(N,1))
    
    Minv = inv(M)
    # For CFE we will need a tridiagonal solver as well
    Soln = np.array(Minv*RHS_T).flatten()
    return Soln

def AdvanceT(dt):
    global ST
    DT = dt
    STL = UpdateT(DT, ST)
    ST = STL
    return ST
