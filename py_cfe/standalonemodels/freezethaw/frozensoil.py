#!/usr/bin/env python
# coding: utf-8

# ## Frozen Soil Model 
# ### Using Freezing-point depression model (temperature and water content dependent) for phase partitioning
# ### Using Crank-Nicolson scheme to solve Diffusion Equation
# 
# #### Algorithm 
# * Get initial Soil Temperature and Soil Moisture (Total and liquid moisture content) profiles
# * Call PhaseChange method, using freezing-point depression equation, to partitioin total water content into liquid and ice
# * Soil heat capacity is updated inside PhaseChange using current timestep liquid/ice contents
# * Thermal conductivities are computed at the end of the timestep for the next timestep
# 
# #### Problem description



import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation, PillowWriter 
import sys, os
import imageio, re
import importlib as imp
import heatequation as HE
import json
import copy
import pandas as pd
#imp.reload(HE) # if needed


class FrozenSoil():
    def __init__(self, cfg_file=None):
        # define domain and soil properties
        self.NSNOW = 0
        self.ISNOW = 0
        self.NSOIL = 0
        self.DSNSO = []
        self.NL = 0
        self.ZBOT = 0
        self.soil_params = {} # dictionary of soil properties
        self.cfg_file = cfg_file #mode input file, contains pointer to forcing data and model params
        self.forcing_file = ""
        self.real_forcing= False
        # dict of forcing variables, if needed
        #self.forcing_vars = {}

        #time step information
        self.starttime = 0
        self.endtime = 0
        self.dt = 0 # time step

        
    def initialize(self):
        
        self.read_config_file()
        self.theta_w = self.get_theta()
        self.SMCT = np.array(self.theta_w) # total soil moisture content
        self.SMCLiq = np.array(self.theta_w) # liquid water content
        self.dz = self.get_dz()
        
        #***************************************************************
        # **** Initial Soil Temperature, Thermal Conductivity **********
        # **** Set boundary conditions ******
        HE.set_Z(self.Z)
        HE.set_initial_conditions(self.ST)
        
        TCOND = self.thermal_conductivity()
        
        HE.set_thermal_conductivity(TCOND)

        if (len(self.BC_T) > 0):
            HE.set_boundary_conditions(self.BC_T[0])
        else:
            HE.set_boundary_conditions(self.Ttop)

        ## Note: we should call PhaseChange here too to get SMCIce initial content, especially it
        ## will be need when soil T in below freezing at the initial time
        self.PhaseChange()
        
    def read_config_file(self):
        with open(self.cfg_file) as file:
            dat = json.load(file)  

        self.forcing_file = dat["forcing_file"]
        self.Z = dat["Z"]
        self.soil_params['vg_n'] = dat["soil_params"]["vg_n"]
        self.soil_params['vg_m'] = 1. - 1./dat["soil_params"]["vg_n"]
        self.soil_params['vg_alpha'] = dat["soil_params"]["vg_alpha"]
        self.soil_params['theta_s'] = dat["soil_params"]["theta_s"]
        self.soil_params['theta_r'] = dat["soil_params"]["theta_r"]
        self.soil_params['psi_top'] = dat["soil_params"]["psi_top"]
        self.soil_params['smcmax'] = dat["soil_params"]["smcmax"]

        self.NL = len(self.Z)
        #Volumetic heat capacity [j/m3/k]
        self.CWAT   = 4.188E06 # water heat capacity
        self.CICE   = 2.094E06 # ice heat capacity
        self.CAIR  = 1004.64  # air heat capacity
        self.CSOIL = 2.00E+6 # rock/soil heat capacity
        self.TKICE  = 2.2  # thermal conductiviyt of ice
        self.HFUS   = 0.3336E06 #latent heat of fusion (j/kg)
        self.PSISAT = 0.759 #Saturated matrix potential for soil type = silt loam
        self.GRAV = 9.86
        self.THKW = 0.57 # TC of water
        self.THKQTZ = 7.7 # TC of Quartz
        self.QUARTZ = 0.6 # loamy sand
        self.THKO  =  2.0 #TC of other mineral 
        self.TFRZ   = 273.15    #freezing/melting point (k)
        self.SMCMAX = 0.4 #porosity (maximum soil moisture)
        self.BEXP = 2.9 # Clap-Honnberger parameter
        self.WDEN = 1000. #[kg/m3]
        
        self.HC = np.zeros(self.NL) #volumetic heat capacity [j/m3/k]
        self.SMCT = np.zeros(self.NL)
        self.SMCLiq = np.zeros(self.NL)
        self.ST = np.array(dat["soil_temperature"]) # initial soil temperature
        self.Ttop = dat["surface_temperature"]

        # Time step and boundary condition
        
        self.bc_change = True

        if self.real_forcing == True:
            self.read_forcing()
            self.BC_T = self.forcing_data['TMP_2maboveground']
            self.Time = self.get_time_in_seconds()
        else:
            self.endtime = dat["end_time_d"] * 86400. #days
            self.dt = dat["dt_s"] # seconds
            self.NTsteps = int((self.endtime)/self.dt)
            self.BC_T = np.ones(self.NTsteps)*260. #default
            self.Time = np.array([self.starttime + i*self.dt for i in range(self.NTsteps+1)])
            #assert (len(self.BC_T) > 0)
            
    def get_Z(self):
        return self.Z
    
    def get_dz(self):
        tmp = np.zeros(self.NL)
        tmp[0] = self.Z[0]
        for i in range(self.NL-1):
            tmp[i+1] = self.Z[i+1] - self.Z[i]
        return tmp

    def get_time_in_seconds(self):
        st0 = pd.Timestamp(self.forcing_data['time'][0])
        Time = []
        for t in self.forcing_data['time']:
            t0 = (pd.Timestamp(t) - t0)
            t3 = pd.Timedelta(t0).total_seconds()
            Time.append(t3)
            
        self.starttime = Time[0]
        self.endtime = Time[-1]
        self.NTsteps = len(Time)
        return np.array(Time)
        
    def get_NTsteps(self):
        return self.NTsteps
    
    def read_forcing(self):
        self.forcing_data = pd.read_csv(self.forcing_file)

    def set_surface_bc(self,bc):
        self.BC_T = bc
    
    # get soil water content
    def get_theta(self):
        n = self.soil_params["vg_n"]
        m = self.soil_params["vg_m"]
        alpha = self.soil_params["vg_alpha"]
        theta_temp = np.zeros(self.NL)
        for i, z in enumerate(self.Z):
            p1 = 101325 + 1000 * 9.8 * z
            if i ==0:
                psi_new = self.soil_params["psi_top"]
            else:
                psi_new = self.soil_params["psi_top"] - p1/(1000 * 9.8)
            fact = 1 + (alpha * psi_new)**n
            theta_temp[i] = self.soil_params["theta_r"] + (self.soil_params["theta_s"] - self.soil_params["theta_r"]) * (fact**(-m))
            
        return np.round(theta_temp,3)



    def soil_heat_capacity(self):
        SMCMAX = self.soil_params["smcmax"]
        for i in range(self.NL):
            sice = self.SMCT[i] - self.SMCLiq[i]
            self.HC[i] = self.SMCLiq[i]*self.CWAT + sice*self.CICE + (1.0-SMCMAX)*self.CSOIL + (SMCMAX-self.SMCT[i])*self.CAIR


    #compute soil bulk thermal properties
    def thermal_conductivity(self):
        
        DF = np.zeros(self.NL)
    
        for i in range(self.NL):
            DF[i] = self.compute_TC(self.SMCT[i], self.SMCLiq[i])
            
        return DF

    ### Thermal Conductivity (Peters-Lidard Model)
    def compute_TC(self, SMC_c,SH2O_c): 

        SATRATIO = SMC_c / self.soil_params['smcmax']
    
        #TC of solids Eq. (10) Peters-Lidard
        THKS = (self.THKQTZ ** self.QUARTZ) * (self.THKO ** (1. - self.QUARTZ))
    
        #SATURATED THERMAL CONDUCTIVITY
        
        #UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
        XUNFROZ = SH2O_c / SMC_c # (phi * Sliq) / (phi * sliq + phi * sice) = sliq/(sliq+sice) 
        
        XU = XUNFROZ * self.SMCMAX # unfrozen volume fraction
        THKSAT = (THKS ** (1. - self.SMCMAX)) * (self.TKICE ** (self.SMCMAX - XU)) * (self.THKW ** XU )
        
        #DRY THERMAL CONDUCTIVITY
        GAMMD = (1. - self.SMCMAX)*2700. #dry density
        THKDRY = (0.135* GAMMD+ 64.7)/ (2700. - 0.947* GAMMD)
    
        # Kersten Number
        # for frozen soil
        if (SH2O_c + 0.001 > SMC_c ):
            KN = SATRATIO
        else:
            if SATRATIO > 0.1:
                KN = np.log10(SATRATIO) + 1.0
            else:
                KN = 0.0
    
        #Thermal conductivity
        DF = KN * (THKSAT - THKDRY) + THKDRY
    
        return DF




    # ### Freezing-point depression Eq. and mass/temprature correction due to phase change

    #REAL, INTENT(IN) :: SNEQVO    !snow mass at last time step(mm)
    #REAL, INTENT(IN) :: SNEQV     !snow water per unit ground area (mm)
    #SnH : snow depth
    
    #SnIce : ice in snow layer [mm] 
    #SnLiq : liquid in snow layer [mm] 
    
    #SMC : total soil water [m3/m3]
    #SH2O : soil liquid water [m3/m3] (liquid saturation??)
    
    #SMCMAX : porosity, saturated value of soil moisture (volumetric)
    
    #Freezing-point depression Eq. and mass/temprature correction due to phase change
    def PhaseChange(self):
        
        SUPERCOOL = np.zeros(self.NL) #supercooled water in soil
        MICE_L = np.zeros(self.NL) #snow/soil ice mass [mm]
        MLIQ_L = np.zeros(self.NL) #snow/soil liquid mass [mm]
        HM_L = np.zeros(self.NL) #energy residual [w/m2]
        XM_L = np.zeros(self.NL) #melting or freezing water [kg/m2]
        IMELT = np.zeros(self.NL)
        SNEQV = 0
        #compute mass of liquid/ice in snow/soil layers in mm
        for i in range(self.NL):
            if i < self.NSNOW: #snow layer
                MICE_L[i] = SNICE[i]
                MLIQ_L[i] = SNLIQ[i]
            else:
                # MICE and MLIQ are in units of [kg/m2]
                MICE_L[i] = (self.SMCT[i] - self.SMCLiq[i]) * self.dz[i] * self.WDEN # [kg/m2]
                MLIQ_L[i] = self.SMCLiq[i] * self.dz[i] * self.WDEN
    
        #set local variables
    
        #create copies of Mice and MLiq
        MICE_IN = np.array(MICE_L)
        MLIQ_IN = np.array(MLIQ_L)
        ##Phase change between ice and liquid water
    
        TMASS_L = np.zeros(self.NL) #snow/soil total (liquid + ice) mass [mm] (WMASS0)
        SUPERCOOL = np.zeros(self.NL)
        for i in range(self.NL):
            IMELT[i] = 0
            TMASS_L[i] = MICE_L[i] + MLIQ_L[i]

        # Soil water potential
        # SUPERCOOL is the maximum liquid water that can exist below (T - TFRZ) freezing point
        for i in range(self.NSNOW,self.NL):
            if (self.ST[i] < self.TFRZ):
                SMP = self.HFUS /(self.GRAV*self.ST[i]) * (self.TFRZ - self.ST[i] )       # [m] Soil Matrix potential
                SUPERCOOL[i] = self.SMCMAX*(SMP/self.PSISAT)**(-1./self.BEXP) #SMCMAX = porsity
                SUPERCOOL[i] = SUPERCOOL[i]*self.dz[i]* self.WDEN #[kg/m2]

        # ****** get layer freezing/melting index ************
        for i in range(self.NL):
            if MICE_L[i] >0 and self.ST[i] > self.TFRZ: #Melting
                IMELT[i] = 1
            elif MLIQ_L[i] > SUPERCOOL[i] and self.ST[i] <= self.TFRZ: #freezing
                IMELT[i] = 2

            # If snow exists, but its thickness is not enough to create a layer
            if (self.ISNOW == 0 and SNEQV >0):
                IMELT[0] = 1

        # ****** get excess or deficit of energy during phase change (use Hm) ********
        #HCPCT = volumetic heat capacity [J/m3/K]
        # Hm = (T- Tref) * HC * DZ /Dt = K * J/(m3 * K) * m * 1/s = (J/s)*m/m3 = W/m2
        self.soil_heat_capacity()
    
        #if HM < 0 --> freezing energy otherwise melting energy
        for i in range(self.NL):
            if IMELT[i] > 0:
                HM_L[i] = (self.ST[i] - self.TFRZ) * (self.HC[i] * self.dz[i])/self.dt
                self.ST[i] = self.TFRZ # Note the temperature does not go below 0 until the entire there is mixture of water and ice
            
            if (IMELT[i] == 1 and HM_L[i] <0):
                HM_L[i] = 0
                IMELT[i] = 0
            
            if (IMELT[i] == 2 and HM_L[i] > 0):
                HM_L[i] = 0
                IMELT[i] = 0
            
            XM_L[i] = HM_L[i]*self.dt/self.HFUS # melting or freezing water [kg/m2] (how much water needs to be melted or freezed for the given energy change)            
        #The rate of melting and freezing for snow without a layer, needs more work.
        if (self.ISNOW == 0 and SNEQV > 0. and XM_L[1] > 0.):
            print ('No snow layer-- TODO? maybe not..')
        
    
        #The rate of melting and freezing for snow and soil
        # here is the mass partition between ice and water and the corresponding adjustment for the next timestep
    
        for i in range(self.NL):
            if IMELT[i] >0 and abs(HM_L[i]) >0:
                if XM_L[i] >0: #melting
                    MICE_L[i] = max(0., MICE_IN[i]-XM_L[i])
                elif XM_L[i] <0: #freezing
                    if i < 0: #snow layers
                        MICE_L[i] = min(TMASS_L[i], MICE_L[i] - XM_L[i])
                    else: #soil layers
                        if TMASS_L[i] < SUPERCOOL[i] :
                            MICE_L[i] = 0
                        else:
                            MICE_L[i] = min(TMASS_L[i] - SUPERCOOL[i], MICE_IN[i] - XM_L[i])
                            MICE_L[i] = max(MICE_L[i],0)
                        
                HEATR = HM_L[i] - self.HFUS*(MICE_IN[i]-MICE_L[i])/self.dt #[W/m2] Energy Residual, last part is the energy due to change in ice mass
                MLIQ_L[i] = max(0.,TMASS_L[i]-MICE_L[i])

                #Correct the Temperature
                if abs(HEATR)>0:
                    f = self.dt/(self.HC[i] * self.dz[i]) # [m2 K/W]
                    self.ST[i] = self.ST[i] + f*HEATR # [K] , this is computed from HM = (T_n+1-T_n) * Heat_capacity * DZ/ DT
                    if i <-3: #snow layer
                        if (MLIQ_L[i]*MICE_L[i] >0): #if snow layer exits 
                            self.ST[i] = self.TFRZ
        
        for i in range(self.NSNOW,self.NL): #soil
            self.SMCLiq[i] =  MLIQ_L[i] / (self.WDEN * self.dz[i]) #[-]
            self.SMCT[i]  = (MLIQ_L[i] + MICE_L[i]) / (self.WDEN * self.dz[i]) # [-]


    def get_solution(self):
        return self.A

    ## Call Advance to advance timestep through while loop
    def run_model(self):
        
        Told = 0
        Tnew = 0
        self.A = []
        cycles = 0

        # write initial coditions to matrix A
        
        self.A.append([Told, copy.deepcopy(self.ST), copy.deepcopy(self.SMCLiq), self.SMCT])

        assert (self.real_forcing != self.bc_change)
        
        while (Tnew <= self.endtime):
            Tnew = Told + self.dt

            # set boundary conditions
            if self.bc_change == True and cycles > 3600:
                HE.set_boundary_conditions(285.)
            if self.real_forcing == True:
                HE.set_boundary_conditions(self.BC_T[cycles])

            
            self.ST = np.array(HE.AdvanceT(Tnew,Told))
            
            self.PhaseChange()
            DF = self.thermal_conductivity()
            HE.set_thermal_conductivity(DF)
            HE.set_temperature(self.ST)
            
            cycles = cycles + 1
            Told = Tnew

            self.A.append([Told, copy.copy(self.ST), copy.deepcopy(self.SMCLiq), copy.deepcopy(self.SMCT)])
        

    ## Call Advance to advance timestep through while loop
    def run_model_forcing(self):
        
        Told = self.Time[0]
        Tnew = self.Time[0]
        self.A = []
        cycles = 0

        # write initial coditions to matrix A
        
        self.A.append([Told, copy.deepcopy(self.ST), copy.deepcopy(self.SMCLiq), self.SMCT])

        assert (self.real_forcing != self.bc_change)

        for i in range(1,len(self.Time)):
            Tnew = self.Time[i]
            self.dt = (self.Time[i] - self.Time[i-1])#.seconds

            # set boundary conditions
            #if self.bc_change == True and cycles > 3600:
            #    HE.set_boundary_conditions(285.)
            #if self.real_forcing == True:
            HE.set_boundary_conditions(self.BC_T[cycles])

            self.ST = np.array(HE.AdvanceT(self.dt))
            
            self.PhaseChange()
            DF = self.thermal_conductivity()
            HE.set_thermal_conductivity(DF)
            HE.set_temperature(self.ST)
            
            cycles = cycles + 1
            Told = Tnew

            self.A.append([Told, copy.copy(self.ST), copy.deepcopy(self.SMCLiq), copy.deepcopy(self.SMCT)])
# Comparison with FrozenSoilModule running through BMI

def read_BMI_data(bmi_outfile, nz):

    BMI_ST = []
    BMI_Water =[]
    BMI_Time = []
    BMI_Ice = []
    nz = 8
    count = 0
    with open(bmi_outfile,'r') as file:
        for f in file.readlines()[9:]:
            sp = f.split()
            if "Time" in f:
                BMI_Time.append(float(sp[2]))
            elif ('Final' not in f and len(sp) == 3):
                BMI_ST.append(float(sp[0][:-1]))
                BMI_Water.append(float(sp[1][:-1]))
                BMI_Ice.append(float(sp[2]))
                count = count + 1


    bmiST = np.reshape(np.array(BMI_ST), (-1,nz))
    bmiW = np.reshape(np.array(BMI_Water), (-1,nz))
    bmiIce = np.reshape(np.array(BMI_Ice), (-1,nz))

    return bmiST, bmiW, bmiIce


def compare_with_bmi_solution(bmi_datafile, outfig_path, nz, Tsteps, step=12):

    bmiST, bmiW, bmiIce = read_BMI_data(bmi_datafile, nz)
    tmpdir = os.path.join(outfig_path, "tempdir")
    os.mkdir(tmpdir)
    
    # cycles = number of time steps
    for i in range(0,Tsteps+1,step):
        fig, axs = plt.subplots(1,2,figsize=(10,6))
        axs[0].plot(bmiST[i],Z,color='r',label='Soil T (BMI)')
        axs[0].plot(A[i][1],Z,color='k',linestyle='dashed',label='Soil T (Python)')
        axs[0].vlines(bmiST[0],ymin=0,ymax=7,color='grey',label='Initial ST',linestyle='dotted')
        axs[0].set_xlabel('Soil temperature [K]',fontsize=12)
        
        axs[1].plot(bmiW[i],Z,color='b',label='Liquid')
        axs[1].plot(bmiIce[i],Z,color='g',label='Ice')
        axs[1].plot(A[i][2],Z,color='k',linestyle='dashed',label='Liquid (Python)')
        axs[1].plot(bmiW[0],Z,color='grey',label='Initial SMC',linestyle='dotted')
        
        axs[1].set_xlabel('Soil moisture [K]',fontsize=12)
        axs[0].text(262.5, 5., 'Time = %s [d]'%(BMI_Time[i]/(3600.*24)), fontsize=12, color="green")
        
        axs[0].set_xlim(262,285.1)
        axs[1].set_xlim(-0.01,0.5)
    
        for j in range(2):
            axs[j].set_ylim(0.,6.9)
            axs[j].invert_yaxis()
            axs[j].xaxis.set_tick_params(labelsize=12)
            axs[j].yaxis.set_tick_params(labelsize=12)
            axs[0].set_ylabel('Depth [m]',fontsize=12)
            axs[0].legend(loc = 'lower left',labelspacing=.2,fontsize=10)
            axs[1].legend(loc = 'lower center',labelspacing=.1,fontsize=10)
        plt.savefig(outfig + 'fig-%s.png'%(BMI_Time[i]/3600.),dpi=100)
        plt.close()

    make_movie(tmpdir)


def sorted_nicely( l ):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

def make_movie(outfile_path):
    images = []
    dir_files = os.listdir(outfile_path)

    Files = sorted_nicely(dir_files)
    Files = [f for f in Files if f.endswith('.png')]
    print (len(dir_files), Files[:3])
    
    for filename in Files:
        images.append(imageio.imread(outfile_path+filename))
    imageio.mimsave(outfile_path+'movie.gif', images)



"""
# domain information such as number of cells(layers) in a column

#initial heat capacity and thermal conductivities
SNICE = [0,0,0] 
SNLIQ = [0, 0, 0]
SNEQV = 0
SNH = 0
QMELT = 0 # [mm/s]

PONDING = 0 #[mm]



#***************************************************************
# **** Initial Soil Temperature, Thermal Conductivity **********
# **** Set boundary conditions ******





# Run model till the end time
run_model(Tend, DT, BC_T, Tsteps=-1)

#SUPERCOOL gives the amount of liquid water that must exist at the given temperature

outfig_path = '/Users/ahmadjan/Core/SimulationsData/postprocessing/frozensoil/'
bmi_datafile='/Users/ahmadjan/codes/bmi/bmi-frozensoil/_build/testing/bmifrozensoilcxx.out'

#compare_with_bmi_solution(bmi_datafile, outfig_path, Nz, Tsteps, step=12)
"""
