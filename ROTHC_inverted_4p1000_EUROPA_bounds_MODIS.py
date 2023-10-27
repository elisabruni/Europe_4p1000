#Matricial form
import sys
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import least_squares
#import numdifftools as ndt
import math
from random import gauss
import xlrd
import pandas as pd
import time
import datetime
from datetime import datetime
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import seaborn as sns

#################
#Type of simulation
##################
optimization_4p1000=True
optimization_4p1000_CALIB=True
print_mode=False
save_mode = True

nppinput_from = '_MODIS'
#nppinput_from = '_ISIMIP'
#nppinput_from = '_DAYCENT'

tstart = time.time()
#######################

###################
#Data

ROOTDIR='/SET/DIRECTORY/'
loc_exp = ROOTDIR+'LUCAS_2009_2015_2018_crop_wetness.csv'
#loc_hum = ROOTDIR+'_lcs_AR_GR30cm_v3.1.xls'
OUTPUT_files=ROOTDIR+'OUTPUTS/ROTHC/'

if(save_mode):
    out_SOC_NCAL = open(OUTPUT_files+"SOCdyn_NCAL_ROTHC_2015_2099_bounds"+nppinput_from+".txt","wb")
    out_SOC_CAL = open(OUTPUT_files+"SOCdyn_CAL_ROTHC_2015_2099_bounds"+nppinput_from+".txt","wb")

##################
#Functions
##################
##############################
#Convert daily data to monthly
##############################
def DailyToMonthly(filein,typeConv):
		c=0
		x=0
		var_month=np.zeros(np.int(len(filein)/365.*12.))
		month_days = [31,28,31,30,31,30,31,31,30,31,30,31]
		lent = np.int(len(var_month))
		for j in range(0,lent):
			for i in month_days:
				if(typeConv==0.): #for precipitation and PET (cumulative)
					var_month[c]=np.sum(filein[x:x+i])
				else: #for Temperature (mean)
					var_month[c]=np.sum(filein[x:x+i])/np.float(i)
				x=x+i
				c+=1
			if c>len(var_month)-12:
				break
		return var_month
	

#########################################################
def DailyToYearly(filein,typeConv):
#########################################################
                c=0
                x=0
        #Set number of years
                var_year=np.zeros(np.int(len(filein)/365))
                lent = np.int(len(var_year))
                for j in range(0,lent):
                        if(typeConv==0.): #for precipitation and PET (cumulative)
                                var_year[c]=np.sum(filein[x:365+x])
                        else: #for Temperature (mean)
                                var_year[c]=np.sum(filein[x:365+x])/365
                        x=x+365
                        c+=1
                return var_year

################################
################################
#ENVIRONMENTAL FUNCTIONS
################################
################################

#############################################
def control_temp_func(temp_mean,param_temp):
#############################################
    #temp mean = monthly average temperature (C)

    #a = 47.91/(1.+np.exp(106.06/(temp_mean+18.27)))
    a = 47.91/(1.+np.exp(106.06/(temp_mean+param_temp)))
    #print 'contro temp',a
    #Test6
    a=np.minimum(a,4.)
    #Test12
    #a=np.minimum(a,4.5)
    return a

#############################################
def control_moist_func(clay, rain, mean_pot_evapot,c_param):
#############################################
	#clay %, rain (mm), mean_pot_evapot (mm)
	#global n_an,n_month,accTSMD, maxTSMD
    global accTSMD, maxTSMD
    
	#Maximum Topsoil Moisture Deficit (TSMD)
    #For a soil depth 0-23cm
    #If soil depth is for ex  X(cm) -> maxTSMD=maxTSMD*X/23
    #For vegetated soils
    if(c_param <1.):
        maxTSMD = -(20.+1.3*clay-0.01*(clay)**2) 
    #For bare soil divide by 1.8
    else:
        maxTSMD = -(20.+1.3*clay-0.01*(clay)**2)/1.8

	#Accumulated TSMD
	#Attenzione formula originale: invece di mean_pot_evapot c'e (open_pan_evapot*0.75). Siccome uso mean_pot_evapot (from Muller o calculated), open_pan_evapot = mean_pot_evapot/0.75

	#number_of_months = n_month

    pot_ev_month = mean_pot_evapot
    rain_month = rain
    if(pot_ev_month < rain_month and accTSMD==0. ):
            accTSMD = 0.
    else:
            accTSMD = accTSMD+(rain_month-pot_ev_month)
            accTSMD = np.maximum(maxTSMD, accTSMD)
            accTSMD = np.minimum(accTSMD,0.)

    #Moisture rate modifier (b)
    if(np.abs(accTSMD)<0.444*np.abs(maxTSMD)):
            b=1.
    else:
            b = 0.2+(1.-0.2)*(np.abs(maxTSMD)-np.abs(accTSMD))/(np.abs(maxTSMD)-0.444*np.abs(maxTSMD))
            #b = np.minimum(1.,b)

    #print 'control moist',b
    #print 'rain',rain
    #print 'mean_pot_evapot',mean_pot_evapot
    #print 'max',maxTSMD
    #print 'acc',accTSMD
    #print 'denom',(maxTSMD-0.444*maxTSMD)

    return b,np.abs(accTSMD)


############################################################
#BIO_HUM partition
############################################################

def BIO_HUM_partition(clay):
	#Partitioning of C between CO2 and BIO+HUM
	X_ratio = 1.67*(1.85+1.6*np.exp(-0.0786*clay))
	BIO_HUM_frac = 1./(X_ratio+1.)
	return BIO_HUM_frac

############################################################
#IOM partition
############################################################

def IOM_partition(SOC):
	#tC/ha
	IOM_frac = 0.049*(SOC**(1.139))
	return IOM_frac

############################################################
#A_matrix
############################################################
def A_matrix(clay,plantin,EOMin,DPM_RPM_HUM_frac_inputs,DPM_RPM_HUM_frac_inputs_EOM):
	
	alpha = BIO_frac_fromRPM*BIO_HUM_partition(clay)
	beta = HUM_frac_fromRPM*BIO_HUM_partition(clay)+plantin[2]*DPM_RPM_HUM_frac_inputs[2]+EOMin[2]*DPM_RPM_HUM_frac_inputs_EOM[2]

	A_matrix = np.zeros((4,4))
	np.fill_diagonal(A_matrix, -1)	

	A_matrix[2,0]=alpha
	A_matrix[2,1]=alpha
	A_matrix[2,2]=(alpha-1.)
	A_matrix[2,3]=alpha

	A_matrix[3,0]=beta
        A_matrix[3,1]=beta
        A_matrix[3,2]=beta
        A_matrix[3,3]=(beta-1.)

	A_out = A_matrix
	
	return A_out

############################################################
#K_matrix
############################################################

def K_matrix(clay,temp_mean,rain, mean_pot_evapot,c_param,t_param,param_temp_in):
	global moist_rate_b, accTSMD_abs, temp_rate_a
	K_matrix=np.zeros((4,4))
	moist_rate_b, accTSMD_abs = control_moist_func(clay, rain, mean_pot_evapot,c_param)
	temp_rate_a = control_temp_func(temp_mean,param_temp_in)
	K_matrix[0,0]=kDPM*temp_rate_a*moist_rate_b*c_param*t_param
	K_matrix[1,1]=kRPM*temp_rate_a*moist_rate_b*c_param*t_param
	K_matrix[2,2]=kBIO*temp_rate_a*moist_rate_b*c_param*t_param
	K_matrix[3,3]=kHUM*temp_rate_a*moist_rate_b*c_param*t_param
	K_out = K_matrix
	return K_out

############################################################
#Spinup
############################################################
def spinup(SOC_data_init,plantin,EOMin,clay,temp,rain, pot_evapot,c_param,t_param,param_temp_in,DPM_RPM_HUM_frac_inputs,DPM_RPM_HUM_frac_inputs_EOM):
    #global ABOVE_mean,BELOW_mean, err2_above, err2_below
    global accTSMD
    global b_rate_array, accTSMD_abs_array
    global a_rate_array, temp_mean_array

    litterin=plantin+EOMin

    matrix_in_mean=np.append(litterin,[0.])

    #Initialize accTSMD
    accTSMD =0.
    b_rate_array=np.zeros(n_month_ss)
    accTSMD_abs_array=np.zeros(n_month_ss)
    a_rate_array=np.zeros(n_month_ss)
    temp_mean_array=np.zeros(n_month_ss)
    for ts in range(0,n_month_ss):
        temp_mean = temp[ts]
        rain_mean = rain[ts]
        pot_evapot_mean = pot_evapot[ts]
	#Calculate mean K_matrix
        if (ts == 0):
            kk_ma_mean=K_matrix(clay,temp_mean,rain_mean, pot_evapot_mean,c_param,t_param,param_temp_in)
        else:
            kk_ma_mean+=K_matrix(clay,temp_mean,rain_mean, pot_evapot_mean,c_param,t_param,param_temp_in)
	b_rate_array[ts]=moist_rate_b
	accTSMD_abs_array[ts]=accTSMD_abs
	a_rate_array[ts]=temp_rate_a
	temp_mean_array[ts]=temp_mean

    kk_ma_mean=kk_ma_mean/n_month_ss

    #print 'litterin'
    #print litterin
    #print 'DPM_RPM_HUM_frac_inputs'
    #print DPM_RPM_HUM_frac_inputs

    a_ma_mean=A_matrix(clay,plantin,EOMin,DPM_RPM_HUM_frac_inputs,DPM_RPM_HUM_frac_inputs_EOM)
    ss_spinup=-np.linalg.solve(np.dot(a_ma_mean,kk_ma_mean),matrix_in_mean)
    IOM_pool = IOM_partition(SOC_data_init)
    ss_spinup_IOM = np.append(ss_spinup,IOM_pool)

    #print 'IN'
    #print matrix_in_mean
    #print 'A'
    #print a_ma_mean
    #print 'K'
    #print kk_ma_mean


    return ss_spinup_IOM


############################################################
#Forward
############################################################
def forward_ACTIVE(SOC_data_init,n_an,init,plantin,EOMin,clay,temp,water,mean_pot_evapot,c_param,t_param,param_temp_in,DPM_RPM_HUM_frac_inputs,DPM_RPM_HUM_frac_inputs_EOM):
    global temp_mean, rain_mean,accTSMD # needed!!
    global dt,n_month,one_month # not really needed

    litterin=plantin+EOMin

    matrix_cpools_tmean = np.zeros((n_an,5))

    length=n_month
    matrix_cpools_t=np.zeros((length+1,4))
    matrix_in = np.zeros(4)
    matrix_cpools_t[0]=init[0:4]

    for i in range(0,len(litterin)):
        matrix_in[i]=litterin[i]

    a_ma=A_matrix(clay,plantin,EOMin,DPM_RPM_HUM_frac_inputs,DPM_RPM_HUM_frac_inputs_EOM)

    #Initialize accTSMD 
    accTSMD=0.
    #Annual average
    for x in range(0,n_an):
        matrix_cpools_ymean = np.zeros(4)
	#Monthly Euler
        for ts in range(x*one_month,(one_month*x)+one_month):
            temp_mean = temp[ts]
            rain_mean = water[ts]
	    evapot_mean = mean_pot_evapot[ts]
            matrix_current= matrix_cpools_t[ts]
            kk_ma = K_matrix(clay,temp_mean,rain_mean,evapot_mean,c_param,t_param,param_temp_in)
            matrix_next = matrix_current + matrix_in + np.dot(a_ma,np.dot(kk_ma,matrix_current))*dt
            matrix_cpools_t[ts+1]=matrix_next
            matrix_cpools_ymean += matrix_next

	#Yearly average C pools		
        matrix_cpools_ymean = matrix_cpools_ymean/one_month

	IOM_pool = IOM_partition(SOC_data_init)
        matrix_cpools_tmean[x] = np.append(matrix_cpools_ymean,IOM_pool)

    return matrix_cpools_tmean

############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):
    global NEW_ITER
    global n_an, c_vegetated, DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM
    global clay_2009,temp_in_monthly_fw,water_in_monthly_fw,PET_in_monthly_fw,c_vegetated,t_param
    global Corg_2009,Corg_in,SOC_lcs_2015
    global print_mode

    #Suppose Corg does not change
    if(Corg_2009<in_new): #if generated in_new is higher than Corg_2009 -> Corg is kept the same
        Corg_new = Corg_2009
        EOM_new = Corg_in     
    else:     #Otherwise in_new is plant only  
        Corg_new = 0
        EOM_new = np.array([0.,0.,0.])

    #plant new is the difference between in_new and Corg
    plant_tot_new = in_new-Corg_new

    #Suddividi C from plant nei different pools
    plant_new = np.zeros(3)
    plant_new[0] = DPM_RPM_HUM_frac[0]*plant_tot_new
    plant_new[1] = DPM_RPM_HUM_frac[1]*plant_tot_new
    plant_new[2] = DPM_RPM_HUM_frac[2]*plant_tot_new    

    #print ' '
    #print '######################'
    #print 'in_new',in_new
    #print 'plant_new',plant_new
    #print 'EOM_new',EOM_new
    #print DPM_RPM_HUM_frac_FYM
    #print ' '
    #print 'DPM_RPM_HUM_frac_mix'
    #print DPM_RPM_HUM_frac_mix
    #print '####################'
    #print ' '
    
    #4 per 1000 between 2015-2099
    if(opt_NONCAL): #Non calibrated
        global spinup_c, predict_c
        predict_c_pools = forward_ACTIVE(SOC_lcs_2015,n_an,spinup_c,plant_new,EOM_new,clay_2009,temp_in_monthly_fw,water_in_monthly_fw,PET_in_monthly_fw,c_vegetated,t_param,non_optimized_t_param,DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM)
        predict_c = np.sum(predict_c_pools,axis=1)
        J_new = abs(target*n_an - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

    elif(opt_CAL): #Calibrated
        global spinup_c_calib, predict_c_calib
        predict_c_pools = forward_ACTIVE(SOC_lcs_2015,n_an,spinup_c_calib,plant_new,EOM_new,clay_2009,temp_in_monthly_fw,water_in_monthly_fw,PET_in_monthly_fw,c_vegetated,t_param,estimate_T_param_calib,DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM)
        predict_c_calib = np.sum(predict_c_pools,axis=1)
        J_new = abs(target*n_an - np.sum(predict_c_calib[predict_c_calib.shape[0]-1]-predict_c_calib[0]))


    #print 'OBJ FUN', J_new
    #print 'predict_c',predict_c
    #print 'predict_c_pools', predict_c_pools

    NEW_ITER+=1
    if (print_mode):
        if ( NEW_ITER % 100 == 0 ):
             print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
    param_opt=in_new
    return J_new


######################################
#Set bounds and constraints for the optimization
#########
bnds=[(0,30)]

#a_constr=1
#b_constr=1
#ab_ratio_var=(1-a_constr)*100
#con1={'type':'ineq','fun':constr1}
#con2={'type':'ineq','fun':constr2}
#cons=[con1,con2]

############
#Parameters

CHI2_PRINT_FREQUENCY=50

#Decomposition rate

kDPM=10.
kRPM=0.3
kBIO=0.66
kHUM=0.02


#Since k is vased on yearly decomposition rates
t_param=1./12. #year->month

#Soil cover
c_vegetated = 0.6
c_bare = 1.

#Pools partitioning for plant residues
#For agricultural crops and improved grassland
DPM_frac_plant = 0.59
RPM_frac_plant = 1-DPM_frac_plant
HUM_frac_plant = 0.

DPM_RPM_HUM_frac = np.array([DPM_frac_plant,RPM_frac_plant,HUM_frac_plant])

#Pools partitioning for FYM
DPM_frac_FYM = 0.49
RPM_frac_FYM = 0.49
HUM_frac_FYM = 0.02

DPM_RPM_HUM_frac_FYM = np.array([DPM_frac_FYM,RPM_frac_FYM,HUM_frac_FYM])

#Partitioning between HUM and BIO from RPM
BIO_frac_fromRPM = 0.46
HUM_frac_fromRPM = 1.-BIO_frac_fromRPM

#Initialization
NEW_ITER = 0
one_month = 12
dt=1.

DPM_RPM_HUM_frac_mix = np.array([0.,0.,0.])

non_optimized_t_param = 18.27

######################################################
#For each site: Set SITE name and experiment duration
#####################################################

C_input_exp = pd.read_csv(loc_exp)
#df_humcoef = pd.read_excel(loc_hum, 'Sheet2',usecols="B:F",nrows=37)

loc_hanpp = ROOTDIR+'NPPin_lucas09.csv'

#HANPP
df_hanpp=pd.read_csv(loc_hanpp)
df_hanpp=df_hanpp[['sample_ID,N,10,0','MEAN,N,19,11']]
df_hanpp=df_hanpp.rename(columns={'MEAN,N,19,11':'MEAN_HANPP'})
#print df_hanpp

#MERGE HANPP AND LUCAS on sample_ID
C_input_exp = pd.merge(C_input_exp,df_hanpp, left_on='sample_ID,N,9,0',  right_on = 'sample_ID,N,10,0', left_index=True)

#Test
#C_input_exp = C_input_exp.loc[(C_input_exp['LC09,C,4'].str.contains('B') & C_input_exp['LC12,C,4'].str.contains('B') &
#    C_input_exp['LC15,C,4'].str.contains('B'))].head(2)

#>>>>>>>>All crop sites both in 2009 and 2015 <<<<<<<<<<<<<<
C_input_exp = C_input_exp.loc[(C_input_exp['LC09,C,4'].str.contains('B') & C_input_exp['LC12,C,4'].str.contains('B') &
        C_input_exp['LC15,C,4'].str.contains('B'))]

#Eliminate trees
list_trees = ['B55','B71','B73','B74','B75','B81','B82','B83','B84']
C_input_exp = C_input_exp[~C_input_exp['LC15,C,4'].isin(list_trees)]

#>>>>>>>> read LUCAS data <<<<<<<<<<<<<<

#Calcul stocks:
#OC=C_input_exp['OC,N,19,13']
#VGR=C_input_exp['VGR,N,19,15']
#BD=C_input_exp['BD,N,19,15']
#SOC = OC* (1-VGR)* BD * 30./10.

#Appends values of variable for each site
SITE_Cin_MOD = C_input_exp['NPPmod,N,19,14'].reset_index(drop=True) #average 2000-2009 NPP from MODIS 
SITE_Cin_hanpp = C_input_exp['MEAN_HANPP'].reset_index(drop=True) #fraction of npp that goes into the soil

SITE_clay_2009 = C_input_exp['clay,N,19,14'].reset_index(drop=True) #%
SITE_SOC_lcs_2009 = C_input_exp['SOC_lcs,N,19,13'].reset_index(drop=True) #OC stock 2009-12 from OC* (1-VRG)* BD  * cm  (if 30 cm, we assume a stratification only for grassland 20 cm *1.27)
SITE_SOC_lcs_sd_2009 = C_input_exp['SOCsd_lcs,N,19,13'].reset_index(drop=True) #confidence interval by Montecarlo
SITE_SOC_tr_2009 = C_input_exp['SOCtr,N,19,15'].reset_index(drop=True) #ratio passive:total SOC pools at 2009
SITE_Corg_2009 = C_input_exp['Corg,N,19,12'].reset_index(drop=True) #possible C applied with manure (kg/ha) Derived from overlay with large-scale database (no direct point information)

#SITE_SOC_lcs_2015 = C_input_exp['SOC_lcs15,N,19,12'].reset_index(drop=True) #SOC stocks 2015
SITE_SOC_lcs_2015 = C_input_exp['SOCstocks_2015'].reset_index(drop=True) #SOC stocks 2015

#Other variables for statistical calibration
#SITE_CaCO3_2009 = C_input_exp['CaCO3,N,19,13'].reset_index(drop=True) #CacO3 g/kg%
#SITE_N_2009 = C_input_exp['N,N,19,14'].reset_index(drop=True) #Nin
SITE_pH_2009 = C_input_exp['pH_in_CaCl,N,19,15'].reset_index(drop=True) #ph

#>>>>>>>>>>>>>#####<<<<<<<<<<<<<<

site_names_all = C_input_exp['POINT_ID,N,11,0'].unique()[0:len(C_input_exp)]
site_names_all = map(str, site_names_all)

N_sites=len(site_names_all)

#Simulations lenght
n_an_ss = 2014-2006+1 #yr ss
n_month_ss = n_an_ss*12 #months ss
n_an = 2099-2015+1 #yr fwd
n_month = n_an*12 #months fwd

#>>>>>>>> CLIMATE FORCING  <<<<<<<<<<<<<<
#print site_names_all
SITE_temp_in_monthly_ss = []
SITE_temp_in_monthly_fw = []
SITE_PET_in_monthly_ss = []
SITE_PET_in_monthly_fw = []
SITE_water_in_monthly_ss = []
SITE_water_in_monthly_fw = []

#For calibration
SITE_temp_in_yearly_fw = []
SITE_PET_in_yearly_fw = []
SITE_water_in_yearly_fw = []

for site in site_names_all:
    #print 'ID site: ',site

    #Soil temperature
    soil_temp =ROOTDIR+'ISIMIP_data/TEMP_RCP26/'+site+'_TAS_2006_2099.txt'
    temp_in_ar=[]
    prev_line=0
    with open(soil_temp) as fileID:
        for _ in range(7): #skip first 7 lines
            next(fileID)
        for line in fileID:
            #Get the info and put it separately in an array
            sub1=line.split(":",1)[1].strip() #Eliminate white spaces

            #MISSING POINTS
            if(sub1=='....'):
                sub1=prev_line
            #---------------

            temp_in = float(sub1)
            prev_line = temp_in

            temp_in_ar.append(temp_in)

    temp_in_ar = np.array(temp_in_ar)
    #K to C
    temp_in_monthly=DailyToMonthly(temp_in_ar,1)-273.15

    #For calibration
    temp_in_yearly=DailyToYearly(temp_in_ar,1)-273.15

    #Potential evapotranspiration
    soil_PET =ROOTDIR+'ISIMIP_data/POTEVAP_RCP26/'+site+'_POTEVAP_2006_2099.txt'
    PET_in_ar=[]
    prev_line=0
    with open(soil_PET) as fileID:
        for _ in range(7): #skip first 7 lines
            next(fileID)
        for line in fileID:
            #Get the info and put it separately in an array
            sub1=line.split(":",1)[1].strip() #Eliminate white spaces

            #MISSING POINTS
            if(sub1=='....'):
                sub1=prev_line
            #---------------

            PET_in = float(sub1)
            prev_line = PET_in

            PET_in_ar.append(PET_in)

    #kg/m2/sec (daily interp) to mm/month (1kg/m2 = 1mm)
    #PET_in_monthly=PET_in_ar*3600.*24.
    PET_in_ar = np.array(PET_in_ar)*3600.*24.
    PET_in_monthly=DailyToMonthly(PET_in_ar,0)

    #For calibration
    PET_in_yearly=DailyToYearly(PET_in_ar,0)


    #Rain
    soil_PR =ROOTDIR+'ISIMIP_data/PR_RCP26/'+site+'_PR_2006_2099.txt'
    PR_in_ar=[]
    prev_line=0
    with open(soil_PR) as fileID:
        for _ in range(7): #skip first 7 lines
            next(fileID)
        for line in fileID:
            #Get the info and put it separately in an array
            sub1=line.split(":",1)[1].strip() #Eliminate white spaces

            #MISSING POINTS
            if(sub1=='....'):
                sub1=prev_line
            #---------------

            PR_in = float(sub1)
            prev_line = PR_in

            PR_in_ar.append(PR_in)

    PR_in_ar = np.array(PR_in_ar)

    #Snowfall
    soil_PRSN =ROOTDIR+'ISIMIP_data/PRSN_RCP26/'+site+'_PRSN_2006_2099.txt'
    PRSN_in_ar=[]
    prev_line=0
    with open(soil_PRSN) as fileID:
        for _ in range(7): #skip first 7 lines
            next(fileID)
        for line in fileID:
            #Get the info and put it separately in an array
            sub1=line.split(":",1)[1].strip() #Eliminate white spaces

            #MISSING POINTS
            if(sub1=='....'):
                sub1=prev_line
            #---------------

            PRSN_in = float(sub1)
            prev_line = PRSN_in

            PRSN_in_ar.append(PRSN_in)

    PRSN_in_ar = np.array(PRSN_in_ar)

    #Rain + snow (kg/m2/sec (daily) -> mm/month)
    water_tot = (PR_in_ar+PRSN_in_ar)*3600.*24.
    water_in_monthly=DailyToMonthly(water_tot,0)

    #For calibration
    water_in_yearly = DailyToYearly(water_tot,0)

    #Calcola forcing for spin-up (2006-2014)
    temp_in_monthly_ss = temp_in_monthly[0:n_month_ss]
    PET_in_monthly_ss = PET_in_monthly[0:n_month_ss]
    water_in_monthly_ss = water_in_monthly[0:n_month_ss]

    #Calcola forcing for simulations (2015-2099)
    #temp_in_monthly_fw = temp_in_monthly[n_month_ss+1:]
    #PET_in_monthly_fw = PET_in_monthly[n_month_ss+1:]
    #water_in_monthly_fw = water_in_monthly[n_month_ss+1:]
    temp_in_monthly_fw = temp_in_monthly[n_month_ss:]
    PET_in_monthly_fw = PET_in_monthly[n_month_ss:]
    water_in_monthly_fw = water_in_monthly[n_month_ss:]

    #For calibration
    temp_in_yearly_fw = temp_in_yearly[n_an_ss:]
    PET_in_yearly_fw = PET_in_yearly[n_an_ss:]
    water_in_yearly_fw = water_in_yearly[n_an_ss:]    

    #Append forcint for all sites
    SITE_temp_in_monthly_ss.append(temp_in_monthly_ss)
    SITE_temp_in_monthly_fw.append(temp_in_monthly_fw)
    SITE_PET_in_monthly_ss.append(PET_in_monthly_ss)
    SITE_PET_in_monthly_fw.append(PET_in_monthly_fw)
    SITE_water_in_monthly_ss.append(water_in_monthly_ss)
    SITE_water_in_monthly_fw.append(water_in_monthly_fw)

    #For calibration
    SITE_temp_in_yearly_fw.append(temp_in_yearly_fw)
    SITE_PET_in_yearly_fw.append(PET_in_yearly_fw)
    SITE_water_in_yearly_fw.append(water_in_yearly_fw)


#>>>>>>>> END_OF_INITIALIZATION <<<<<<<<<<<<<<<<<<<<<<<<<<<<

tinit = time.time()


#>>>>>>>> MAIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<
SOC_init_model = np.zeros(N_sites)
SOC_init_model_calib = np.zeros(N_sites)
C_add_4p1000 = np.zeros(N_sites)
C_add_4p1000_calib = np.zeros(N_sites)
SOC_2018 = np.zeros(N_sites)
SOC_2018_calib = np.zeros(N_sites)
ACT_SOC_2015 = np.zeros(N_sites)
ACT_SOC_2015_calib = np.zeros(N_sites)
SITE_estimate_T_param_calib = np.zeros(N_sites)

for j in range(0,N_sites):
    print ' '
    print 'START'
    print '---'
    print 'Site: ',site_names_all[j]

    #Upload forcing 
    temp_in_monthly_ss = SITE_temp_in_monthly_ss[j]
    temp_in_monthly_fw = SITE_temp_in_monthly_fw[j]
    PET_in_monthly_ss = SITE_PET_in_monthly_ss[j]
    PET_in_monthly_fw = SITE_PET_in_monthly_fw[j]
    water_in_monthly_ss = SITE_water_in_monthly_ss[j]
    water_in_monthly_fw = SITE_water_in_monthly_fw[j]

    #Upload soil data
    Cin_MOD = SITE_Cin_MOD.iloc[j]/12. #MODIS tC/ha/yr to tC/ha/month
    npp_in = SITE_Cin_hanpp.iloc[j] #fraction
    Cin_MOD = Cin_MOD*npp_in #actual npp that goes into the soil

    Corg_2009 = SITE_Corg_2009.iloc[j]/1000./12. #kg/ha/yr to t/ha/month

    #Dividi Cin in different pools
    plant_in = np.array([Cin_MOD*DPM_frac_plant,Cin_MOD*RPM_frac_plant,Cin_MOD*HUM_frac_plant])
    Corg_in = np.array([Corg_2009*DPM_frac_FYM,Corg_2009*RPM_frac_FYM,Corg_2009*HUM_frac_FYM])

    Cin_array = plant_in+Corg_in # Array of C in total (plant + EOM)
    C_mean = np.sum(Cin_array) #C in total

    clay_2009 = SITE_clay_2009.iloc[j]
    SOC_lcs_2009 = SITE_SOC_lcs_2009.iloc[j]
    SOC_lcs_sd_2009 = SITE_SOC_lcs_sd_2009.iloc[j]
    SOC_tr_2009 = SITE_SOC_tr_2009.iloc[j]
    SOC_lcs_2015 = SITE_SOC_lcs_2015.iloc[j]


    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #Run spinup
    spinup_c=spinup(SOC_lcs_2015,plant_in,Corg_in,clay_2009,temp_in_monthly_ss,water_in_monthly_ss,
            PET_in_monthly_ss,c_vegetated,t_param,non_optimized_t_param,DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM)
    spinup_c_sum=np.sum(spinup_c) #tC/ha


    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<   
    #Run forward
    fwd_c=forward_ACTIVE(SOC_lcs_2015,n_an,spinup_c,plant_in, Corg_in,clay_2009,temp_in_monthly_fw,water_in_monthly_fw,PET_in_monthly_fw,c_vegetated,t_param,non_optimized_t_param,DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM)

    #Dynamics_all=np.concatenate((np.array([spinup_c]),fwd_c))
    Dynamics_sum = np.sum(fwd_c,axis=1)
    #print 'Dynamics non calib'
    #print Dynamics_sum

    #2015
    SITE_SOC_model = Dynamics_sum[0]
    SOC_init_model[j] = SITE_SOC_model

    #2015 active (year 0)
    ACT_SOC_2015[j] = fwd_c[0,2] #BIO pool

    #2018
    SOC_2018[j] = Dynamics_sum[3]


    #print '2015: '
    #print SITE_SOC_model
    #print 'SOC in 2009:'
    #print SOC_lcs_2009
    #print 'SOC in 2015:'
    #print SOC_lcs_2015
    #print '---'


    if(save_mode):
        np.save(out_SOC_NCAL,Dynamics_sum)

    #Plot SOC dynamics over 2015-20100
    #plt.plot(np.arange(n_an),Dynamics_sum)

    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #Run inverse modelling to 4p1000

    if(optimization_4p1000):
        opt_NONCAL=True
        opt_CAL=False

        target = SITE_SOC_model*0.004

        #litter income prior (1D) ->splitted in J_new
        in_opt = C_mean*(1+0.004)

        print ' '

        opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt = opt_mean.x[0]
        #print "SLSQP: Optimum solution:", litter_opt

        #calculate percentage increase/decrease of inputs
        delta_input = (litter_opt-C_mean)*12. #tC/ha/yr
        input_change=(litter_opt-C_mean)*100./C_mean

        C_add_4p1000[j]=delta_input

        if(print_mode):
            #litter_opt_array=np.array([litter_opt*DPM_RPM_HUM_frac_mix[0],litter_opt*DPM_RPM_HUM_frac_mix[1],litter_opt*DPM_RPM_HUM_frac_mix[2]])
            print ' '
            print 'NON CALIBRATED SIMULATIONS'
            print 'Initial litter (tC/ha/month):'
            print C_mean
            print '4p1000 litter (tC/ha/month):'
            print litter_opt

            print "% change of litter inputs:",input_change

            print ' '

            #Check 4p1000 target is reached
            END=predict_c.shape[0]-1
            #print 'END',END
            C_fin=np.sum(predict_c[END])
            #C_init=np.sum(predict_c[0])
            C_init=SITE_SOC_model
            SUMO=C_fin-C_init
            Target_reached=(C_fin-C_init)/(C_init*n_an)
            if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
            else:
                print "Target not reached"
                print 'C init', C_init
                #print predict_c[0]
                print SITE_SOC_model
                print 'C_fin', C_fin
            print ' '



    #>>>>>>>>>>>>>>>> 4p1000 with optimized param <<<<<<<<<<<<<

    if(optimization_4p1000_CALIB):
        opt_NONCAL=False
        opt_CAL=True

        #########################################################
        #Derive calibrated parameters from statistical functions
        #########################################################
        #CaCO3 = SITE_CaCO3_2009.iloc[j]
        #N_in = SITE_N_2009.iloc[j]
        ##print 'N_in',N_in
        #if N_in > 0:
        #    N_in_cat = 1
        #else:
        #    N_in_cat = 0

        pH_2009 = SITE_pH_2009.iloc[j]

        loc_functions = '/home/surface7/ebruni/MULTIMODEL_EUROPA/CALIB_PARAM_FUNC/V8_wetness2/'
        func_df = pd.read_csv(loc_functions+'T_param_calib_func.csv')
        interc_coef=func_df.loc['(Intercept)']
        SOC_init_coef=func_df.loc['Initial_SOC']
        temp_coef=func_df.loc['Temp']
        prec_coef=func_df.loc['Prec']
        PET_coef=func_df.loc['PET']
        Cin_coef=func_df.loc['Litter_in'] #litter in in tC/ha/yr
        #N_coef=func_df.loc['c(Nin)']
        #CaCO3_coef=func_df.loc['Carbonate']
        Clay_coef=func_df.loc['Clay']
        pH_coef = func_df.loc['pH']

        #Yearly climate var for calibration
        PET_in_yearly_fw = SITE_PET_in_yearly_fw[j]
        water_in_yearly_fw = SITE_water_in_yearly_fw[j]
        temp_in_yearly_fw = SITE_temp_in_yearly_fw[j]
        
        #V7
        #estimate_T_param_calib = interc_coef+SOC_init_coef*SOC_lcs_2015+temp_coef*np.mean(temp_in_yearly_fw)+prec_coef*np.mean(water_in_yearly_fw)+PET_coef*np.mean(PET_in_yearly_fw)+(Cin_MOD*12.)*Cin_coef+N_in_cat*N_coef+CaCO3*CaCO3_coef
        #V8
        estimate_T_param_calib = interc_coef+SOC_init_coef*SOC_lcs_2015+temp_coef*np.mean(temp_in_yearly_fw)+prec_coef*np.mean(water_in_yearly_fw)+PET_coef*np.mean(PET_in_yearly_fw)+(Cin_MOD*12.)*Cin_coef+clay_2009/100.*Clay_coef+pH_2009*pH_coef
        estimate_T_param_calib = np.float(estimate_T_param_calib)
        #Test 1
        #T_param bounded to be > - temp_min and Func temp bounded to be <4.5
        #min_temp_ss = np.min(temp_in_monthly_ss)
        #min_temp_fw = np.min(temp_in_monthly_fw)
        #min_temp = np.minimum(min_temp_ss,min_temp_fw)
        #estimate_T_param_calib = np.maximum(estimate_T_param_calib,-min_temp)
        #Test 2
        #T_param bounded to be > 8  and Func temp bounded to be <4.5
        estimate_T_param_calib = np.maximum(estimate_T_param_calib,15.)
        #test3
        estimate_T_param_calib = np.minimum(estimate_T_param_calib,30.)

        SITE_estimate_T_param_calib[j]=estimate_T_param_calib
        print 'estimated T_param_calib'
        print estimate_T_param_calib


        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
        #Run spinup (calibrated param)
        spinup_c_calib=spinup(SOC_lcs_2015,plant_in,Corg_in,clay_2009,temp_in_monthly_ss,water_in_monthly_ss,
                PET_in_monthly_ss,c_vegetated,t_param,estimate_T_param_calib,DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM)
        spinup_c_sum_calib=np.sum(spinup_c_calib) #tC/ha


        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<   
        #Run forward (calibrated param)
        fwd_c_calib=forward_ACTIVE(SOC_lcs_2015,n_an,spinup_c_calib,plant_in, Corg_in,clay_2009,temp_in_monthly_fw,water_in_monthly_fw,PET_in_monthly_fw,c_vegetated,t_param,estimate_T_param_calib,DPM_RPM_HUM_frac,DPM_RPM_HUM_frac_FYM)

        #Dynamics_all_calib=np.concatenate((np.array([spinup_c_calib]),fwd_c_calib))
        Dynamics_sum_calib = np.sum(fwd_c_calib,axis=1)
        #print 'Dynamics calib'
        #print Dynamics_sum_calib

        #2015 active (year 0)
        ACT_SOC_2015_calib[j] = fwd_c_calib[0,2] #BIO pool

        #SOC 2018
        SOC_2018_calib[j] = Dynamics_sum_calib[3]


        SITE_SOC_model_calib = Dynamics_sum_calib[0]
        SOC_init_model_calib[j]=SITE_SOC_model_calib

        if(save_mode):
            np.save(out_SOC_CAL,Dynamics_sum_calib)

        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
        #Run inverse modelling to 4p1000

        target = SITE_SOC_model_calib*0.004

        #litter income prior (1D) ->splitted in J_new
        in_opt = C_mean*(1+0.004)

        print ' '

        opt_mean_calib=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt_calib = opt_mean_calib.x[0]
        #print "SLSQP: Optimum solution:", litter_opt_calib

        #calculate percentage increase/decrease of inputs
        delta_input_calib = (litter_opt_calib-C_mean)*12. #tC/ha/yr
        input_change_calib=(delta_input_calib)*100./C_mean

        C_add_4p1000_calib[j]=delta_input_calib

        if(print_mode):
            print ' '
            print 'CALIBRATED SIMULATIONS'
            print 'Initial litter (tC/ha/month):'
            print C_mean
            print '4p1000 litter (tC/ha/month):'
            print litter_opt_calib

            print "% change of litter inputs:",input_change_calib

            print ' '

            #Check 4p1000 target is reached
            END=predict_c_calib.shape[0]-1
            #print 'END',END
            C_fin=np.sum(predict_c_calib[END])
            #C_init=np.sum(predict_c_calib[0])
            C_init=SITE_SOC_model_calib
            SUMO=C_fin-C_init
            Target_reached=(C_fin-C_init)/(C_init*n_an)
            if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
            else:
                print "Target not reached"
                print 'C init', C_init
                #print predict_c_calib[0]
                print SITE_SOC_model_calib
                print 'C_fin', C_fin
            print ' '


    if(print_mode):
        print ' '
        print '---'
        print 'Versions comparison'
        print 'Non calibrated Stationary solution (2015): '
        print SITE_SOC_model
        print 'Calibrated Stationary solution (2015): '
        print SITE_SOC_model_calib

        print 'Observed SOC in 2009:'
        print SOC_lcs_2009
        print 'Observed SOC in 2015:'
        print SOC_lcs_2015
        print ' '
        print 'Error - % non explained - (non calib):', (SITE_SOC_model-SOC_lcs_2015)/SOC_lcs_2015*100.
        print 'Error - % non explained - (calib):', (SITE_SOC_model_calib-SOC_lcs_2015)/SOC_lcs_2015*100.
        if(optimization_4p1000 and optimization_4p1000_CALIB):
            print ' '
            print 'Additional C input to the 4p1000:'
            print 'Non calib:',litter_opt
            print 'Calib:',litter_opt_calib
            
        print '---'


df_outputs = pd.DataFrame({'Site':site_names_all,'SOC_2015_NCAL':SOC_init_model,'SOC_2015_CAL':SOC_init_model_calib,'ACT_2015_NCAL':ACT_SOC_2015,'ACT_2015_CAL':ACT_SOC_2015_calib,'SOC_2018_NCAL':SOC_2018,'SOC_2018_CAL':SOC_2018_calib,'Add_C_4p1000_NCAL':C_add_4p1000,'Add_C_4p1000_CAL':C_add_4p1000_calib,'Calib T_param':SITE_estimate_T_param_calib})


if(save_mode):
    df_outputs.to_csv(OUTPUT_files+'ROTHC_outputs_bounds'+nppinput_from+'.csv')
    out_SOC_NCAL.close()
    out_SOC_CAL.close()

tend = time.time()
init_time = tinit - tstart
tot_time=tend-tstart
print '--'
print 'END'
print 'Initialization time:',init_time
print 'Total time', tot_time
print 'Number of sites simulated:',N_sites
exit()

