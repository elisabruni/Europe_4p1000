############################################################
import sys
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import least_squares
from random import gauss
import xlrd
import pandas as pd
import time
import datetime
from datetime import datetime


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
OUTPUT_files=ROOTDIR+'OUTPUTS/ICBM/'

if(save_mode):
        out_SOC_NCAL = open(OUTPUT_files+"SOCdyn_NCAL_ICBM_2015_2099_k1_bounds"+nppinput_from+".txt","wb")
        out_SOC_CAL = open(OUTPUT_files+"SOCdyn_CAL_ICBM_2015_2099_k1_bounds"+nppinput_from+".txt","wb")


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

#################################################
#Convert monthly data to daily for NON cumulative var
#################################################
def MonthlyToDaily(filein):
        c=0
        x=0
        var_days=np.zeros(np.int(len(filein)/12.*365))
        month_days = [31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
        rep_month_days = np.tile(month_days,len(filein)/12)
        daily_var_monthlyrep = filein
        lent = np.int(len(rep_month_days))
        for j in range(0,len(daily_var_monthlyrep)):
            days_in_month = np.int(rep_month_days[j])
            var_days[c:c+days_in_month]=np.repeat(daily_var_monthlyrep[x],days_in_month)
            c+=days_in_month
            x+=1
        return var_days

################################
################################
#ENVIRONMENTAL FUNCTIONS
################################
###################
#Moisture function
###################
def rwat2(clay, water_in_m3):

        # mcw = wilting point
        # mcfc = field capacity

        #From katterer et al 2006
        mcfc=0.2699+0.3247*clay
        mcw=0.0284+0.3790*clay

        #water_in_m3 = np.maximum(water_in_m3,alpha*mcw)

        if(water_in_m3<alpha*mcw):
                re_wat=0.18*(water_in_m3/(alpha*mcw))

        elif((water_in_m3>alpha*mcw or water_in_m3==alpha*mcw) and (water_in_m3<gamma*mcfc or water_in_m3==gamma*mcfc)):
                re_wat=0.18+(1.-0.18)*((water_in_m3-alpha*mcw)/(gamma*mcfc-alpha*mcw))

        elif(water_in_m3>gamma*mcfc):
                re_wat=1+(1-rmin)*((water_in_m3-gamma*mcfc)/(gamma*mcfc-mcfc))

        return re_wat

###################
#Temp function
###################
def func_temp(Tsoil,Tmin):
        temp_func=np.maximum(0.,((Tsoil-Tmin)**2)/((30-Tmin)**2))
        return temp_func

##################
#Equations
##################

def Young_SS(litterin,k1,r):
        #litterin_0 = litterin[0]
        litterin_0 = litterin
        Y_0 = litterin_0/(k1*r)
        return Y_0

def Young_C(n_an,litterin,k1,r,Y_SS):
        t = np.arange(1,n_an+1)
        Y_t = litterin/(k1*r) + (Y_SS - litterin/(k1*r))*np.exp(-k1*r*t)
        return Y_t

def Old_SS(plantin,EOMin,k2,r,hplant,hEOM):
        litterin_0 = np.array([plantin,EOMin])  #C in from plant and EOM
        h_0 = np.array([hplant,hEOM])   #hum coeff plant and EOM
        hum_C_0 = np.sum(litterin_0*h_0) #humified C
        O_0 = hum_C_0/(k2*r)
        return O_0

def Old_C(n_an,plantin,EOMin,k1,k2,r,hplant,hEOM,Y_SS,O_SS):
        t = np.arange(1,n_an+1)
        litterin_t = np.array([plantin,EOMin])  #C in from plant and EOM
        h_t = np.array([hplant,hEOM])   #hum coeff plant and EOM
        h_t_av = np.mean(h_t[np.nonzero(h_t)]) #mean among hum coeff (exclude if h=0)
        hum_C_t = np.sum(litterin_t*h_t) #humified C
        #O_t = h*litterin/(k2*r) + (O_SS - h*litterin/(k2*r)-h*(k1*r*Y_SS-litterin)/(r*(k2-k1)))*np.exp(-k2*r*t)+(h*(k1*r*Y_SS-litterin)/(r*(k2-k1)))*np.exp(-k1*r*t)
        O_t = hum_C_t/(k2*r) + (O_SS - hum_C_t/(k2*r)- (h_t_av*k1*r*Y_SS - hum_C_t)/(r*(k2-k1)))*np.exp(-k2*r*t)+((h_t_av*k1*r*Y_SS-hum_C_t)/(r*(k2-k1)))*np.exp(-k1*r*t)
        return O_t




############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):

    global NEW_ITER
    global n_an,k1, k2

    global C_mean,Corg_2009
    global print_mode


    #Suppose Corg does not change
    if(Corg_2009<np.sum(in_new)): #if generated in_new is higher than Corg_2009 -> Corg is kept the same
        Corg_new = Corg_2009
    else:   #Otherwise in_new is plant only 
        Corg_new = 0.

    #plant new is the difference between in_new and Corg
    plant_tot_new = in_new-Corg_new


    #4 per 1000 between 2015-2099
    if(opt_NONCAL): #Non calibrated
        global young_steady, old_steady, predict_c
        global re_func_rate_fw
        #Young
        predict_young_c = Young_C(n_an,in_new,k1,re_func_rate_fw,young_steady)
        #Old
        predict_old_c = Old_C(n_an,plant_tot_new,Corg_new,k1,k2,re_func_rate_fw,hstraw,h_fym,young_steady,old_steady)
        #Total
        predict_c = predict_young_c+predict_old_c

        J_new = abs(target*n_an - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

    elif(opt_CAL): #Calibrated
        global young_steady_calib, old_steady_calib, predict_c_calib
        global estimate_k1_param_calib

        #Young
        predict_young_c_calib = Young_C(n_an,in_new,estimate_k1_param_calib,re_func_rate_fw,young_steady_calib)
        #Old
        predict_old_c_calib = Old_C(n_an,plant_tot_new,Corg_new,estimate_k1_param_calib,k2,re_func_rate_fw,hstraw,h_fym,young_steady_calib,old_steady_calib)
        #Total
        predict_c_calib = predict_young_c_calib+predict_old_c_calib

        J_new = abs(target*n_an - np.sum(predict_c_calib[predict_c_calib.shape[0]-1]-predict_c_calib[0]))

    #print 'OBJ FUN', J_new
    #print 'predict_c',predict_c

    NEW_ITER+=1
    if(print_mode):
        if ( NEW_ITER % 100 == 0 ):
            print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new

    param_opt=in_new
    return J_new

############################################################
#Initialization
NEW_ITER = 0
CHI2_PRINT_FREQUENCY=2

#############################################
#PARAMETERS

#--------
#Moisture
#--------
gamma=0.85 #From graph in Karlsson et al 2011
rmin=0.55 #From graph in Karlsson et al 2011
alpha=0.5 #From Fortin et al. 2011
#----------
#Temperature
#----------
Tmin=-3.8 #From Karlsson et al 2011
#----------

#Decomposition rates
k1 = 0.8
k2 = 0.00605
########################
#humification coefficient
#and enviromental parameter
#Steady state
h0 = 0.125
#r0 = 1.
#Straw (without N)
hstraw = 0.125
h_fym=0.31
h_ss=0.47


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

#Appends values of variable for each site
SITE_Cin_MOD = C_input_exp['NPPmod,N,19,14'].reset_index(drop=True) #average 2000-2009 NPP from MODIS 
SITE_Cin_hanpp = C_input_exp['MEAN_HANPP'].reset_index(drop=True) #fraction of npp that goes into the soil

SITE_clay_2009 = C_input_exp['clay,N,19,14'].reset_index(drop=True) #%
SITE_CaCO3_2009 = C_input_exp['CaCO3,N,19,13'].reset_index(drop=True) #CacO3 g/kg
SITE_pH_2009 = C_input_exp['pH_in_CaCl,N,19,15'].reset_index(drop=True)
SITE_SOC_lcs_2009 = C_input_exp['SOC_lcs,N,19,13'].reset_index(drop=True) #OC stock 2009-12 from OC* (1-VRG)* BD  * cm  (if 30 cm, we assume a stratification only for grassland 20 cm *1.27)
SITE_SOC_lcs_sd_2009 = C_input_exp['SOCsd_lcs,N,19,13'].reset_index(drop=True) #confidence interval by Montecarlo
SITE_SOC_tr_2009 = C_input_exp['SOCtr,N,19,15'].reset_index(drop=True) #ratio passive:total SOC pools at 2009
SITE_Corg_2009 = C_input_exp['Corg,N,19,12'].reset_index(drop=True) #possible C applied with manure (kg/ha) Derived from overlay with large-scale database (no direct point information)

#SITE_SOC_lcs_2015 = C_input_exp['SOC_lcs15,N,19,12'].reset_index(drop=True) #SOC stocks 2015
SITE_SOC_lcs_2015 = C_input_exp['SOCstocks_2015'].reset_index(drop=True) #SOC stocks 2015

#Other variables for statistical calibration
SITE_OC_2015 = C_input_exp['OC_15,N,19,13'].reset_index(drop=True) #g/kg
SITE_N_2015  = C_input_exp['N_15,N,19,14'].reset_index(drop=True)  #g/kg
SITE_CN_2015 = SITE_OC_2015/SITE_N_2015

SITE_WET2 = C_input_exp['wetness'].reset_index(drop=True)

#>>>>>>>>>>>>>#####<<<<<<<<<<<<<<

site_names_all = C_input_exp['POINT_ID,N,11,0'].unique()[0:len(C_input_exp)]
site_names_all = map(str, site_names_all)

N_sites=len(site_names_all)

#Simulations lenght
n_an_ss = 2014-2006+1 #yr ss
n_an = 2099-2015+1 #yr fwd
n_day_ss = n_an_ss*365
n_day = n_an*365
#>>>>>>>> CLIMATE FORCING  <<<<<<<<<<<<<<
SITE_re_func_rate_ss = []
SITE_re_func_rate_fw = []

#For calibration
SITE_temp_in_yearly_fw = []
SITE_PET_in_yearly_fw = []
SITE_water_in_yearly_fw = []
j=0
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
    #K to C (daily temp)
    temp_in_ar = np.array(temp_in_ar)-273.15

    ################
    #For calibration
    temp_in_yearly=DailyToYearly(temp_in_ar,1)

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

    #kg/m2/sec (daily interp) to mm/year (1kg/m2 = 1mm)
    PET_in_ar = np.array(PET_in_ar)*3600.*24.
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

    #Rain + snow (kg/m2/sec (daily) -> mm/year)
    water_tot = (PR_in_ar+PRSN_in_ar)*3600.*24.
    water_in_yearly=DailyToYearly(water_tot,0)

    #----

    ###############
    #Soil moisture
    ###############
    soil_moist =ROOTDIR+'ISIMIP_data/SOILMOIST_RCP26/'+site+'_SOILMOIST_2006_2099.txt'
    moist_in_ar_month=[]
    prev_line=0
    with open(soil_moist) as fileID:
        for _ in range(8): #skip first 8 lines
            next(fileID)
        for line in fileID:
            #Get the info and put it separately in an array
            sub1=line.split(":",1)[1].strip() #Eliminate white spaces

            #MISSING POINTS
            if(sub1=='....'):
                sub1=prev_line
            #---------------

            moist_in = float(sub1)
            prev_line = moist_in

            moist_in_ar_month.append(moist_in)

    moist_in_ar_month = np.array(moist_in_ar_month) #m3 m3 in top 18cm
    hum_in_ar = MonthlyToDaily(moist_in_ar_month) #daily moisture

    #Read clay value
    clay_2009 = SITE_clay_2009.iloc[j]/100. #% to g/g

    #Calcola forcing for spin-up (2006-2014)
    temp_in_daily_ss = temp_in_ar[0:n_day_ss]
    hum_in_daily_ss = hum_in_ar[0:n_day_ss]

    #Calcola forcing for simulations (2015-2099)
    temp_in_daily_fw = temp_in_ar[n_day_ss:]
    hum_in_daily_fw = hum_in_ar[n_day_ss:]


    #Calculate moist and temp functions (with daily values)
    re_day_ss = np.zeros(len(hum_in_daily_ss))
    for i in range(0,len(hum_in_daily_ss)):
            #Steady state
            moist_day_ss=hum_in_daily_ss[i]
            temp_day_ss=temp_in_daily_ss[i]

            re_wat_ss = rwat2(clay_2009, moist_day_ss)
            re_temp_ss = func_temp(temp_day_ss,Tmin)

            re_dayS = re_wat_ss*re_temp_ss
            re_day_ss[i]=re_dayS


    re_day_fw = np.zeros(len(hum_in_daily_fw))
    for i in range(0,len(hum_in_daily_fw)):
            #Forward
            moist_day_fw=hum_in_daily_fw[i]
            temp_day_fw=temp_in_daily_fw[i]

            re_wat_fw = rwat2(clay_2009, moist_day_fw)
            re_temp_fw = func_temp(temp_day_fw,Tmin)

            re_day = re_wat_fw*re_temp_fw
            re_day_fw[i]=re_day

    #Average re_day
    av_re_day_ss = np.mean(re_day_ss)
    av_re_day_fw = np.mean(re_day_fw)

    #Append forcint for all sites
    SITE_re_func_rate_ss.append(av_re_day_ss)
    SITE_re_func_rate_fw.append(av_re_day_fw)

    #For calibration
    temp_in_yearly_fw = temp_in_yearly[n_an_ss:]
    PET_in_yearly_fw = PET_in_yearly[n_an_ss:]
    water_in_yearly_fw = water_in_yearly[n_an_ss:]

    #For calibration
    SITE_temp_in_yearly_fw.append(temp_in_yearly_fw)
    SITE_PET_in_yearly_fw.append(PET_in_yearly_fw)
    SITE_water_in_yearly_fw.append(water_in_yearly_fw)

    j+=1

#Calculate average SITE_temp_func_rate and SITE_hum_func_rate
#Normalize against site in Sweden
#SITE id : 47304102 (coordinates: 59.82 N 17.28)
SITE_ID = C_input_exp['POINT_ID,N,11,0'].reset_index(drop=True)
swedish_site_index = SITE_ID[SITE_ID==47304102].index[0]
swedish_site_re_ss = SITE_re_func_rate_ss[swedish_site_index]
swedish_site_re_fw = SITE_re_func_rate_fw[swedish_site_index]

SITE_re_func_rate_norm_ss = SITE_re_func_rate_ss/swedish_site_re_ss
SITE_re_func_rate_norm_fw = SITE_re_func_rate_fw/swedish_site_re_fw


#>>>>>>>> END_OF_INITIALIZATION <<<<<<<<<<<<<<<<<<<<<<<<<<<<

tinit = time.time()

################################################
#Set bounds and constraints for the optimization
#########
bnds=[(0,30)]


#>>>>>>>> MAIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<
SOC_init_model = np.zeros(N_sites)
SOC_init_model_calib = np.zeros(N_sites)
SOC_2018 = np.zeros(N_sites)
SOC_2018_calib = np.zeros(N_sites)
ACT_SOC_2015 = np.zeros(N_sites)
ACT_SOC_2015_calib = np.zeros(N_sites)
C_add_4p1000 = np.zeros(N_sites)
C_add_4p1000_calib = np.zeros(N_sites)

SITE_estimate_k1_param_calib = np.zeros(N_sites)

for j in range(0,N_sites):
    print ' '
    print 'START'
    print '---'
    print 'Site: ',site_names_all[j]

    #Upload forcing 
    re_func_rate_ss = SITE_re_func_rate_norm_ss[j] #Normalized temp and moist func rate
    re_func_rate_fw = SITE_re_func_rate_norm_fw[j] #Normalized temp and moist func rate

    #Upload soil data
    Cin_MOD = SITE_Cin_MOD.iloc[j]/10. #MODIS tC/ha/yr to kg/m2/yr
    npp_in = SITE_Cin_hanpp.iloc[j] #fraction
    Cin_MOD = Cin_MOD*npp_in #actual npp that goes into the soil
    
    Corg_2009 = SITE_Corg_2009.iloc[j]/10000. #kg/ha/yr to kg/m2/yr

    C_mean = Cin_MOD+Corg_2009 #C in total

    SOC_lcs_2009 = SITE_SOC_lcs_2009.iloc[j]
    SOC_lcs_sd_2009 = SITE_SOC_lcs_sd_2009.iloc[j]
    SOC_tr_2009 = SITE_SOC_tr_2009.iloc[j]
    SOC_lcs_2015 = SITE_SOC_lcs_2015.iloc[j]


    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #Run spinup
    young_steady = Young_SS(C_mean,k1,re_func_rate_ss)
    old_steady = Old_SS(Cin_MOD,Corg_2009,k2,re_func_rate_ss,hstraw,h_fym)
    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #Run forward
    young_fw = Young_C(n_an,C_mean,k1,re_func_rate_fw,young_steady)
    old_fw = Old_C(n_an,Cin_MOD,Corg_2009,k1,k2,re_func_rate_fw,hstraw,h_fym,young_steady,old_steady)

    Dynamics_sum = young_fw+old_fw

    #print Dynamics_sum
    #print len(Dynamics_sum)

    #2015
    SITE_SOC_model = Dynamics_sum[0]
    SOC_init_model[j] = SITE_SOC_model

    #2015 active (year 0)
    ACT_SOC_2015[j] = young_fw[0]

    #SOC 2018
    SOC_2018[j] = Dynamics_sum[3]

    if(save_mode):
        np.save(out_SOC_NCAL,Dynamics_sum*10.) #tC/ha


    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #Run inverse modelling to 4p1000


    if(optimization_4p1000):
        opt_NONCAL=True
        opt_CAL=False
        target = SITE_SOC_model*0.004

        #litter income prior (2D)
        in_opt = C_mean*(1+0.004)

        print ' '

        opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt = opt_mean.x[0]
        #print "SLSQP: Optimum solution:", litter_opt

        #calculate percentage increase/decrease of inputs
        delta_input = (litter_opt-C_mean)*10. #tC/ha/yr
        input_change=(litter_opt-C_mean)*100./C_mean

        C_add_4p1000[j]=delta_input

        if(print_mode):
            #litter_opt_array=np.array([litter_opt*DPM_RPM_HUM_frac_mix[0],litter_opt*DPM_RPM_HUM_frac_mix[1],litter_opt*DPM_RPM_HUM_frac_mix[2]])
            print ' '
            print 'NON CALIBRATED SIMULATIONS'
            print 'Initial litter (tC/ha/yr):'
            print C_mean*10.
            print '4p1000 litter (tC/ha/yr):'
            print litter_opt*10.

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
        temp_in_yearly_fw = SITE_temp_in_yearly_fw[j]
        PET_in_yearly_fw = SITE_PET_in_yearly_fw[j]
        water_in_yearly_fw = SITE_water_in_yearly_fw[j]
        pH_2009 = SITE_pH_2009.iloc[j]
        CaCO3_2009 = SITE_CaCO3_2009.iloc[j]
        CN_2015 = SITE_CN_2015.iloc[j]
        WET2 = SITE_WET2.iloc[j]
        clay_2009 = SITE_clay_2009.iloc[j]/100. #% to g/g

        loc_functions = '/home/surface7/ebruni/MULTIMODEL_EUROPA/CALIB_PARAM_FUNC/V8_wetness2/'

        ###
        #k1
        ###
        func_df_k1 = pd.read_csv(loc_functions+'k1_calib_func.csv')
        interc_coef_k1=func_df_k1.loc['(Intercept)']
        SOC_init_coef_k1=func_df_k1.loc['Initial_SOC']
        temp_coef_k1=func_df_k1.loc['Temp']
        PET_coef_k1=func_df_k1.loc['PET']
        clay_coef_k1= func_df_k1.loc['Clay']
        CaCO3_coef_k1=func_df_k1.loc['Carbonate']
        pH_coef_k1 = func_df_k1.loc['pH']
        soil_CN_coef_k1 = func_df_k1.loc['Soil.C.N']
        WET2_coef_k1 = func_df_k1.loc['Wetness_class2']
        estimate_k1_param_calib = interc_coef_k1+SOC_init_coef_k1*SOC_lcs_2015+temp_coef_k1*np.mean(temp_in_yearly_fw)+PET_coef_k1*np.mean(PET_in_yearly_fw)+clay_coef_k1*clay_2009+CaCO3_2009*CaCO3_coef_k1+pH_2009*pH_coef_k1+CN_2015*soil_CN_coef_k1+WET2_coef_k1*WET2
        estimate_k1_param_calib = np.float(estimate_k1_param_calib)


        #Bounds
        #0<k1<20
        #Bounds1
        #estimate_k1_param_calib = np.minimum(estimate_k1_param_calib,50.)
        #Bounds2
        #estimate_k1_param_calib = np.minimum(estimate_k1_param_calib,30.)
        #Bounds3
        estimate_k1_param_calib = np.minimum(estimate_k1_param_calib,15.)
        estimate_k1_param_calib = np.maximum(estimate_k1_param_calib,0.1)


        print 'k1',estimate_k1_param_calib 

        SITE_estimate_k1_param_calib[j]=estimate_k1_param_calib


        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
        #Run spinup (calibrated param)
        young_steady_calib = Young_SS(C_mean,estimate_k1_param_calib,re_func_rate_ss)
        old_steady_calib = Old_SS(Cin_MOD,Corg_2009,k2,re_func_rate_ss,hstraw,h_fym)
        print 'young_steady_calib',young_steady_calib
        print 'old_steady_calib',old_steady_calib
        print 'n_an',n_an
        print 'C_mean',C_mean
        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
        #Run forward (calibrated param)
        young_fw_calib = Young_C(n_an,C_mean,estimate_k1_param_calib,re_func_rate_fw,young_steady_calib)
        old_fw_calib = Old_C(n_an,Cin_MOD,Corg_2009,estimate_k1_param_calib,k2,re_func_rate_fw,hstraw,h_fym,young_steady_calib,old_steady_calib)

        Dynamics_sum_calib = young_fw_calib+old_fw_calib


        print Dynamics_sum_calib
        #print len(Dynamics_sum)

        #2015
        SITE_SOC_model_calib = Dynamics_sum_calib[0]
        SOC_init_model_calib[j] = SITE_SOC_model_calib


        #2015 active (year 0)
        ACT_SOC_2015_calib[j] = young_fw_calib[0]

        #SOC 2018
        SOC_2018_calib[j] = Dynamics_sum_calib[3]

        if(save_mode):
            np.save(out_SOC_CAL,Dynamics_sum_calib*10.) #tC/ha

        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
        #Run inverse modelling to 4p1000
        
        target = SITE_SOC_model_calib*0.004

        #litter income prior (1D)
        in_opt = C_mean*(1+0.004)

        print ' '

        opt_mean_calib=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt_calib = opt_mean_calib.x[0]
        #print "SLSQP: Optimum solution:", litter_opt

        #calculate percentage increase/decrease of inputs
        delta_input_calib = (litter_opt_calib-C_mean)*10. #tC/ha/yr
        input_change_calib=(litter_opt_calib-C_mean)*100./C_mean

        C_add_4p1000_calib[j]=delta_input_calib

        if(print_mode):
            print ' '
            print 'NON CALIBRATED SIMULATIONS'
            print 'Initial litter (tC/ha/yr):'
            print C_mean*10.
            print '4p1000 litter (tC/ha/yr):'
            print litter_opt_calib*10.

            print "% change of litter inputs:",input_change_calib

            print ' '

            #Check 4p1000 target is reached
            END=predict_c_calib.shape[0]-1
            #print 'END',END
            C_fin=np.sum(predict_c_calib[END])
            #C_init=np.sum(predict_c[0])
            C_init=SITE_SOC_model_calib
            SUMO=C_fin-C_init
            Target_reached=(C_fin-C_init)/(C_init*n_an)
            if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
            else:
                print "Target not reached"
                print 'C init', C_init
                #print predict_c[0]
                print SOC_lcs_2015
                print 'C_fin', C_fin
            print ' '

df_outputs = pd.DataFrame({'Site':site_names_all,'SOC_2015_NCAL':SOC_init_model*10.,'SOC_2015_CAL':SOC_init_model_calib*10.,'ACT_2015_NCAL':ACT_SOC_2015*10.,'ACT_2015_CAL':ACT_SOC_2015_calib*10.,'SOC_2018_NCAL':SOC_2018*10.,'SOC_2018_CAL':SOC_2018_calib*10.,'Add_C_4p1000_NCAL':C_add_4p1000,'Add_C_4p1000_CAL':C_add_4p1000_calib,'Calib k1':SITE_estimate_k1_param_calib})


if(save_mode):
    df_outputs.to_csv(OUTPUT_files+'ICBM_outputs_k1_bounds'+nppinput_from+'.csv')
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


