############################################################
#AMGv1 and AMGv2
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


#############
#Type of simulation
AMGv2=False

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
loc_hum = ROOTDIR+'_lcs_AR_GR30cm_v3.1.xls'
OUTPUT_files=ROOTDIR+'OUTPUTS/AMG/'

if(save_mode):
        out_SOC_NCAL = open(OUTPUT_files+"SOCdyn_NCAL_AMG_2015_2099"+nppinput_from+".txt","wb")
        out_SOC_CAL = open(OUTPUT_files+"SOCdyn_CAL_AMG_2015_2099"+nppinput_from+".txt","wb")


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
def control_temp_func(temp_mean):
#############################################
	#temp_mean = mean annual temperature

	aT = 25. 
	cT = 0.120 #1/Kelvin 
	TRef = 15. #Celcius degrees

	bT = (aT - 1)*np.exp(cT*TRef)

	if(temp_mean>0):
		temp_func = aT/(1+bT*np.exp(-cT*temp_mean))
	else:
		temp_func = 0.

 	return temp_func


		
#############################################
def control_moist_func(water_mean,pot_evapot):
#############################################
	#water_mean = cumulative annual water inputs (precipitation and irrigation water)
	#pot_evapot = potential evapotranspiration (mm)

	aH = 3.0/100
	bH = 5.247 #1/m

	moist_func = 1/(1+ aH*np.exp(-bH*(water_mean-pot_evapot)/1000))

	return moist_func


#############################################
def control_clay_func(clay):
#############################################
	#clay (gclay/kgsoil)
	global AMGv2
	if(AMGv2):
		aM =2.519/1000 #AMGv2
	else:
		aM = 2.720/1000 #g/kg2
	
	clay_func = np.exp(-aM*clay)

	return clay_func

#############################################
def control_carbonate_func(carbonate):
#############################################	
	#carbonate (CaCO3) content (gCaCo3/kgsoil)
	global AMGv2
	if(AMGv2):
		cM = 1.50/1000 #AMGv2
	else:
		cM = 1.67/1000  #g/kg
	
	carbonate_func = 1/(1+cM*carbonate)
	return carbonate_func

#In AMGv2
#############################################
def control_pH_func(pH):
#############################################
	apH = 0.112 
	bpH = 8.5
	pH_func = np.exp(-apH*(pH-bpH)**2)
	return pH_func

#In AMGv2
#############################################
def control_CN_ratio_func(CN_ratio):
#############################################
	aCN = 0.060 
	bCN = 11.
	CN_ratio_func = 0.8*np.exp(-aCN*(CN_ratio-bCN)**2) + 0.2
	return CN_ratio_func

#############################################
#Mineralization rate function
#############################################
#k0 is the potential mineralization rate set at 0.165 and 0.290 for AMGv1 and AMGv2
#To be optimized
#############################################
if(AMGv2):
	def k_rate_func(k0, temp_mean, water_mean, pot_evapot, clay, carbonate, pH, CN_ratio):
		k_rate = k0*control_temp_func(temp_mean)*control_moist_func(water_mean,pot_evapot)*control_clay_func(clay)*control_carbonate_func(carbonate)*control_pH_func(pH)*control_CN_ratio_func(CN_ratio)
		return k_rate
else:
	def k_rate_func(k0, temp_mean, water_mean, pot_evapot, clay, carbonate):
		k_rate = k0*control_temp_func(temp_mean)*control_moist_func(water_mean,pot_evapot)*control_clay_func(clay)*control_carbonate_func(carbonate)
		return k_rate

#############################################
def humcoef_func(CN_crop):
#############################################
	ecoef = 0.69
	fcoef = 11.2
	humcoef_AG = 1- (ecoef*CN_crop)/(fcoef+CN_crop)
	humcoef_BG = 0.39
	#add other humcoef for OM
	humcoef = np.array([humcoef_AG,humcoef_BG])
	return humcoef

#############################################
#If humcoefficients already defined
def forward_ACT(n_an,init,plantin,EOMin,temp_mean,water_mean,pot_evapot,k0,history_land,humcoef_plantin,humcoef_EOMin,clay,carbonate):
#############################################
	#litterin 2x array (AG,BG)
	#humcoef 2x array (AG,BG) humus coefficient
	global AMGv2

	#number of years = n_an because you don't have initialization
	ACT_cpools_tyear = np.zeros(n_an)
	#Set first year of active pool
	ACT_cpools_tyear[0]=init*(1-history_land)
	dt=1
	for x in range(1,n_an):
		temp_t = temp_mean[x]
		water_t = water_mean[x]
		pot_evapot_t = pot_evapot[x]
		ACT_current= ACT_cpools_tyear[x-1]
		if(AMGv2):
			kk_ma = k_rate_func(k0,temp_t, water_t, pot_evapot_t, clay, carbonate, pH, CN_ratio)
		else:
			kk_ma = k_rate_func(k0,temp_t, water_t, pot_evapot_t, clay, carbonate)

		ACT_next = ACT_current + (np.sum(plantin*humcoef_plantin+EOMin*humcoef_EOMin) - kk_ma*ACT_current)*dt
		ACT_cpools_tyear[x]=ACT_next

	return ACT_cpools_tyear	

#############################################
def SOC_stable_func(init,history_land):
#############################################

	SOC_stable = init*history_land
	
	return SOC_stable

#############################################
def SOC_tot_func(SOC_stable, SOC_act):
#############################################
	SOC_tot = SOC_stable + SOC_act
	return SOC_tot


############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):

    global NEW_ITER
    global n_an
    global clay_2009,CaCO3_2009, temp_in_yearly_fw,water_in_yearly_fw,PET_in_yearly_fw,arable_scoef
    global hum_Corg_av,humcoef_wa

    global C_mean,Corg_2009,Corg_in,SOC_lcs_2015
    global print_mode


    #Suppose Corg does not change
    if(Corg_2009<np.sum(in_new)): #if generated in_new is higher than Corg_2009 -> Corg is kept the same
        Corg_new = Corg_2009
        EOM_new = Corg_in
    else:   #Otherwise in_new is plant only 
        Corg_new = 0.
        EOM_new = np.array([0.,0.])

    #plant new is the difference between total (in_new) and Corg (EOM_new)
    plant_tot_new = in_new-EOM_new


    #4 per 1000 between 2015-2099
    if(opt_NONCAL): #Non calibrated
        global predict_c, k0_prior
        #Active
        predict_act=forward_ACT(n_an,SOC_lcs_2015,plant_tot_new,EOM_new,temp_in_yearly_fw,water_in_yearly_fw,PET_in_yearly_fw,k0_prior,arable_scoef,humcoef_wa,hum_Corg_av,clay_2009,CaCO3_2009)
        #Stable
        predict_stable = SOC_stable_func(SOC_lcs_2015,arable_scoef)
        #Total
        predict_c = predict_act+predict_stable

        #print 'predict_act'
        #print predict_act
        #print 'k0:', k0_prior
        

        J_new = abs(target*n_an - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

    elif(opt_CAL): #Calibrated
        global predict_c_calib, estimate_k0_param_calib
        #Active
        predict_act_calib=forward_ACT(n_an,SOC_lcs_2015,plant_tot_new,EOM_new,temp_in_yearly_fw,water_in_yearly_fw,PET_in_yearly_fw,estimate_k0_param_calib,arable_scoef,humcoef_wa,hum_Corg_av,clay_2009,CaCO3_2009)
        #print 'predict_act_calib'
        #print predict_act_calib
        #print 'k0 calib:', estimate_k0_param
        #Stable
        predict_stable_calib = SOC_stable_func(SOC_lcs_2015,arable_scoef)
        #Total
        predict_c_calib = predict_act_calib+predict_stable_calib

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
#############################################
#Fraction of initial SOC in the stable part
#############################################
arable_scoef=0.65
grassland_scoef=0.4


k0_prior=0.165 #need to be optimized

######################################################
#For each site: Set SITE name and experiment duration
#####################################################

C_input_exp = pd.read_csv(loc_exp)
df_humcoef = pd.read_excel(loc_hum, 'Sheet2',usecols="B:F",nrows=37)

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

#Reset index
#C_input_exp = C_input_exp.reset_index(drop=True)
#SITE_ID = C_input_exp['POINT_ID,N,11,0'].reset_index(drop=True)
SITE_LC15 = C_input_exp['LC15,C,4'].reset_index(drop=True)

#>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<
#Calculate weighted average humcoef and S:R ratios
#>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<

#List of crop species codes
species_indata = C_input_exp['LC15,C,4'].unique()[0:len(C_input_exp)]
species_indata = map(str, species_indata)

sum_SR_exp = 0
sum_AG_exp = 0
for i in species_indata:
        #Count number of sites for each crop species
        count_exp_sp  = np.float(C_input_exp[C_input_exp['LC15,C,4'] == i].shape[0])
        #Extract S:R for species i
        SR_ratio_sp = np.float(df_humcoef[df_humcoef['Code']== i]['SR_ratio'])
        #Extract AGhumcoef for species i
        hum_AG_sp = np.float(df_humcoef[df_humcoef['Code']== i]['humcoefAG'])
        #Nominatore weighted av
        sum_SR_exp+=SR_ratio_sp*count_exp_sp
        sum_AG_exp+=hum_AG_sp*count_exp_sp

#Weighted average S:R
weight_av_SR = sum_SR_exp/len(C_input_exp)

#Weighted average humcoefAG
hum_AG_plant = sum_AG_exp/len(C_input_exp)

#Hum coefficient BG
hum_BG_plant = 0.4

humcoef_wa = np.array([hum_AG_plant,hum_BG_plant])

#Average humcoefficient for FYM (from Levavasseur et al. 2020, optimized h)
hum_Corg_av = 0.548


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
SITE_SOC_lcs_2009 = C_input_exp['SOC_lcs,N,19,13'].reset_index(drop=True) #OC stock 2009-12 from OC* (1-VRG)* BD  * cm  (if 30 cm, we assume a stratification only for grassland 20 cm *1.27)
SITE_SOC_lcs_sd_2009 = C_input_exp['SOCsd_lcs,N,19,13'].reset_index(drop=True) #confidence interval by Montecarlo
SITE_SOC_tr_2009 = C_input_exp['SOCtr,N,19,15'].reset_index(drop=True) #ratio passive:total SOC pools at 2009
SITE_Corg_2009 = C_input_exp['Corg,N,19,12'].reset_index(drop=True) #possible C applied with manure (kg/ha) Derived from overlay with large-scale database (no direct point information)

#SITE_SOC_lcs_2015 = C_input_exp['SOC_lcs15,N,19,12'].reset_index(drop=True) #SOC stocks 2015
SITE_SOC_lcs_2015 = C_input_exp['SOCstocks_2015'].reset_index(drop=True) #SOC stocks 2015

#Other variables for statistical calibration
#SITE_pH_2009 = C_input_exp['pH_in_CaCl,N,19,15'].reset_index(drop=True)
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

#>>>>>>>> CLIMATE FORCING  <<<<<<<<<<<<<<
#print site_names_all
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

    #Calcola forcing for simulations (2015-2099)
    temp_in_yearly_fw = temp_in_yearly[n_an_ss:]
    PET_in_yearly_fw = PET_in_yearly[n_an_ss:]
    water_in_yearly_fw = water_in_yearly[n_an_ss:]

    #Append forcint for all sites
    SITE_temp_in_yearly_fw.append(temp_in_yearly_fw)
    SITE_PET_in_yearly_fw.append(PET_in_yearly_fw)
    SITE_water_in_yearly_fw.append(water_in_yearly_fw)


#>>>>>>>> END_OF_INITIALIZATION <<<<<<<<<<<<<<<<<<<<<<<<<<<<

tinit = time.time()

################################################
# OPTIMIZATION CONSTRAINTS
#ordine AG,BG
################################################
def constr1(x):
        global a_constr,ab_ratio
        return x[0]-a_constr*x[1]*ab_ratio

def constr2(x):
        global b_constr,ab_ratio
        return b_constr*x[1]*ab_ratio-x[0]

######################################
#Set bounds and constraints for the optimization
#########
bnds=[(0,15),(0,15)]

a_constr=1
b_constr=1
ab_ratio_var=(1-a_constr)*100
con1={'type':'ineq','fun':constr1}
con2={'type':'ineq','fun':constr2}
cons=[con1,con2]

#>>>>>>>> MAIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C_add_4p1000 = np.zeros(N_sites)
C_add_4p1000_calib = np.zeros(N_sites)
SOC_2018 = np.zeros(N_sites)
SOC_2018_calib = np.zeros(N_sites)
SITE_estimate_k0_param_calib = np.zeros(N_sites)

for j in range(0,N_sites):
    print ' '
    print 'START'
    print '---'
    print 'Site: ',site_names_all[j]

    #Upload forcing 
    temp_in_yearly_fw = SITE_temp_in_yearly_fw[j]
    PET_in_yearly_fw = SITE_PET_in_yearly_fw[j]
    water_in_yearly_fw = SITE_water_in_yearly_fw[j]

    #Upload soil data
    Cin_MOD = SITE_Cin_MOD.iloc[j] #MODIS tC/ha/yr
    npp_in = SITE_Cin_hanpp.iloc[j] #fraction
    Cin_MOD = Cin_MOD*npp_in #actual npp that goes into the soil
    Corg_2009 = SITE_Corg_2009.iloc[j]/1000. #kg/ha/yr to t/ha/yr

    #Code species
    code_crop_site_2015 = np.str(SITE_LC15.iloc[j])
    #Hum coefficient AG 2015
    hum_AG_plant_2015 = np.float(df_humcoef[df_humcoef['Code']==code_crop_site_2015]['humcoefAG'])

    #Convert total C in from crops to above and below pools from S:R ratios
    #below = Cin/(1+S:R)
    #above = Cin - below
    SR_plant_2015 = np.float(df_humcoef[df_humcoef['Code']==code_crop_site_2015]['SR_ratio'])
    BG_2015 = Cin_MOD/(1+SR_plant_2015)
    AG_2015 = Cin_MOD - BG_2015

    #print 'Cin_MOD',Cin_MOD
    #print 'BG_2015',BG_2015
    #print 'AG_2015',AG_2015
    #print 'BG_2015+AG_2015',BG_2015+AG_2015

    #above and below weighted averages over the database
    BG_wa = Cin_MOD/(1+weight_av_SR)
    AG_wa = Cin_MOD - BG_wa

    #Dividi Cin in different pools (weighted average S:R)
    plant_in = np.array([AG_wa,BG_wa])
    #Allocate 50% above and 50% below ?
    #alloc_AG_Corg = 0.5
    #Test2
    alloc_AG_Corg = 0.9
    alloc_BG_Corg = 1-alloc_AG_Corg
    Corg_in = np.array([Corg_2009*alloc_AG_Corg,Corg_2009*alloc_BG_Corg])

    Cin_array = plant_in+Corg_in # Array of C in total (plant + EOM)
    C_mean = np.sum(Cin_array) #C in total

    clay_2009 = SITE_clay_2009.iloc[j]*10.#% to g/kg
    CaCO3_2009 = SITE_CaCO3_2009.iloc[j]
    SOC_lcs_2009 = SITE_SOC_lcs_2009.iloc[j]
    SOC_lcs_sd_2009 = SITE_SOC_lcs_sd_2009.iloc[j]
    SOC_tr_2009 = SITE_SOC_tr_2009.iloc[j]
    SOC_lcs_2015 = SITE_SOC_lcs_2015.iloc[j]


    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #No spinup in AMG -> SOC init = SOC_lcs_2015
    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<   
    #Run forward
    SOC_stable = SOC_stable_func(SOC_lcs_2015,arable_scoef)
    SOC_act = forward_ACT(n_an,SOC_lcs_2015,plant_in,Corg_in,temp_in_yearly_fw,water_in_yearly_fw,PET_in_yearly_fw,k0_prior,arable_scoef,humcoef_wa,hum_Corg_av,clay_2009,CaCO3_2009)
    Dynamics_sum = SOC_tot_func(SOC_stable,SOC_act)

    print Dynamics_sum
    #print len(Dynamics_sum)

    #SOC 2018
    SOC_2018[j] = Dynamics_sum[3]

    if(save_mode):
        np.save(out_SOC_NCAL,Dynamics_sum)


    #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
    #Run inverse modelling to 4p1000


    if(optimization_4p1000):
        opt_NONCAL = True
        opt_CAL = False
        target = SOC_lcs_2015*0.004

        #litter income prior (2D)
        in_opt = Cin_array*(1+0.004)

        print ' '

        opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt = opt_mean.x
        #print "SLSQP: Optimum solution:", litter_opt

        total_opt_in=np.sum(litter_opt)

        #calculate percentage increase/decrease of inputs
        delta_input = (total_opt_in-C_mean) #tC/ha/yr
        print '***********'
        print 'total_opt_in',total_opt_in
        print 'delta_input',delta_input
        input_change=(total_opt_in-C_mean)*100./C_mean

        C_add_4p1000[j]=delta_input

        if(print_mode):
            #litter_opt_array=np.array([litter_opt*DPM_RPM_HUM_frac_mix[0],litter_opt*DPM_RPM_HUM_frac_mix[1],litter_opt*DPM_RPM_HUM_frac_mix[2]])
            print ' '
            print 'NON CALIBRATED SIMULATIONS'
            print 'k0:', k0_prior
            print 'Initial litter (tC/ha/yr):'
            print C_mean
            print '4p1000 litter (tC/ha/yr):'
            print total_opt_in

            print "% change of litter inputs:",input_change

            print ' '

            #Check 4p1000 target is reached
            END=predict_c.shape[0]-1
            #print 'END',END
            C_fin=np.sum(predict_c[END])
            #C_init=np.sum(predict_c[0])
            C_init=SOC_lcs_2015
            SUMO=C_fin-C_init
            Target_reached=(C_fin-C_init)/(C_init*n_an)
            if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
                print 'C init', C_init
                print 'C_fin', C_fin
            else:
                print "Target not reached"
                print 'C init', C_init
                #print predict_c[0]
                print SOC_lcs_2015
                print 'C_fin', C_fin
            print ' '

    #>>>>>>>>>>>>>>>> 4p1000 with optimized param <<<<<<<<<<<<<

    if(optimization_4p1000_CALIB):
        opt_NONCAL = False
        opt_CAL = True

        #########################################################
        #Derive calibrated parameters from statistical functions
        #########################################################
        CN_2015 = SITE_CN_2015.iloc[j]
        WET2 = SITE_WET2.iloc[j]
        #--
        loc_functions = '/home/surface7/ebruni/MULTIMODEL_EUROPA/CALIB_PARAM_FUNC/V8_wetness2/'
        func_df = pd.read_csv(loc_functions+'k0_calib_func.csv')
        interc_coef=func_df.loc['(Intercept)']
        SOC_init_coef=func_df.loc['Initial_SOC']
        temp_coef=func_df.loc['Temp']
        PET_coef=func_df.loc['PET']
        Cin_coef=func_df.loc['Litter_in'] #litter in in tC/ha/yr
        soil_CN_coef=func_df.loc['Soil.C.N']
        WET2_coef = func_df.loc['Wetness_class2']
        estimate_k0_param_calib = interc_coef+SOC_init_coef*SOC_lcs_2015+temp_coef*np.mean(temp_in_yearly_fw)+PET_coef*np.mean(PET_in_yearly_fw)+Cin_MOD*Cin_coef+soil_CN_coef*CN_2015+WET2_coef*WET2
        estimate_k0_param_calib = np.float(estimate_k0_param_calib)
        #0<k0<1
        estimate_k0_param_calib = np.minimum(estimate_k0_param_calib,1.)
        estimate_k0_param_calib = np.maximum(estimate_k0_param_calib,0.)

        SITE_estimate_k0_param_calib[j]=estimate_k0_param_calib


        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<   
        #Run forward (calibrated param)
        SOC_stable_calib = SOC_stable_func(SOC_lcs_2015,arable_scoef)
        SOC_act_calib = forward_ACT(n_an,SOC_lcs_2015,plant_in,Corg_in,temp_in_yearly_fw,water_in_yearly_fw,PET_in_yearly_fw,estimate_k0_param_calib,arable_scoef,humcoef_wa,hum_Corg_av,clay_2009,CaCO3_2009)
        Dynamics_sum_calib = SOC_tot_func(SOC_stable_calib,SOC_act_calib)


        print Dynamics_sum_calib
        #print len(Dynamics_sum)

        #SOC 2018
        SOC_2018_calib[j] = Dynamics_sum_calib[3]

        if(save_mode):
            np.save(out_SOC_CAL,Dynamics_sum_calib)

        #>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<
        #Run inverse modelling to 4p1000

        target = SOC_lcs_2015*0.004

        #litter income prior (2D)
        in_opt = Cin_array*(1+0.004)

        print ' '

        opt_mean_calib=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt_calib = opt_mean_calib.x
        #print "SLSQP: Optimum solution:", litter_opt

        total_opt_in_calib=np.sum(litter_opt_calib)

        #calculate percentage increase/decrease of inputs
        delta_input_calib = (total_opt_in_calib-C_mean) #tC/ha/yr
        print 'total_opt_in_calib',total_opt_in_calib
        print 'delta_input_calib',delta_input_calib
        print '***********'
        input_change_calib=(total_opt_in_calib-C_mean)*100./C_mean

        C_add_4p1000_calib[j]=delta_input_calib

        if(print_mode):
            print ' '
            print 'CALIBRATED SIMULATIONS'
            print 'k0 calib:'
            print estimate_k0_param_calib
            print 'Initial litter (tC/ha/yr):'
            print C_mean
            print '4p1000 litter (tC/ha/yr):'
            print total_opt_in_calib

            print "% change of litter inputs:",input_change_calib

            print ' '

            #Check 4p1000 target is reached
            END=predict_c_calib.shape[0]-1
            #print 'END',END
            C_fin=np.sum(predict_c_calib[END])
            #C_init=np.sum(predict_c[0])
            C_init=SOC_lcs_2015
            SUMO=C_fin-C_init
            Target_reached=(C_fin-C_init)/(C_init*n_an)
            if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
                print 'C init', C_init
                print 'C_fin', C_fin
            else:
                print "Target not reached"
                print 'C init', C_init
                #print predict_c[0]
                #print SOC_lcs_2015
                print 'C_fin', C_fin
            print ' '


df_outputs = pd.DataFrame({'Site':site_names_all,'SOC_2018_NCAL':SOC_2018,'SOC_2018_CAL':SOC_2018_calib,'Add_C_4p1000_NCAL':C_add_4p1000,'Add_C_4p1000_CAL':C_add_4p1000_calib,'Calib k0':SITE_estimate_k0_param_calib})


if(save_mode):
    df_outputs.to_csv(OUTPUT_files+'AMG_outputs'+nppinput_from+'.csv')
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


