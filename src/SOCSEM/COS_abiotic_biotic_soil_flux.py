# -*- coding: utf-8 -*-
################################################################################
#####            Whelan et al 2021 Soil OCS Exchange Model           ###########
#####            SOCSEM version 8.0.1                                ###########
################################################################################

import numpy as np

################################################################################
#####            version 8.0 notes                                   ###########
################################################################################
# SOCSSEM v8.0 was updated to include a greater emphasis on new field data.
# For example, though the shape of the OCS flux response curves with temperature
# are typically based on soil incubations, the overall relationships are often 
# scaled to field observations. Further notes are provided by the equations 
# of interest. 

#Changes from v8.0.0
#Stream-lined representation of error
#Cleaned up formatting
#Used flux_at_all_theta function instead of flux_theta_g_constant
#to prevent errors with declaring "CONSTANT_PARAM" in unexpected places
#For grasslands, Grass at other theta is at 26.9 VWC instead of 25

################################################################################
#####                 Soil COS Emission to Atmosphere                ###########
################################################################################
# By convention, emission to the atmosphere is positive,
# uptake from the atmosphere is negative. 

# The exponential shape of OCS production with temperature
# was first noted by Liu et al. 2010 in lab incubations
# and observed in the field by Maseyk et al. 2014
# I have subsequently observed this relationship
# with soils from many different biomes, e.g. Whelan et al., 2016
# This should vary with microbial biomass and N content
# however, those data are hard to come by at a workable resolution
# Here we use biome as a proxy

# These fits are from lab-based soil incubation experiments
# Formula cos_flux = a / (1 + b e-kx ), k > 0
# y = a to the right, set to the max flux observed in the field
# This maximum should not be achieved under environmental conditions
# y = 0 to the left, the y axis intercept passes through (0, a/(1+b))

# Max Fluxes are assigned as the largest flux observed 
# in the field (when available) or from incubation experiements
# in pmol OCS / m2 / sec
# Temperate Forest: 20. (Commane et al., 2015)
# Tropical Forest: Meredith et al 2018 incubations found highest Q10 numbers 
# for tropical sites (~3)
# and Whelan et al., 2016 incubated soil samples from the Peruvian Amazon
# Unfortunately, no additional rainforest measurements were performed
# so we arbitrarily triple the maximum incubation observation as 2.7*3
# Grassland 3.9 (Whelan and Rhew 2016)
# Agriculture 83. max scaled incubation flux observed with US-Bo1 flux site soil
# from Whelan et al.,  2016
# Commented out number after k and b are 1 standard deviation error

#r2 = 0.96 for PRU, Los Amigos Research Station, Peru [Whelan et al. 2016]
def OCS_rainforest_production(soil_temp):
    pru_a = 2.7*3 
    pru_k = 0.123581 # 1.46458979e-02
    pru_b = 205. # 1.01737470e+02
    return pru_a/(1+pru_b*np.exp(-pru_k*soil_temp)) #in pmol OCS m^-2 sec^-1

#We are assuming (likely incorrectly) that temperate and boreal curves are the same
#r2 = 0.997 for WRC, Willow Creek Forest Tall Tower Site [Whelan et al. 2016]
def OCS_forest_production(soil_temp):
    wrc_a = 20.
    wrc_k = 0.160745 # 5.96007774e-03
    wrc_b = 644.7 # 1.32871427e+02
    return wrc_a/(1+wrc_b*np.exp(-wrc_k*soil_temp)) #in pmol OCS m^-2 sec^-1

#r2 = 0.79 for DOE ARM SGP, a wheat field [Maseyk et al. 2014]
def OCS_ag_production(soil_temp):
    bond_a = 83.
    ok_k = .087689 # 5.60537742e-03
    ok_b = 146.9 # 3.18441017e+01
    return bond_a/(1+ok_b*np.exp(-ok_k*soil_temp)) #in pmol OCS m^-2 sec^-1
    
#r2 = 0.98 for STUNT, Stunt Ranch UC Reserve, CA, savannah [Whelan et al. 2016]
def OCS_grass_production(soil_temp):
    stunt_a = 3.9
    stunt_k = 0.115481 #1.50019288e-02
    stunt_b = 286. # 1.58117367e+02
    return  stunt_a/(1+stunt_b*np.exp(-stunt_k*soil_temp)) #in pmol OCS m^-2 sec^-1

#r2 = 0.64, Mollie Beattie Habitat Community, Port Aransas, TX, [Whelan et al., 2013]
def OCS_wetland_production(soil_temp):
    #y = a at high temperature
    tx_a = 295. #max value recorded by DeLaune et al. (2002)
    tx_k = 0.07855407
    tx_b = 41.16705972
    return tx_a/(1+tx_b*np.exp(-tx_k*soil_temp)) #in pmol OCS m^-2 sec^-1
    
################################################################################
#####                  Soil COS Uptake from Atmosphere               ###########
################################################################################
# This is an NO production model borrowed from Behrendt et al. 2014
# Repurposed for COS consumption
# Van Diest and Kesselmeier [2008] re-purposed a different NO production model, 
# but the spirit of the pursuit is the same

# To characterize OCS uptake, we first subtract OCS production from 
# lab soil incubations and field observations.
# Then estimate
# (1) maximum COS uptake, opt_flux
# (2) soil volumetric water content (theta_opt in VWC, between 0 and 100) at max uptake
# (3) the "other" OCS uptake flux_at_theta_g at a separate soil moisture, theta_g, where theta_g > theta_opt
# (4) how temperature alters the opt_flux and flux_at_theta_g

# This model won't support freezing temperatures
# We approximate fluxes at 0 C as 0 COS pmol m^-2 sec^-1, 
# though we know this is incorrect for some systems
# The role of lichen and bryophytes is woefully ignored

# For each temperature (bin or individual value), this calculates the soil moisture response curve
# The first equation is a curve "shape" value, a 
# Note that a is typically 1 to 3.  
# a increases to ~30 when optimum soil moisture approaches the "other" moisture, theta_g

def curve_shape_a(opt_flux, flux_at_theta_g, theta_g, theta_opt):
    Rj = opt_flux/flux_at_theta_g
    outs = np.log(Rj)
    inversed = (np.log(theta_opt/theta_g)) + (theta_g/theta_opt - 1)
    ins = 1/inversed
    return outs*ins 

#This curve is used for optimizing when theta_g is allowed to change
def flux_at_all_theta(all_theta, opt_flux, flux_at_theta_g, theta_g, theta_opt):
    a = curve_shape_a(opt_flux, flux_at_theta_g, theta_g, theta_opt)
    power_fct_increasing = (all_theta/theta_opt) ** a
    exp_fct_decreasing = np.exp(-1 * a * ((all_theta/theta_opt)- 1)) 
    return opt_flux * power_fct_increasing * exp_fct_decreasing

#This curve is used for optimizing when theta_g is set to CONSTANT PARAM
#When you have observations you'd like to fit for a specific theta_g
CONSTANT_PARAM = 35 #Arbitrary value to start with
def flux_theta_g_constant(all_theta, opt_flux, flux_at_theta_g,theta_opt):
    #theta g is set to CONSTANT_PARAM
    a = curve_shape_a(opt_flux, flux_at_theta_g, CONSTANT_PARAM, theta_opt)
    power_fct_increasing = (all_theta/theta_opt) ** a
    exp_fct_decreasing = np.exp(-1 * a * ((all_theta/theta_opt)- 1)) 
    return opt_flux * power_fct_increasing * exp_fct_decreasing
    
################################################################################
######            empirically derived coefficients                     #########
######           for uptake over varying temperatures                  #########
################################################################################
# The shape and behavior of OCS uptake was first characterized by
# Kesselmeier et al., 1999 with further data in van Diest and Kesselmeier (2008). 
# OCS uptake has an optimum temperature and soil moisture, with less uptake
# as conditions deviate from the optimum.

# Here, we subtract out abiotic production based on the equations above
# and search for the maximum uptake flux observed. 
# We declare that this is as close to the "optimum" uptake that we can estimate

# It appears that for boreal soils, the optimum is at a much higher temperature
# (25 to 30 C) versus other soils (around 15 C)
# Van Diest and Kesselmeier (2008) suggest that perhaps this was because of
# the presence of extromophiles in soil with large temperature swings

# First, we figure out the highest expected flux at various **temperatures**
# Second, we will feed these values into the relationship with **soil moisture**

# Commented numbers by fits are 1 standard deviation errors between model and observations

#GRASSLAND SOIL
#Based on Stunt Ranch field and lab data
#Sun et al. (2016) and Whelan et al. (2016)
def grass_opt_sw():
    return 12.5 #based on fits to field and lab Stunt Ranch Data, error ~1.9
def grass_opt_uptake(T):#gives us optimum uptake for different temperature curves
    return flux_at_all_theta(T, -4.5,   -1.48268657,  25., 10.86745456) #error 0.52,  0.21,  1.0, 1.8 
def grass_other_uptake(T):
    return flux_at_all_theta(T,  -2.33809598,  -1.27719641,  25., 14.75202332)#errors  0.44,  0.50,  1.0, 2.7
def grass_uptake(T,SW):
    opt_uptake_for_T = grass_opt_uptake(T)
    uptake_at_27_sw = grass_other_uptake(T)
    return flux_at_all_theta(SW,  opt_uptake_for_T,  uptake_at_27_sw,  26.9, grass_opt_sw()) #error 0.3

#BOREAL FOREST SOIL
#Based on "Siberian" soil data from van Diest and Kesselmeier (2008) and
#Field data from Hyytiala, Finland (Sun et al. 2018)
def boreal_opt_sw():
    return 12.5 #error 1.3
def boreal_opt_uptake(T):
    return flux_at_all_theta(T, -18.24779932, -12.,  35., 28.05488082) #error 2.3, 6.8, 2.5, 2.5
def boreal_other_uptake(T):
    return flux_at_all_theta(T, -5.89511395,  -3.76628476,  35., 28.05488082)#25.4438112) #error 1.1, 2.5, 3.5, 3.5 
def boreal_uptake(T,SW):
    opt_uptake_for_T = boreal_opt_uptake(T)
    uptake_at_20_sw = boreal_other_uptake(T)
    return flux_at_all_theta(SW,  opt_uptake_for_T,  uptake_at_20_sw,  19.3, boreal_opt_sw()) #error 0.6

#TEMPERATE FOREST SOIL
#Based on soil incubations from Willow Creek FLUXNET site (US-WCr) (Whelan et al. 2016)
#informed by field data from Harvard Forest (Wehr et al. 2017) and Wind River (Rastogi et al.2018 a,b)
def temperate_opt_sw():
    return 24.6 #average of 3 values, standard deviation 0.6
def temperate_opt_uptake(T):
    return -12.6 #Not enough data, upper limit declared to be the highest uptake recorded at Wind River, error 1.0
def temperate_other_uptake(T):
    return -0.17629655*T + 0.47914552 #r2 of 0.9 with 4 data points
def temperate_uptake(T,SW):
    opt_uptake_for_T = temperate_opt_uptake(T)
    uptake_at_51_sw = temperate_other_uptake(T)
    return flux_at_all_theta(SW,  opt_uptake_for_T,  uptake_at_51_sw,  51., temperate_opt_sw()) #error 21.6

#TROPICAL FOREST SOIL
#Based on soil incubations from Los Amigos Research Station, Peru (Whelan et al. 2016)
def tropical_opt_sw():
    return 24.6 #Incubations of *temperate* forest soil, average of 3 values, standard deviation 0.6
def tropical_opt_uptake(T):#gives us optimum uptake for different temperature curves
    return -2.7 #Not enough data, the upper limit is declared to be the highest uptake via soil incubations
def tropical_other_uptake(T):
    return -0.86 # 0.74 standard deviation error of 7 incubation measurements from 10 to 40 C
    # other uptake recorded at soil moisture 31 +/- 1.0 VWC
def tropical_uptake(T,SW):
    CONSTANT_PARAM = 31.
    opt_uptake_for_T = tropical_opt_uptake(T)
    uptake_at_31_sw = tropical_other_uptake(T)
    return flux_at_all_theta(SW,  opt_uptake_for_T,  uptake_at_31_sw,  31., tropical_opt_sw()) #error 1.0

#AGRICULTURAL SOIL
#Based on soil incubations from Bondville FLUXNET site (US-Bo1) (Whelan et al. 2016)
#and field data from Oklahoma wheatfield (Maseyk et al. 2014)
def ag_opt_sw(): 
    return 17.7 #average of 3 values, standard deviation 2.13
def ag_opt_uptake():#gives us optimum uptake for different temperature curves
    return -9.7 #error 1.8
    #the upper limit is declared to be the highest uptake recorded at Oklahoma
    #less the abiotic production function found by fitting obsevations
    #with <20 VWC
    #it appears that the max uptake increases over increasing temperature curves
    #evidence that we might be over-correcting for OCS production
    #though the agricultural soil OCS production equation is based on field data alone
def ag_other_uptake():
    return -5.36 #average of bondville soil incubations from 20.9 to 22.2 VWC
    #at various temperatures, less OCS production calculated above, error 0.78
def ag_uptake(T,SW):
    CONSTANT_PARAM = 22. #error 1.1
    opt_uptake_for_T = ag_opt_uptake()
    uptake_at_22_sw = ag_other_uptake()
    return flux_at_all_theta(SW,  opt_uptake_for_T,  uptake_at_22_sw,  22., ag_opt_sw()) #error 1.1

################################################################################
######            soil uptake and emission together                   ##########
################################################################################
#Now, add abiotic and biotic components together
#For each temperature bin/individual temperature
#Temperature in C (to make it obvious when we're assuming 0 soil fluxes)
#Soil moisture in %VWC, between 0 and 100 (more like 2 and 45)
#We have not been able to detect soil OCS uptake in wetland soils
#Wetland emissions are so large that it seems like uptake is a mild 
#uncertainty in production

def grass_soil_OCS(temperature, soilw):
    if np.isscalar(temperature) and temperature > 0.0:
        biotic_flux = grass_uptake(temperature,soilw)
        abiotic_flux = OCS_grass_production(temperature)
        return biotic_flux + abiotic_flux
    if np.isscalar(temperature)==False: 
        out = np.array(temperature)*0.0
        justright = np.where(temperature>0.0)
        biotic_flux = grass_uptake(temperature[justright],soilw[justright])
        abiotic_flux = OCS_grass_production(temperature[justright])
        out[justright] = biotic_flux + abiotic_flux
        return out
    else: 
        return 0.0
    
def bforest_soil_OCS(temperature, soilw, avg_T = 20, avg_SW = 20):
    if np.isscalar(temperature) and temperature > 0.0:
        biotic_flux = boreal_uptake(temperature,soilw)
        abiotic_flux = OCS_forest_production(temperature)
        return biotic_flux + abiotic_flux
    if np.isscalar(temperature)==False:
        out = np.array(temperature)*0.0
        justright = np.where(temperature>0.0)
        biotic_flux = boreal_uptake(temperature[justright],soilw[justright])
        abiotic_flux = OCS_forest_production(temperature[justright])
        out[justright] = biotic_flux + abiotic_flux
        return out
    else:
        return 0.0

def tforest_soil_OCS(temperature, soilw, avg_T = 20, avg_SW = 20):
    if np.isscalar(temperature) and temperature > 0.0:
        biotic_flux = temperate_uptake(temperature,soilw)
        abiotic_flux = OCS_forest_production(temperature)
        return biotic_flux + abiotic_flux
    if np.isscalar(temperature)==False:
        out = np.array(temperature)*0.0
        justright = np.where(temperature>0.0)
        biotic_flux = temperate_uptake(temperature[justright],soilw[justright])
        abiotic_flux = OCS_forest_production(temperature[justright])
        out[justright] = biotic_flux + abiotic_flux
        return out   
    else:
        return 0.0

def tropforest_soil_OCS(temperature, soilw, avg_T = 20, avg_SW = 20):
    if np.isscalar(temperature) and temperature > 0.0:
        biotic_flux = tropical_uptake(temperature,soilw)
        abiotic_flux = OCS_rainforest_production(temperature)
        return biotic_flux + abiotic_flux
    if np.isscalar(temperature)==False:
        out = np.array(temperature)*0.0
        justright = np.where(temperature>0.0)
        biotic_flux = tropical_uptake(temperature[justright],soilw[justright])
        abiotic_flux = OCS_rainforest_production(temperature[justright])
        out[justright] = biotic_flux + abiotic_flux
        return out   
    else:
        return 0.0
        
def ag_soil_OCS(temperature, soilw):
    if np.isscalar(temperature) and temperature > 0.0:
        biotic_flux = ag_uptake(temperature,soilw)
        abiotic_flux = OCS_ag_production(temperature)
        return biotic_flux + abiotic_flux
    if np.isscalar(temperature)==False:
        out = np.array(temperature)*0.0
        justright = np.where(temperature>0.0)
        biotic_flux = ag_uptake(temperature[justright],soilw[justright])
        abiotic_flux = OCS_ag_production(temperature[justright])
        out[justright] = biotic_flux + abiotic_flux
        return out
    else:
        return 0.0

def wetland_OCS(temperature, soilw=0):
    return OCS_wetland_production(temperature)

'''
Works Cited
Whelan et al. 2016 https://doi.org/10.5194/acp-16-3711-2016
Whelan and Rhew 2015 https://doi.org/10.1007/s10533-016-0207-7
Whelan et al. 2013 https://doi.org/10.1016/j.atmosenv.2013.02.048
Commane et al. https://doi.org/10.1073/pnas.1504131112
Liu et al., 2010 https://doi.org/10.5194/bg-7-753-2010
Maseyk et al., 2014 https://doi.org/10.1073/pnas.1319132111
Meredith et al., 2018 https://doi.org/10.3390/soilsystems2030037
Sun et al. 2018 https://doi.org/10.5194/acp-18-1363-2018
Van Diest and Kesselmeier 2008 https://doi.org/10.5194/bg-5-475-2008
'''