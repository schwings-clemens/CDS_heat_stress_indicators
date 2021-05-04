import sys
import numpy as np
import xarray as xr


#Check temperature units and convert temperature if required
#Input:
# - T: temperature in K
# - conv: Unit to which T in K should be converted to
#Output:
# - T_out: Converted temperature
def convert_T(T, conv):
    
    #Check that T is in Kelvin
    check = np.mean(T)<150
    if check==True:
        print('Temperature is not given in Kelvin!')
        sys.exit(1)
    
    #Convert T if required
    if conv=='Celsius':
        T_out = T - 273.15
    elif conv=='Fahrenheit':
        T_out = 9/5 * (T - 273.15) + 32
    else:
        T_out = T
    
    return T_out
    
    
# Calculate vapour pressure and relative humidity from specific humidity
#Input:
# - huss: specific humidity in kg/kg
# - p: pressure in Pa
# - T: temperature in K
#Output:
# - e: vapour pressure in Pa
# - RH: relative humidity
def get_humidity(huss, p, T):

    T_C = convert_T(T, 'Celsius')

    M_H2O = 18.01528/1000 # kg/mol
    M_air = 28.964/1000 # kg/mol for dry air

    # Calculate vapour pressure according to August-Roche-Magnus equation (in Pa)
    e = huss * p * M_air / M_H2O

    # Saturation vapor pressure (in Pa)
    e_s = 610.94 * np.exp(17.625 * T_C / (T_C + 243.04) )
    
    # Relative humidty
    RH = 100 * e / e_s

    #Set RH larger then 100 to 100
    RH = RH.where(RH<=100, 100)
    
    return e, RH


#Wet-bulb globe temperature
#Input:
# - T: temperature in K
# - RH: relative humidity
# - p: pressure in Pa
#Output:
# - WBGT: wet-bulb globe temperature in °C
def WBGT(T, RH, p):

    #Convert T to degCelsius
    T_C = convert_T(T, 'Celsius')
    
    #Calculate wet bulb temperature 
    T_w = WBT_DavJon(T, RH, p)

    #Wet-bulb globe temperature indoors (Knutson 2016, Epstein 2006)
    WBGT = 0.7 * T_w + 0.3 * T_C
    
    return WBGT


#Heat index according to NOAA (Rothfusz 1990, Steadman 1979)
#Input:
# - T: temperature in K
# - RH: relative humidity
#Output:
# - HI_C: NOAA heat index in °C
def HI_NOAA(T, RH):

    #Convert T to Fahrenheit
    T_F = convert_T(T, 'Fahrenheit')

    #Calculate heat index
    c    = [-42.379, 2.04901523, 10.14333127, -0.22475541, -0.00683783, -0.05481717, 0.00122874, 0.00085282, -0.00000199]
    HI_F = c[0] + c[1]*T_F + c[2]*RH + c[3]*T_F*RH +c[4]*T_F*T_F + c[5]*RH*RH + \
           c[6]*T_F*T_F*RH + c[7]*T_F*RH*RH + c[8]*T_F*T_F*RH*RH

    #Adjustments to HI
    radicand = (17 - np.abs(T_F - 95)) / 17
    radicand = xr.where(radicand<0, np.NaN, radicand)
    HI_adj1 =  (13 - RH) / 4 * np.sqrt(radicand)
    HI_adj2 =  (RH - 85) * (87 - T_F) / 50
    HI_simple = 0.5 * (T_F + 61.0 + ((T_F - 68.0) * 1.2) + (RH * 0.094))    
    HI_F = xr.where((RH<13) & (T_F>80) & (T_F<112), HI_F - HI_adj1, HI_F) #Adjustment 1
    HI_F = xr.where((RH>85) & (T_F>80) & (T_F<87), HI_F + HI_adj2, HI_F)  #Adjustment 2
    HI_F = xr.where(HI_simple<80, HI_simple, HI_F) #Adjustment 3

    #Convert to degCelsius
    HI_C = (HI_F - 32) * 5/9

    return HI_C


#Humidex (Masterson and Richardson 1979)
#Input:
# - T: temperature in K
# - e: vapour pressure in Pa
#Output:
# - Humidex: Humidex in °C
def Humidex(T, e):

    #Convert T to Celsius
    T_C = convert_T(T, 'Celsius')

    Humidex = T_C + 5/9 * (e/100 - 10)

    return Humidex


#Wet bulb temperature (Davies-Jones 2008, Bolton 1980, Dunne 2013)
#Input:
# - T: temperature in K
# - RH: relative humidity
# - p: pressure in Pa
#Output:
# - WBT: wet-bulb temperature in °C
def WBT_DavJon(T, RH, p):
    
    C    = 273.15 # K
    lamb = 3.504  # lambda = c_pd/R_d, Davies-Jones (2008)
    
    #Saturation vapor pressure (Dunne 2013, Bolton 1980)
    e_s = np.exp(-2991.2729/T**2 - 6017.0128/T + 18.87643854 - 0.028354721 * T + 1.7838301e-5 * T**2 - \
                  8.4150417e-10 * T**3 + 4.4412543e-13 * T**4 + 2.858487 * np.log(T))  

    #Saturation water vapor mixing ratio
    w_s = 621.97 * e_s / (p - e_s)
    
    #Water vapor mixing ratio and vapor pressure
    w = RH/100 * w_s    
    e = RH/100 * e_s
    
    # Temperature at condensation level
    T_L = (1 / (T - 55) - np.log(RH/100) / 2840)**(-1) + 55

    #Equivalent potential temperature (simpler formula)
    theta_E = T * (100000 / p)**(0.2854 * (1 - 0.28e-3 * w)) * np.exp( (3.376 / T_L - 0.00254) * w * (1 + 0.81e-3 * w))

    #Wet bulb temperature
    WBT = 45.114 - 51.489 * (theta_E/C)**(-lamb)

    #Set data type to float32
    WBT = WBT.astype('float32')
    
    return WBT


# Universal Thermal Climate Index (from comfort_models.py on github; https://gist.github.com/chriswmackey)
#Input:
# - T: temperature in K
# - e: vapour pressure in Pa
#Output:
# - UTCI: Universal Thermal Climate Index in °C
def UTCI(T, e):

    Pa = e / 1000 #convert to kPa
    Ta = convert_T(T, 'Celsius')
    va = 0.5 # m/s
    D_Tmrt = 0 # Assume that mean radiant temperature = near-surface air temperature
    
    #UTCI calculation split into chunks to avoid memory overload
    UTCI1  = get_UTCI1(Ta, Pa, va, D_Tmrt)
    UTCI2  = get_UTCI2(Ta, Pa, va, D_Tmrt)
    UTCI3  = get_UTCI3(Ta, Pa, va, D_Tmrt)
    UTCI4  = get_UTCI4(Ta, Pa, va, D_Tmrt)
    UTCI5  = get_UTCI5(Ta, Pa, va, D_Tmrt)
    UTCI6  = get_UTCI6(Ta, Pa, va, D_Tmrt)
    UTCI7  = get_UTCI7(Ta, Pa, va, D_Tmrt)
    UTCI8  = get_UTCI8(Ta, Pa, va, D_Tmrt)
    UTCI9  = get_UTCI9(Ta, Pa, va, D_Tmrt)
    UTCI10 = get_UTCI10(Ta, Pa, va, D_Tmrt)
    UTCI11 = get_UTCI11(Ta, Pa, va, D_Tmrt)
    UTCI12 = get_UTCI12(Ta, Pa, va, D_Tmrt)
    UTCI13 = get_UTCI13(Ta, Pa, va, D_Tmrt)
    
    #Sum up all UTCI terms
    UTCI = UTCI1 + UTCI2 + UTCI3 + UTCI4 + UTCI5 + UTCI6 + UTCI7 + UTCI8 + UTCI9 + UTCI10 + UTCI11 + UTCI12 + UTCI13

    return UTCI

#Chunk 1 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI1: Universal Thermal Climate Index chunk 1 in °C
def get_UTCI1(Ta, Pa, va, D_Tmrt):
    UTCI1 = (Ta +
           (0.607562052) +
           (-0.0227712343) * Ta +
           (8.06470249 * (10 ** (-4))) * Ta**2 +
           (-1.54271372 * (10 ** (-4)))* Ta**3 +
           (-3.24651735 * (10 ** (-6)))* Ta**4 +
           (7.32602852 * (10 ** (-8))) * Ta**5 +
           (1.35959073 * (10 ** (-9))) * Ta**6 +
           (-2.25836520) * va +
           (0.0880326035) * Ta * va +
           (0.00216844454) * Ta**2 * va +
           (-1.53347087 * (10 ** (-5))) * Ta**3 * va +
           (-5.72983704 * (10 ** (-7))) * Ta**4 * va +
           (-2.55090145 * (10 ** (-9))) * Ta**5 * va)
    return UTCI1.load()
    
#Chunk 2 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI2: Universal Thermal Climate Index chunk 2 in °C
def get_UTCI2(Ta, Pa, va, D_Tmrt):
    UTCI2 = ((-0.751269505) * va**2 +
           (-0.00408350271) * Ta * va**2 +
           (-5.21670675 * (10 ** (-5))) * Ta**2 * va**2 +
           (1.94544667 * (10 ** (-6))) * Ta**3 * va**2 +
           (1.14099531 * (10 ** (-8))) * Ta**4 * va**2 +
           (0.158137256) * va**3 +
           (-6.57263143 * (10 ** (-5))) * Ta * va**3 +
           (2.22697524 * (10 ** (-7))) * Ta**2 * va**3 +
           (-4.16117031 * (10 ** (-8))) * Ta**3 * va**3 +
           (-0.0127762753) * va**4 +
           (9.66891875 * (10 ** (-6))) * Ta * va**4 +
           (2.52785852 * (10 ** (-9))) * Ta**2 * va**4 +
           (4.56306672 * (10 ** (-4))) * va**5 +
           (-1.74202546 * (10 ** (-7))) * Ta * va**5 +
           (-5.91491269 * (10 ** (-6))) * va**6)
    return UTCI2.load()

#Chunk 3 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI3: Universal Thermal Climate Index chunk 3 in °C
def get_UTCI3(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI3 = 0
    else:
        UTCI3 = ((0.398374029) * D_Tmrt +
               (1.83945314 * (10 ** (-4))) * Ta * D_Tmrt +
               (-1.73754510 * (10 ** (-4))) * Ta**2 * D_Tmrt +
               (-7.60781159 * (10 ** (-7))) * Ta**3 * D_Tmrt +
               (3.77830287 * (10 ** (-8))) * Ta**4 * D_Tmrt +
               (5.43079673 * (10 ** (-10))) * Ta**5 * D_Tmrt +
               (-0.0200518269) * va * D_Tmrt +
               (8.92859837 * (10 ** (-4))) * Ta * va * D_Tmrt +
               (3.45433048 * (10 ** (-6))) * Ta**2 * va * D_Tmrt +
               (-3.77925774 * (10 ** (-7))) * Ta**3 * va * D_Tmrt +
               (-1.69699377 * (10 ** (-9))) * Ta**4 * va * D_Tmrt +
               (1.69992415 * (10 ** (-4))) * va**2 * D_Tmrt +
               (-4.99204314 * (10 ** (-5))) * Ta * va**2 * D_Tmrt +
               (2.47417178 * (10 ** (-7))) * Ta**2 * va**2 * D_Tmrt +
               (1.07596466 * (10 ** (-8))) * Ta**3 * va**2 * D_Tmrt +
               (8.49242932 * (10 ** (-5))) * va**3 * D_Tmrt +
               (1.35191328 * (10 ** (-6))) * Ta * va**3 * D_Tmrt +
               (-6.21531254 * (10 ** (-9))) * Ta**2 * va**3 * D_Tmrt +
               (-4.99410301 * (10 ** (-6))) * va**4 * D_Tmrt +
               (-1.89489258 * (10 ** (-8))) * Ta * va**4 * D_Tmrt +
               (8.15300114 * (10 ** (-8))) * va**5 * D_Tmrt)
        UTCI3.load()
    return UTCI3

#Chunk 4 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI4: Universal Thermal Climate Index chunk 4 in °C
def get_UTCI4(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI4 = 0
    else:
        UTCI4 = ((7.55043090 * (10 ** (-4))) * D_Tmrt**2 +
               (-5.65095215 * (10 ** (-5))) * Ta * D_Tmrt**2 +
               (-4.52166564 * (10 ** (-7))) * Ta**2 * D_Tmrt**2 +
               (2.46688878 * (10 ** (-8))) * Ta**3 * D_Tmrt**2 +
               (2.42674348 * (10 ** (-10))) * Ta**4 * D_Tmrt**2 +
               (1.54547250 * (10 ** (-4))) * va * D_Tmrt**2 +
               (5.24110970 * (10 ** (-6))) * Ta * va * D_Tmrt**2 +
               (-8.75874982 * (10 ** (-8))) * Ta**2 * va * D_Tmrt**2 +
               (-1.50743064 * (10 ** (-9))) * Ta**3 * va * D_Tmrt**2 +
               (-1.56236307 * (10 ** (-5))) * va**2 * D_Tmrt**2 +
               (-1.33895614 * (10 ** (-7))) * Ta * va**2 * D_Tmrt**2 +
               (2.49709824 * (10 ** (-9))) * Ta**2 * va**2 * D_Tmrt**2 +
               (6.51711721 * (10 ** (-7))) * va**3 * D_Tmrt**2 +
               (1.94960053 * (10 ** (-9))) * Ta * va**3 * D_Tmrt**2 +
               (-1.00361113 * (10 ** (-8))) * va**4 * D_Tmrt**2)
        UTCI4.load()
    return UTCI4

#Chunk 5 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI5: Universal Thermal Climate Index chunk 5 in °C
def get_UTCI5(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI5 = 0
    else:
        UTCI5 = ((-1.21206673 * (10 ** (-5))) * D_Tmrt**3 +
               (-2.18203660 * (10 ** (-7))) * Ta * D_Tmrt**3 +
               (7.51269482 * (10 ** (-9))) * Ta**2 * D_Tmrt**3 +
               (9.79063848 * (10 ** (-11))) * Ta**3 * D_Tmrt**3 +
               (1.25006734 * (10 ** (-6))) * va * D_Tmrt**3 +
               (-1.81584736 * (10 ** (-9))) * Ta * va * D_Tmrt**3 +
               (-3.52197671 * (10 ** (-10))) * Ta**2 * va * D_Tmrt**3 +
               (-3.36514630 * (10 ** (-8))) * va**2 * D_Tmrt**3 +
               (1.35908359 * (10 ** (-10))) * Ta * va**2 * D_Tmrt**3 +
               (4.17032620 * (10 ** (-10))) * va**3 * D_Tmrt**3 +
               (-1.30369025 * (10 ** (-9))) * D_Tmrt**4 +
               (4.13908461 * (10 ** (-10))) * Ta * D_Tmrt**4 +
               (9.22652254 * (10 ** (-12))) * Ta**2 * D_Tmrt**4 +
               (-5.08220384 * (10 ** (-9))) * va * D_Tmrt**4 +
               (-2.24730961 * (10 ** (-11))) * Ta * va * D_Tmrt**4 +
               (1.17139133 * (10 ** (-10))) * va**2 * D_Tmrt**4 +
               (6.62154879 * (10 ** (-10))) * D_Tmrt**5 +
               (4.03863260 * (10 ** (-13))) * Ta * D_Tmrt**5 +
               (1.95087203 * (10 ** (-12))) * va * D_Tmrt**5 +
               (-4.73602469 * (10 ** (-12))) * D_Tmrt**6)
        UTCI5.load()
    return UTCI5

#Chunk 6 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI6: Universal Thermal Climate Index chunk 6 in °C
def get_UTCI6(Ta, Pa, va, D_Tmrt):
    UTCI6 = ((5.12733497) * Pa +
           (-0.312788561) * Ta * Pa +
           (-0.0196701861) * Ta**2 * Pa +
           (9.99690870 * (10 ** (-4))) * Ta**3 * Pa +
           (9.51738512 * (10 ** (-6))) * Ta**4 * Pa +
           (-4.66426341 * (10 ** (-7))) * Ta**5 * Pa +
           (0.548050612) * va * Pa +
           (-0.00330552823) * Ta * va * Pa +
           (-0.00164119440) * Ta**2 * va * Pa +
           (-5.16670694 * (10 ** (-6))) * Ta**3 * va * Pa +
           (9.52692432 * (10 ** (-7))) * Ta**4 * va * Pa +
           (-0.0429223622) * va**2 * Pa +
           (0.00500845667) * Ta * va**2 * Pa +
           (1.00601257 * (10 ** (-6))) * Ta**2 * va**2 * Pa +
           (-1.81748644 * (10 ** (-6))) * Ta**3 * va**2 * Pa +
           (-1.25813502 * (10 ** (-3))) * va**3 * Pa +
           (-1.79330391 * (10 ** (-4))) * Ta * va**3 * Pa +
           (2.34994441 * (10 ** (-6))) * Ta**2 * va**3 * Pa +
           (1.29735808 * (10 ** (-4))) * va**4 * Pa +
           (1.29064870 * (10 ** (-6))) * Ta * va**4 * Pa +
           (-2.28558686 * (10 ** (-6))) * va**5 * Pa)
    return UTCI6.load()
    
#Chunk 7 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI7: Universal Thermal Climate Index chunk 7 in °C
def get_UTCI7(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI7 = 0
    else:
        UTCI7 = ((-0.0369476348) * D_Tmrt * Pa +
               (0.00162325322) * Ta * D_Tmrt * Pa +
               (-3.14279680 * (10 ** (-5))) * Ta**2 * D_Tmrt * Pa +
               (2.59835559 * (10 ** (-6))) * Ta**3 * D_Tmrt * Pa +
               (-4.77136523 * (10 ** (-8))) * Ta**4 * D_Tmrt * Pa +
               (8.64203390 * (10 ** (-3))) * va * D_Tmrt * Pa +
               (-6.87405181 * (10 ** (-4))) * Ta * va * D_Tmrt * Pa +
               (-9.13863872 * (10 ** (-6))) * Ta**2 * va * D_Tmrt * Pa +
               (5.15916806 * (10 ** (-7))) * Ta**3 * va * D_Tmrt * Pa +
               (-3.59217476 * (10 ** (-5))) * va**2 * D_Tmrt * Pa +
               (3.28696511 * (10 ** (-5))) * Ta * va**2 * D_Tmrt * Pa +
               (-7.10542454 * (10 ** (-7))) * Ta**2 * va**2 * D_Tmrt * Pa +
               (-1.24382300 * (10 ** (-5))) * va**3 * D_Tmrt * Pa +
               (-7.38584400 * (10 ** (-9))) * Ta * va**3 * D_Tmrt * Pa +
               (2.20609296 * (10 ** (-7))) * va**4 * D_Tmrt * Pa +
               (-7.32469180 * (10 ** (-4))) * D_Tmrt**2 * Pa)
        UTCI7.load()
    return UTCI7

#Chunk 8 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI8: Universal Thermal Climate Index chunk 8 in °C
def get_UTCI8(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI8 = 0
    else:
        UTCI8 = ((-1.87381964 * (10 ** (-5))) * Ta * D_Tmrt**2 * Pa +
               (4.80925239 * (10 ** (-6))) * Ta**2 * D_Tmrt**2 * Pa +
               (-8.75492040 * (10 ** (-8))) * Ta**3 * D_Tmrt**2 * Pa +
               (2.77862930 * (10 ** (-5))) * va * D_Tmrt**2 * Pa +
               (-5.06004592 * (10 ** (-6))) * Ta * va * D_Tmrt**2 * Pa +
               (1.14325367 * (10 ** (-7))) * Ta**2 * va * D_Tmrt**2 * Pa +
               (2.53016723 * (10 ** (-6))) * va**2 * D_Tmrt**2 * Pa +
               (-1.72857035 * (10 ** (-8))) * Ta * va**2 * D_Tmrt**2 * Pa +
               (-3.95079398 * (10 ** (-8))) * va**3 * D_Tmrt**2 * Pa +
               (-3.59413173 * (10 ** (-7))) * D_Tmrt**3 * Pa +
               (7.04388046 * (10 ** (-7))) * Ta * D_Tmrt**3 * Pa +
               (-1.89309167 * (10 ** (-8))) * Ta**2 * D_Tmrt**3 * Pa +
               (-4.79768731 * (10 ** (-7))) * va * D_Tmrt**3 * Pa +
               (7.96079978 * (10 ** (-9))) * Ta * va * D_Tmrt**3 * Pa +
               (1.62897058 * (10 ** (-9))) * va**2 * D_Tmrt**3 * Pa +
               (3.94367674 * (10 ** (-8))) * D_Tmrt**4 * Pa +
               (-1.18566247 * (10 ** (-9))) * Ta * D_Tmrt**4 * Pa +
               (3.34678041 * (10 ** (-10))) * va * D_Tmrt**4 * Pa +
               (-1.15606447 * (10 ** (-10))) * D_Tmrt**5 * Pa)
        UTCI8.load()
    return UTCI8
            
#Chunk 9 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI9: Universal Thermal Climate Index chunk 9 in °C
def get_UTCI9(Ta, Pa, va, D_Tmrt):
    UTCI9 = ((-2.80626406) * Pa**2 +
           (0.548712484) * Ta * Pa**2 +
           (-0.00399428410) * Ta**2 * Pa**2 +
           (-9.54009191 * (10 ** (-4))) * Ta**3 * Pa**2 +
           (1.93090978 * (10 ** (-5))) * Ta**4 * Pa**2 +
           (-0.308806365) * va * Pa**2 +
           (0.0116952364) * Ta * va * Pa**2 +
           (4.95271903 * (10 ** (-4))) * Ta**2 * va * Pa**2 +
           (-1.90710882 * (10 ** (-5))) * Ta**3 * va * Pa**2 +
           (0.00210787756) * va**2 * Pa**2 +
           (-6.98445738 * (10 ** (-4))) * Ta * va**2 * Pa**2 +
           (2.30109073 * (10 ** (-5))) * Ta**2 * va**2 * Pa**2 +
           (4.17856590 * (10 ** (-4))) * va**3 * Pa**2 +
           (-1.27043871 * (10 ** (-5))) * Ta * va**3 * Pa**2 +
           (-3.04620472 * (10 ** (-6))) * va**4 * Pa**2)
    return UTCI9.load()
        
#Chunk 10 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI10: Universal Thermal Climate Index chunk 10 in °C
def get_UTCI10(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI10 = 0
    else:
        UTCI10 = ((0.0514507424) * D_Tmrt * Pa**2 +
               (-0.00432510997) * Ta * D_Tmrt * Pa**2 +
               (8.99281156 * (10 ** (-5))) * Ta**2 * D_Tmrt * Pa**2 +
               (-7.14663943 * (10 ** (-7))) * Ta**3 * D_Tmrt * Pa**2 +
               (-2.66016305 * (10 ** (-4))) * va * D_Tmrt * Pa**2 +
               (2.63789586 * (10 ** (-4))) * Ta * va * D_Tmrt * Pa**2 +
               (-7.01199003 * (10 ** (-6))) * Ta**2 * va * D_Tmrt * Pa**2 +
               (-1.06823306 * (10 ** (-4))) * va**2 * D_Tmrt * Pa**2 +
               (3.61341136 * (10 ** (-6))) * Ta * va**2 * D_Tmrt * Pa**2 +
               (2.29748967 * (10 ** (-7))) * va**3 * D_Tmrt * Pa**2 +
               (3.04788893 * (10 ** (-4))) * D_Tmrt**2 * Pa**2 +
               (-6.42070836 * (10 ** (-5))) * Ta * D_Tmrt**2 * Pa**2 +
               (1.16257971 * (10 ** (-6))) * Ta**2 * D_Tmrt**2 * Pa**2 +
               (7.68023384 * (10 ** (-6))) * va * D_Tmrt**2 * Pa**2 +
               (-5.47446896 * (10 ** (-7))) * Ta * va * D_Tmrt**2 * Pa**2 +
               (-3.59937910 * (10 ** (-8))) * va**2 * D_Tmrt**2 * Pa**2 +
               (-4.36497725 * (10 ** (-6))) * D_Tmrt**3 * Pa**2 +
               (1.68737969 * (10 ** (-7))) * Ta * D_Tmrt**3 * Pa**2 +
               (2.67489271 * (10 ** (-8))) * va * D_Tmrt**3 * Pa**2 +
               (3.23926897 * (10 ** (-9))) * D_Tmrt**4 * Pa**2)
        UTCI10.load()
    return UTCI10
            
#Chunk 11 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI11: Universal Thermal Climate Index chunk 11 in °C
def get_UTCI11(Ta, Pa, va, D_Tmrt):
    UTCI11 = ((-0.0353874123) * Pa**3 +
           (-0.221201190) * Ta * Pa**3 +
           (0.0155126038) * Ta**2 * Pa**3 +
           (-2.63917279 * (10 ** (-4))) * Ta**3 * Pa**3 +
           (0.0453433455) * va * Pa**3 +
           (-0.00432943862) * Ta * va * Pa**3 +
           (1.45389826 * (10 ** (-4))) * Ta**2 * va * Pa**3 +
           (2.17508610 * (10 ** (-4))) * va**2 * Pa**3 +
           (-6.66724702 * (10 ** (-5))) * Ta * va**2 * Pa**3 +
           (3.33217140 * (10 ** (-5))) * va**3 * Pa**3)
    return UTCI11.load()
        
#Chunk 12 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI12: Universal Thermal Climate Index chunk 12 in °C
def get_UTCI12(Ta, Pa, va, D_Tmrt):
    if D_Tmrt==0:
        UTCI12 = 0
    else:
        UTCI12 = ((-0.00226921615) * D_Tmrt * Pa**3 +
               (3.80261982 * (10 ** (-4))) * Ta * D_Tmrt * Pa**3 +
               (-5.45314314 * (10 ** (-9))) * Ta**2 * D_Tmrt * Pa**3 +
               (-7.96355448 * (10 ** (-4))) * va * D_Tmrt * Pa**3 +
               (2.53458034 * (10 ** (-5))) * Ta * va * D_Tmrt * Pa**3 +
               (-6.31223658 * (10 ** (-6))) * va**2 * D_Tmrt * Pa**3 +
               (3.02122035 * (10 ** (-4))) * D_Tmrt**2 * Pa**3 +
               (-4.77403547 * (10 ** (-6))) * Ta * D_Tmrt**2 * Pa**3 +
               (1.73825715 * (10 ** (-6))) * va * D_Tmrt**2 * Pa**3 +
               (-4.09087898 * (10 ** (-7))) * D_Tmrt**3 * Pa**3)
        UTCI12.load()
    return UTCI12

#Chunk 13 of UTCI calculation
#Input:
# - Ta: temperature in °C
# - Pa: pressure in kPa
# - va: wind speed in m/s
# - D_Tmrt: difference mean radiant temperature - air temperature
#Output:
# - UTCI13: Universal Thermal Climate Index chunk 13 in °C
def get_UTCI13(Ta, Pa, va, D_Tmrt):
    UTCI13 = ((0.614155345) * Pa**4 +
           (-0.0616755931) * Ta * Pa**4 +
           (0.00133374846) * Ta**2 * Pa**4 +
           (0.00355375387) * va * Pa**4 +
           (-5.13027851 * (10 ** (-4))) * Ta * va * Pa**4 +
           (1.02449757 * (10 ** (-4))) * va**2 * Pa**4 +
           (-0.00148526421) * D_Tmrt * Pa**4 +
           (-4.11469183 * (10 ** (-5))) * Ta * D_Tmrt * Pa**4 +
           (-6.80434415 * (10 ** (-6))) * va * D_Tmrt * Pa**4 +
           (-9.77675906 * (10 ** (-6))) * D_Tmrt**2 * Pa**4 +
           (0.0882773108) * Pa**5 +
           (-0.00301859306) * Ta * Pa**5 +
           (0.00104452989) * va * Pa**5 +
           (2.47090539 * (10 ** (-4))) * D_Tmrt * Pa**5 +
           (0.00148348065) * Pa**6)    
    return UTCI13.load()
