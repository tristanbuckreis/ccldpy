# ccldpy v2.0.1
# August 27, 2024
# Tristan E. Buckreis (tristanbuckreis@ucla.edu)
# https://github.com/tristanbuckreis/ccldpy/

###########################################################################################################################################

import numpy as np
import pandas as pd

###########################################################################################################################################

def WellsCoppersmith1994(magnitude, eqType):
    """
    Magnitude-Scaling relationships defined by Wells & Coppersmith (1994) for all style-of-faultings.
    
    Input Arguments:
    - magnitude = moment magnitude (Mw)
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    
    Output Arguments: 
    - A = rupture area (km^2)
    - AR = rupture aspect ratio (L/W) - median constrained to be >= 1.0 when extrapilated to small magnitudes
    - L = rupture along-strike length (km)
    - W = rupture down-dip width (km)
    """
    if eqType in ['intraslab','interface','stable']:
        raise ValueError("eqType: WellsCoppersmith1994 cannot be used for subduction-type or stable-continental events ('intraslab', 'interface', or 'stable')")
        return(np.nan, np.nan, np.nan, np.nan)

    elif eqType == 'crustal':
        # Rupture Area:
        a1 = -3.49; b1 = 0.91; s1 = 0.24
        A = 10**(a1 + b1*magnitude + np.random.normal(0, s1))
        # Rupture Length:
        a2 = -2.44; b2 =0.59; s2 = 0.16
        L = 10**(a2 + b2*magnitude + np.random.normal(0, s2))
        # Aspect Ratio:
        AR = L**2 / A
        # Check if the Aspect Ratio < 1.0, if so constrain and derive L & W:
        if AR < 1.0:
            # AR = 1.0
            s_cy08 = 0.16
            AR = np.random.normal(1, s_cy08)
            L = np.sqrt(A * AR)
            W = np.sqrt(A / AR)
        else:
            W = L / AR
        return(A, AR, L, W)

def Leonard2014(magnitude, mechanism, eqType):
    """
    Magnitude-Scaling relationships defined by Leonard (2014) for different tectonic regimes and style-of-faultings.
    
    Input Arguments:
    - magnitude = moment magnitude (Mw)
    - mechanism = known or preferred style-of-faulting; can be inferred from rake angle (e.g., Ancheta et al. 2013)
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    
    Output Arguments: 
    - A = rupture area (km^2)
    - AR = rupture aspect ratio (L/W) - median constrained to be >= 1.0 when extrapilated to small magnitudes
    - L = rupture along-strike length (km)
    - W = rupture down-dip width (km)
    """
    if eqType in ['intraslab','interface']:
        raise ValueError("eqType: Leonard2014 cannot be used for subduction-type events ('intraslab' or 'interface')")
        return(np.nan, np.nan, np.nan, np.nan)
        
    elif eqType == 'crustal':
        if mechanism == 'SS':
            # Rupture Area:
            a1 = 3.99; b1 = 1.00; s1 = 0.13
            A = 10**((magnitude - a1 - np.random.normal(0, s1))/b1)
            # Rupture Length:
            a2 = 4.170; b2 = 1.667; s2 = 0.19
            L = 10**((magnitude - a2 - np.random.normal(0, s2))/b2)
            if L > 45.0:
                a3 = 5.27; b3 = 1.000; 
                L = 10**((magnitude - a3 - np.random.normal(0, s2))/b3)
        elif mechanism in ['NM','RV']:
            # Rupture Area:
            a1 = 4.00; b1 = 1.00; s1 = 0.15
            A = 10**((magnitude - a1 - np.random.normal(0, s1))/b1)
            # Rupture Length:
            a2 = 4.000; b2 = 2.000; s2 = 0.23
            L = 10**((magnitude - a2 - np.random.normal(0, s2))/b2)
            if L > 5.4:
                a3 = 4.240; b3 = 1.667;
                L = 10**((magnitude - a3 - np.random.normal(0, s2))/b3)
        # Aspect Ratio:
        AR = L**2 / A
        # Check if the Aspect Ratio < 1.0, if so constrain and derive L & W:
        if AR < 1.0:
            # AR = 1.0
            s_cy08 = 0.16
            AR = np.random.normal(1, s_cy08)
            L = np.sqrt(A * AR)
            W = np.sqrt(A / AR)
        else:
            W = L / AR
        return(A, AR, L, W)
    
    elif eqType == 'stable':
        if mechanism == 'SS':
            # Rupture Area:
            a1 = 4.18; b1 = 1.00; s1 = 0.09
            A = 10**((magnitude - a1 - np.random.normal(0, s1))/b1)
            # Rupture Length:
            a2 = 4.250; b2 = 1.667; s2 = 0.18
            L = 10**((magnitude - a2 - np.random.normal(0, s2))/b2)
            if L > 60.0:
                a3 = 5.44; b3 = 1.000; 
                L = 10**((magnitude - a3 - np.random.normal(0, s2))/b3)
        elif mechanism in ['NM','RV']:
            # Rupture Area:
            a1 = 4.19; b1 = 1.00; s1 = 0.10
            A = 10**((magnitude - a1 - np.random.normal(0, s1))/b1)
            # Rupture Length:
            a2 = 4.320; b2 = 1.667; s2 = 0.19
            L = 10**((magnitude - a2 - np.random.normal(0, s2))/b2)
        # Aspect Ratio:
        AR = L**2 / A
        # Check if the Aspect Ratio < 1.0, if so constrain and derive L & W:
        if AR < 1.0:
            # AR = 1.0
            s_cy08 = 0.16
            AR = np.random.normal(1, s_cy08)
            L = np.sqrt(A * AR)
            W = np.sqrt(A / AR)
        else:
            W = L / AR
        return(A, AR, L, W)

def ThingbaijamEtAl2017(magnitude, mechanism, eqType):
    """
    Magnitude-Scaling relationships defined by Thingbaijam et al. (2017) for different tectonic regimes and style-of-faultings.
    
    Input Arguments:
    - magnitude = moment magnitude (Mw)
    - mechanism = known or preferred style-of-faulting; can be inferred from rake angle (e.g., Ancheta et al. 2013)
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    
    Output Arguments: 
    - A = rupture area (km^2)
    - AR = rupture aspect ratio (L/W) - median constrained to be >= 1.0 when extrapilated to small magnitudes
    - L = rupture along-strike length (km)
    - W = rupture down-dip width (km)
    """
    if eqType in ['intraslab','stable']:
        raise ValueError("eqType: ThingbaijamEtAl2017 cannot be used for instra-slab or stable-continental events ('intraslab' or 'stable')")
        return(np.nan, np.nan, np.nan, np.nan)
    
    elif eqType == 'crustal':
        if mechanism == 'SS':
            # Rupture Area:
            a1 = -3.486; b1 = 0.942; s1 = 0.184
            A = 10**(a1 + b1*magnitude + np.random.normal(0, s1))
            # Rupture Length:
            a2 = -2.943; b2 = 0.681; s2 = 0.151
            L = 10**(a2 + b2*magnitude + np.random.normal(0, s2))
        elif mechanism == 'NM':
            # Rupture Area:
            a1 = -2.551; b1 = 0.808; s1 = 0.181
            A = 10**(a1 + b1*magnitude + np.random.normal(0, s1))
            # Rupture Length:
            a2 = -1.722; b2 = 0.485; s2 = 0.128
            L = 10**(a2 + b2*magnitude + np.random.normal(0, s2))
        elif mechanism == 'RV':
            # Rupture Area:
            a1 = -4.362; b1 = 1.049; s1 = 0.121
            A = 10**(a1 + b1*magnitude + np.random.normal(0, s1))
            # Rupture Length:
            a2 = -2.693; b2 = 0.614; s2 = 0.083
            L = 10**(a2 + b2*magnitude + np.random.normal(0, s2))
        # Aspect Ratio:
        AR = L**2 / A
        # Check if the Aspect Ratio < 1.0, if so constrain and derive L & W:
        if AR < 1.0:
            # AR = 1.0
            s_cy08 = 0.16
            AR = np.random.normal(1, s_cy08)
            L = np.sqrt(A * AR)
            W = np.sqrt(A / AR)
        else:
            W = L / AR
        return(A, AR, L, W)
    
    elif eqType == 'interface':
        # Rupture Area:
        a1 = -3.292; b1 = 0.949; s1 = 0.150
        A = 10**(a1 + b1*magnitude + np.random.normal(0, s1))
        # Rupture Length:
        a2 = -2.412; b2 = 0.583; s2 = 0.107
        L = 10**(a2 + b2*magnitude + np.random.normal(0, s2))
        # Aspect Ratio:
        AR = L**2 / A
        # Check if the Aspect Ratio < 1.0, if so constrain and derive L & W:
        if AR < 1.0:
            # AR = 1.0
            s_cy08 = 0.16
            AR = np.random.normal(1, s_cy08)
            L = np.sqrt(A * AR)
            W = np.sqrt(A / AR)
        else:
            W = L / AR
        return(A, AR, L, W)

def ChiouYoungs2008(magnitude, mechanism, eqType, model):
    """
    Magnitude-Scaling aspect ratio relationship defined by Chiou & Youngs (2008) for different style-of-faultings
    in combination with other published magnitude-scaling area relationships.
    
    Input Arguments:
    - magnitude = moment magnitude (Mw)
    - mechanism = known or preferred style-of-faulting; can be inferred from rake angle (e.g., Ancheta et al. 2013)
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    - model = magntiude-scaling area relationship to use with Chiou & Youngs (2008) aspect ratio relation
        "WellsCoppersmith1994" = Wells & Coppersmith (1994)
        "Leonard2014" = Leonard (2014)
        "ThingbaijamEtAl2017" = Thingbaijam et al. (2017)
    
    Output Arguments: 
    - A = rupture area (km^2)
    - AR = rupture aspect ratio (L/W) - median constrained to be >= 1.0 when extrapilated to small magnitudes
    - L = rupture along-strike length (km)
    - W = rupture down-dip width (km)
    """
    if eqType in ['intraslab','interface','stable']:
        raise ValueError("eqType: ChiouYoungs2008 cannot be used for subduction-type or stable-continental events ('intraslab', 'interface', or 'stable')")
        return(np.nan, np.nan, np.nan, np.nan)
    
    elif eqType == 'crustal':
        if model == 'WellsCoppersmith1994':
            # Rupture Area:
            a1 = -3.49; b1 = 0.91; s1 = 0.24
            A = 10**(a1 + b1*magnitude + np.random.normal(0,s1))
        elif model == 'Leonard2014':
            if mechanism == 'SS':
                # Rupture Area:
                a1 = 3.99; b1 = 1.00; s1 = 0.13
                A = 10**((magnitude - a1 - np.random.normal(0,s1))/b1)
            elif mechanism in ['NM','RV']:
                # Rupture Area:
                a1 = 4.00; b1 = 1.00; s1 = 0.15
                A = 10**((magnitude - a1 - np.random.normal(0,s1))/b1)
        elif model == 'ThingbaijamEtAl2017':
            if mechanism == 'SS':
                # Rupture Area:
                a1 = -3.486; b1 = 0.942; s1 = 0.184
                A = 10**(a1 + b1*magnitude + np.random.normal(0,s1))
            elif mechanism == 'NM':
                # Rupture Area:
                a1 = -2.551; b1 = 0.808; s1 = 0.181
                A = 10**(a1 + b1*magnitude + np.random.normal(0,s1))
            elif mechanism == 'RV':
                # Rupture Area:
                a1 = -4.362; b1 = 1.049; s1 = 0.121
                A = 10**(a1 + b1*magnitude + np.random.normal(0,s1))
        # Aspect Ratio:
        s = 0.16
        if magnitude < 4.0:
            # AR = 1.0
            AR = np.random.normal(1, s)
        else:
            FNM = 0.0; FRV = 0.0
            if mechanism == 'NM':
                FNM = 1.0
            elif mechanism == 'RV':
                FRV = 1.0
            ass = 0.01752; anm = -0.00472; arv = -0.01099
            b = 3.097; 
            AR = 10**( (ass + anm*FNM + arv*FRV)*(magnitude - 4.0)**b + np.random.normal(0,s))
        # Compute Length and Width:
        L = np.sqrt(A * AR)
        W = np.sqrt(A / AR)
        return(A, AR, L, W)

def ContrerasEtAl2022(magnitude, eqType):
    """
    Magnitude-Scaling aspect ratio relationship defined by Contreras et al. (222) for subduction-type earthquakes.
    
    Input Arguments:
    - magnitude = moment magnitude (Mw)
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    
    Output Arguments: 
    - A = rupture area (km^2)
    - AR = rupture aspect ratio (L/W) - median constrained to be >= 1.0 when extrapilated to small magnitudes
    - L = rupture along-strike length (km)
    - W = rupture down-dip width (km)
    """  
    if eqType in ['crustal','stable']:
        raise ValueError("eqType: ContrerasEtAl2022 cannot be used for non-subduction-type events ('crustal' or 'stable')")
        return(np.nan, np.nan, np.nan, np.nan)
    
    elif eqType == 'interface':
        # Rupture Area:
        a1 = -3.8290; b1 = 1.0; s1 = 0.270
        A = 10.0**(a1 + b1*magnitude + np.random.normal(0,s1))
        # Aspect Ratio:
        if magnitude > 7.25:
            AR = 10.0**(0.2759 * (magnitude - 7.25) + np.random.normal(0,0.192)) # std = 0.441 in ln unit
        else:
            AR = 10.00**(0.0 + np.random.normal(0,0.0717)) # std = 0.165 in ln unit
        # Compute Length and Width:
        L = np.sqrt(A * AR)
        W = np.sqrt(A / AR)
        return(A, AR, L, W)
    
    elif eqType == 'intraslab':
        # Rupture Area:
        a1 = -3.251; b1 = 0.890; s1 = 0.184
        A = 10.0**(a1 + b1*magnitude + np.random.normal(0,s1))
        # Aspect Ratio:
        if magnitude > 6.5:
            AR = 10.0**(0.0938 * (magnitude - 6.5) + np.random.normal(0,0.164)) # std = 0.378in ln unit
        else:
            AR = 10.00**(0.0 + np.random.normal(0,0.104)) # std = 0.239 in ln unit
        # Compute Length and Width:
        L = np.sqrt(A * AR)
        W = np.sqrt(A / AR)
        return(A, AR, L, W)
    
def get_mechanism_based_on_rake(rake):
    """
    Infer style-of-faulting (mechanism) from rake angle (e.g., Ancheta et al. 2013)
    
    Input Arguments:
    rake = rake angle (degrees)
    
    Output Arguments:
    mechanism = indicator string for style-of-faulting (mechanism)
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    """
    if (-180 <= rake < -150) | (-30 <= rake < 30) | (150 <= rake <= 180):
        return('SS')
    elif (-120 <= rake < -60) | (-150 <= rake < -120) | (-60 <= rake < -30):
        return('NM')
    elif (60 <= rake < 120) | (30 <= rake < 60) | (120 <= rake < 150):
        return('RV')
    
def discrete(n,x,p,tp):
    for i in range(10):
        if (tp >= p[i]) and (tp <= p[i+1]):
            tx = x[i] + (x[i+1] - x[i])*(tp - p[i])/(p[i+1] - p[i])
            return(tx)
    return

def get_hyp_down_dip_position(eqtype, region):
    """
    Simulate hypocenter's down-dip relative position on the rupture surface based on Chiou & Youngs (2008)
    
    Input arguments:
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    - region = geographic region where the earthquake occured
        "japan"
        "chile"
        "other"
        
    Output Arguments:
    - fd = down-dip relative position of hypocenter on rupture surface [0, 1]
        0 ~ along the top-depth of the rupture
        1 ~ along the bottom-depth of the rupture
    """
    nxf = 11
    xdf = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00])
    ssd = np.array([0, 0.025, 0.05, 0.1, 0.175, 0.275, 0.4, 0.55, 0.7, 0.85, 1.00])
    iad = np.array([0.00, 0.012, 0.051, 0.139, 0.294, 0.500, 0.706, 0.861, 0.949, 0.988, 1.00]) 
    ied0 = np.array([0.000, 0.024, 0.085, 0.206, 0.389, 0.599, 0.783, 0.906, 0.969, 0.993, 1.00])
    ied1 = np.array([0.0, 0.002, 0.012, 0.044, 0.121, 0.262, 0.460, 0.671, 0.843, 0.950, 1.00])
    ied2 = np.array([0.00, 0.013, 0.053, 0.143, 0.297, 0.500, 0.703, 0.857, 0.947, 0.987, 1.00])
    if eqtype in ['crustal', 'stable']:
        fd = discrete(nxf, xdf, ssd, np.random.rand())
    elif eqtype == 'intraslab':
        fd = discrete(nxf, xdf, iad, np.random.rand())
    elif eqtype == 'interface':
        if region == "japan": 
            fd = discrete(nxf, xdf, ied0, np.random.rand())
        elif region == "chile": 
            fd = discrete(nxf, xdf, ied1, np.random.rand())
        elif region == "other":
            fd = discrete(nxf, xdf, ied2, np.random.rand())
    return(fd)

def get_hyp_along_strike_position(eqtype):
    """
    Simulate hypocenter's along-strike relative position on the rupture surface based on Chiou & Youngs (2008)
    
    Input arguments:
    - eqType = type of earthquake (tectonic regime)
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
        
    Output Arguments:
    - fl = along-strike relative position of hypocenter on rupture surface [0, 1]
        0 ~ along the "left" edge of the rupture
        1 ~ along the "right" edge of the rupture
    """
    nxf = 11
    xdf = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00])
    hypx = np.array([0, 0.05, 0.125, 0.225, 0.35, 0.5, 0.65, 0.775, 0.875, 0.95, 1.00])  
    hypxe = np.array([0.000, 0.007, 0.034, 0.112, 0.272, 0.500, 0.728, 0.888, 0.966, 0.993, 1.00])
    hypxa = np.array([0.00, 0.015, 0.057, 0.148, 0.301, 0.500, 0.699, 0.852, 0.943, 0.985, 1.00])
    if eqtype in ['crustal', 'stable']:
        fl = discrete(nxf, xdf, hypx, np.random.rand())
    elif eqtype == 'intraslab':
        fl = discrete(nxf, xdf, hypxa, np.random.rand())
    elif eqtype == 'interface':
        fl = discrete(nxf, xdf, hypxe, np.random.rand())
    return(fl)

def get_median_index(array):
    index = np.argpartition(array, len(array) // 2)[len(array) // 2]
    return(index)

def LL2XY(zlon, zlat, lon, lat):
    """
    Convert latitude/longitude to x/y coordinates on the basis of a spherical earth 
    (ignoring ellipsoidal effects) using the haversine formula.
    
    Reference: http://www.movable-type.co.uk/scripts/latlong.html
    
    Input Arguments:
    - zlon = centered longitude (degrees): corresponding to x/y coordinate of (0,0)
    - zlat = centered latitude (degrees): corresponding to x/y coordinate of (0,0)
    - lon = numpy array of longitudes (degrees)
    - lat = numpy array of latitudes (degrees)
    
    Output Arguments:
    -x = numpy array of x-coordinate
    -y = numpy array of y-coordinates
    
    """
    Rearth = 6371.0
    dtr = np.arcsin(1.0)/90
    lam1 = (zlon + 360)*dtr if zlon < 0 else zlon*dtr
    phi1 = zlat*dtr
        
    #lam2 = (lon + 360)*dtr if lon < 0 else lon*dtr
    if type(lon) == 'float':
        lam2 = (lon + 360)*dtr if lon < 0 else lon*dtr
    else:
        lam2 = np.where(lon < 0, (lon + 360)*dtr, lon*dtr)
        
    phi2 = lat*dtr
    dphi = phi2 - phi1
    dlam = lam2 - lam1
    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlam/2)**2
    c = 2*np.arctan2( np.sqrt(a), np.sqrt(1-a) )
    dist = Rearth*c
    if dist == 0:
        theta = 0
    else:
        theta = np.arctan2( np.sin(dlam)*np.cos(phi2), 
                            np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(dlam) )
    x = dist*np.sin(theta)
    y = dist*np.cos(theta)
        
    return(x,y)

def XY2LL(zlon, zlat, x, y):
    """
    Convert x/y coordinates to latitude/longitude on the basis of a spherical earth 
    (ignoring ellipsoidal effects) using the haversine formula.
    
    Reference: http://www.movable-type.co.uk/scripts/latlong.html
    
    Input Arguments:
    - zlon = centered longitude (degrees): corresponding to x/y coordinate of (0,0)
    - zlat = centered latitude (degrees): corresponding to x/y coordinate of (0,0)
    - x = numpy array of x-coordinate
    - y = numpy array of y-coordinates
    
    Output Arguments:
    - lon = numpy array of longitudes (degrees)
    - lat = numpy array of latitudes (degrees)
    
    """
    Rearth = 6371.0
    dtr = np.arcsin(1.0)/90
    lam1 = (zlon + 360)*dtr if zlon < 0 else zlon*dtr
    phi1 = zlat*dtr

    d = np.sqrt(x**2 + y**2)
    delta = d/Rearth
    theta = np.arctan2(x, y) 
    phi2 = np.arcsin( np.sin(phi1)*np.cos(delta) + np.cos(phi1)*np.sin(delta)*np.cos(theta) )
    lam2 = lam1 + np.arctan2( np.sin(theta)*np.sin(delta)*np.cos(phi1),
                              np.cos(delta) - np.sin(phi1)*np.sin(phi2) )
    lat = phi2/dtr
    lon = lam2/dtr
    if type(lon) == 'float':
        lon = lon - 360 if lon > 180 else lon
    else:
        lon = np.where(lon > 180, lon - 360, lon)
        
    return(lon, lat) 

def pointTriangleDistance(TRI_xyz, P_xyz):
    """
    Compute shortest 3D distance between a triangle and point (Eberly 1999)
    
    Reference: https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    
    Coded in Python by Meera Kota (UCLA)
    
    Input Arguments:
    - TRI_xyz = numpy arrays of 3D triangles in XY space [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]
    - P_xyz = numpy array of 3D points in XY space [[x1,y1,z1],...,[xn,yn,zn]] where n = number of points    
    """
    E0 = TRI_xyz[:, 1] - TRI_xyz[:, 0]
    E1 = TRI_xyz[:, 2] - TRI_xyz[:, 0]
    a = np.sum(np.multiply(E0, E0), axis=1)
    b = np.sum(np.multiply(E0, E1), axis=1)
    c = np.sum(np.multiply(E1, E1), axis=1)

    det = a * c - b * b
    det[det==0] = 1.0e-8

    D = TRI_xyz[:, 0] - P_xyz

    d = np.sum(np.multiply(E0, D), axis=2)
    e = np.sum(np.multiply(E1, D), axis=2)
    f = np.sum(np.multiply(D, D), axis=2)

    s = b * e - c * d
    t = b * d - a * e

    sqrdistance = np.empty((len(P_xyz), len(TRI_xyz)), dtype=float)

    # Region 4
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d < 0.0) & (-d >= a)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d < 0.0) & (-d < a)
    sqrdistance[cond] = (-d * d / a + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e < 0.0) & (-e >= c)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e < 0.0) & (-e < c)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 3
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e < 0.0) & (-e >= c)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e < 0.0) & (-e < c)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 5
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d < 0.0) & (-d >= a)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d < 0.0) & (-d < a)
    sqrdistance[cond] = (-d * d / a + f)[cond]

    # Region 0
    invDet = 1.0 / det
    stemp = s * invDet
    ttemp = t * invDet
    cond = (s + t <= det) & (s >= 0) & (t >= 0)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]

    # Region 2
    tmp0 = b + d
    tmp1 = c + e
    numer = tmp1 - tmp0
    denom = a - 2.0 * b + c
    denom[denom==0] = 1.0e-8
    stemp = numer / denom
    ttemp = 1.0 - stemp
    cond = (s + t > det) & (s < 0.0) & (tmp1 > tmp0) & (numer >= denom)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 > tmp0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 <= 0.0)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 > 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 > 0.0) & (e < 0.0)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 6
    tmp0 = b + e
    tmp1 = a + d
    numer = tmp1 - tmp0
    denom = a - 2.0 * b + c
    denom[denom==0] = 1.0e-8
    ttemp = numer / denom
    stemp = 1.0 - ttemp
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 > tmp0) & (numer >= denom)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 > tmp0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 <= 0)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 > 0) & (d >= 0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 > 0) & (d < 0)
    sqrdistance[cond] = (-d * d / a + f)[cond]

    # Region 1
    numer = c + e - b - d
    denom = a - 2.0 * b + c
    stemp = numer / denom
    ttemp = 1.0 - stemp
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer <= 0)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer > 0) & (numer >= denom)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer > 0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]

    # account for numerical round-off error
    sqrdistance[sqrdistance <= 0] = 0
    return np.sqrt(sqrdistance).T

def check_input_arguments(eqType, method, nsims, strike, dip, rake, strike2, dip2, rake2):
    # Check that required arguments are defined for each simulation method  -------------------------------------------------
    if method in ['A','D']:
        if None in np.array([strike, dip, rake]):
            raise ValueError("'strike', 'dip', and 'rake' must be defined when using 'method' = '%s'"%method)
    if method == 'B':
        if None in np.array([strike2, dip2, rake2]):
            raise ValueError("'strike2', 'dip2', and 'rake2' must be defined when using 'method' = 'B'")
    if method == 'C':
        if None in np.array([strike, dip, rake, strike2, dip2, rake2]):
            raise ValueError("'strike', 'dip', 'rake', 'strike2', 'dip2', and 'rake2' must be defined when using 'method' = 'C'")

    # Check that number of simulations are compatible with sepcified eqType  ------------------------------------------------
    if (eqType == 'crustal') & (nsims[6] > 0):
        raise ValueError("nsims[6] = %i: ContrerasEtAl2022 cannot be used for 'crustal' events"%nsims[6])
    elif (eqType == 'stable'):
        if (nsims[0] > 0):
            raise ValueError("nsims[0] = %i: WellsCoppersmith1994 cannot be used for 'stable' events"%nsims[0])
        elif (nsims[2] > 0):
            raise ValueError("nsims[2] = %i: ThingbaijamEtAl2017 cannot be used for 'stable' events"%nsims[2])
        elif (sum(nsims[3:6]) > 0):
            raise ValueError("nsims[3:6] = [%i,%i,%i]: ChiouYoungs2008 cannot be used for 'stable' events"%(nsims[3], nsims[4], nsims[5]))
        elif (nsims[6] > 0):
            raise ValueError("nsims[6] = %i: ContrerasEtAl2022 cannot be used for 'stable' events"%nsims[6])
    elif (eqType in ['interface','intraslab']):
        if (nsims[0] > 0):
            raise ValueError("nsims[0] = %i: WellsCoppersmith1994 cannot be used for '%s' events"%(nsims[0],eqType))
        elif (nsims[1] > 0):
            raise ValueError("nsims[1] = %i: Leonard2014 cannot be used for '%s' events"%(nsims[1],eqType))
        elif (sum(nsims[3:6]) > 0):
            raise ValueError("nsims[3:6] = [%i,%i,%i]: ChiouYoungs2008 cannot be used for '%s' events"%(nsims[3], nsims[4], nsims[5], eqType))

    # Check that there are only positive integers for the number of simulations  --------------------------------------------
    if np.any(np.array(nsims), where= np.array(nsims) < 0):
        raise ValueError("nsims = %s: can not specify a negative number of simulatiuons for any scaling relationship"%str(list(nsims)))

    # Check that the total number of simulations is odd and greater than zero  ----------------------------------------------
    if sum(nsims) % 2 == 0:
        if eqType == 'crustal':
            if nsims[0] > 0:
                nsims[0] = nsims[0] + 1
            elif nsims[1] > 0:
                nsims[1] = nsims[1] + 1
            elif nsims[2] > 0:
                nsims[2] = nsims[2] + 1
            elif nsims[3] > 0:
                nsims[3] = nsims[3] + 1
            elif nsims[4] > 0:
                nsims[4] = nsims[4] + 1
            elif nsims[5] > 0:
                nsims[5] = nsims[5] + 1
            else:
                raise ValueError("nsims = %s: must specify the number of simulations partitioned to each scaling relationship"%str(list(nsims)))
        if eqType == 'stable':
            if nsims[1] > 0:
                nsims[1] = nsims[1] + 1
            else:
                raise ValueError("nsims = %s: must specify the number of simulations partitioned to each scaling relationship"%str(list(nsims)))
        if eqType in ['interface','intraslab']:
            if nsims[6] > 0:
                nsims[6] = nsims[6] + 1
            elif nsims[2] > 0:
                nsims[2] = nsims[2] + 1
            else:
                raise ValueError("nsims = %s: must specify the number of simulations partitioned to each scaling relationship"%str(list(nsims)))

###########################################################################################################################################
                
v_WellsCoppersmith1994 = np.vectorize(WellsCoppersmith1994)
v_Leonard2014 = np.vectorize(Leonard2014)
v_ThingbaijamEtAl2017 = np.vectorize(ThingbaijamEtAl2017)
v_ChiouYoungs2008 = np.vectorize(ChiouYoungs2008)
v_ContrerasEtAl2022 = np.vectorize(ContrerasEtAl2022)
v_get_mechanism_based_on_rake = np.vectorize(get_mechanism_based_on_rake)
v_get_hyp_down_dip_position = np.vectorize(get_hyp_down_dip_position)
v_get_hyp_along_strike_position = np.vectorize(get_hyp_along_strike_position)
v_get_median_index = np.vectorize(get_median_index, signature='(n)->()')
v_LL2XY = np.vectorize(LL2XY, signature='(),(),(n),(n)->(n),(n)')
v_XY2LL = np.vectorize(XY2LL, signature='(),(),(n),(n)->(n),(n)')

###########################################################################################################################################

def get_rupture_surface_simulations(eqid,
                                    eqType, region,
                                    elat, elon, hypd,
                                    magnitude,
                                    method, nsims,
                                    mechanism = None,
                                    strike = None, dip = None, rake = None,
                                    strike2 = None, dip2 = None, rake2 = None):
    """
    Simulate earthquake rupture surface that minimizes the difference between the median distance of a pseudo-grid of sites 
    and a stochastic set of possible ruptures. Randomized fault attributes are dependent on the fault category. Fault scaling 
    models depend on earthquake type. Adapted from the CCLD routine originally coded in Frortran by Chiou and Youngs (2008).
    
    Input Arguments:
    - eqid = unique integer identifier for the event (to set the random seed)
    - eqType = type of earthquake:
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    - region = geographic region where the earthquake occured"
        "japan"
        "chile"
        "other"
    - elat = hypocenter latitude (degrees)
    - elon = hypocenter longiude (degrees)
    - hypd = hypocenter depth (km); positive into the ground
    - magnitude = earthquake moment magnitude (Mw)
    - method = code for rupture simulation constraints
        "A" = strike, dip, and rake for first nodal plane solution preferred 
              (optional "strike", "dip", and "rake" arguments are required)
              Warning: not recommended
        "B" = strike, dip, and rake for second nodal plane solution preferred
              (optional "strike2", "dip2", and "rake2" arguments are required)
              Warning: not recommended
        "C" = strike, dip, and rake are known for two nodal planes, and neither is preferred
              (optional "strike", "dip", "rake", "strike2", "dip2", and "rake2" arguments are required)
        "D" = One nodal plane solution for strike, dip, and rake; randomize the strike and dip
              (optional "strike", "dip", and "rake" arguments are required)
              Warning: not recommended
        "E" = No nodal plane solutions; randomize strike, dip, and rake
              (dip and rake are assigned based on faulting mechanism)
              (if optional "mechanism" argument is not speficied, simulations randomly assign one)
    - nsims = Number of simulations assigned to each M-scaling relationship. Total number of simulations should be odd.
        nsims[0] = Wells & Coppersmith (1994) - [recommended 334]
        nsims[1] = Leonard (2014) - [recommended 333]
        nsims[2] = Thingbaijam et al. (2017) [recommended 333]
        nsims[3] = Chiou & Youngs (2008) aspect ratio model with Wells & Coppersmith (1994) area relationship [recommended 111]
        nsims[4] = Chiou & Youngs (2008) aspect ratio model with Leonard (2014) area relationship [recommended 111]
        nsims[5] = Chiou & Youngs (2008) aspect ratio model with Thingbaijam et al. (2017) area relationship [recommended 111]
        nsims[6] = Contreras et al. (2022) - [recommended 333]
    - mechanim = known or preferred style-of-faulting [default None]
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    - strike = strike-angle (degrees) of the first nodal plane solution [default None]
    - dip = dip-angle (degrees) of the first nodal plane solution [default None]
    - rake = rake-angle (degrees) of the first nodal plane solution [default None]
    - strike2 = strike-angle (degrees) of the second nodal plane solution [default None]
    - dip2 = dip-angle (degrees) of the second nodal plane solution [default None]
    - rake2 = rake-angle (degrees) of the second nodal plane solution [default None]
    
    Output Arguments:
    - median_simulation = index of selected (median) rupture
    - along_strikes = numpy array of the hypocenter along-strike relative position [0, 1]
    - down_dips = numpy array of the hypocenter down-dip relative position [0, 1]
    - lon1s = numpy array of longitude (degrees) for the ULC rupture vertex
    - lat1s = numpy array of latitude (degrees) for the ULC rupture vertex
    - rpz1s = numpy array of depth (km) for the ULC rupture vertex
    - lon2s = numpy array of longitude (degrees) for the URC rupture vertex
    - lat2s = numpy array of latitude (degrees) for the URC rupture vertex
    - rpz2s = numpy array of depth (km) for the URC rupture vertex
    - lon3s = numpy array of longitude (degrees) for the LLC rupture vertex
    - lat3s = numpy array of latitude (degrees) for the LLC rupture vertex
    - rpz3s = numpy array of depth (km) for the LLC rupture vertex
    - lon4s = numpy array of longitude (degrees) for the LRC rupture vertex
    - lat4s = numpy array of latitude (degrees) for the LRC rupture vertex
    - rpz4s = numpy array of depth (km) for the LRC rupture vertex
    - strikes = numpy array of strike angle (degrees) 
    - dips = numpy array of dip angle (degrees)
    - rakes = numpy array of rake angle (degrees)
    - models = numpy array of magnitude-scaling relationship used for each realization
    - areas = numpy array of area (km^2)
    - aspect_ratios = numpy array of aspect ratio
    - lengths = numpy array of rupture length (km)
    - widths = numpy array of rupture width (km)
    - top_depths = numpy array of top depth (km) of rupture
    - bottom_depths = numpy array of bottom depth (km) of rupture
    """
    
    check_input_arguments(eqType, method, nsims, strike, dip, rake, strike2, dip2, rake2)
    try:
        strike = float(strike); dip = float(dip); rake = float(rake)
    except:
        pass
    try:
        strike2 = float(strike2); dip2 = float(dip2); rake2 = float(rake2)
    except:
        pass

    # Set random seed -------------------------------------------------------------------------------------------------------
    np.random.seed(seed=eqid)
    total_sims = np.sum(nsims)

    # Set earthquake parameters  --------------------------------------------------------------------------------------------
    if method == 'A':
        # Prefer the first nodal plane solution:
        strikes = np.repeat([strike], total_sims)
        dips = np.repeat([dip], total_sims)
        rakes = np.repeat([rake], total_sims)
        mechanisms = np.repeat(get_mechanism_based_on_rake(rake), total_sims)
    elif method == 'B':
        # Prefer the second nodal plane solution:
        strikes = np.repeat([strike2], total_sims)
        dips = np.repeat([dip2], total_sims)
        rakes = np.repeat([rake2], total_sims)
        mechanisms = np.repeat(get_mechanism_based_on_rake(rake2), total_sims)
    elif method == 'C':
        # Randomize between the two nodal plane solutions:
        nodal_planes = np.random.choice([1,2], total_sims)
        strikes = np.where(nodal_planes == 1, strike, strike2)
        dips = np.where(nodal_planes == 1, dip, dip2)
        rakes = np.where(nodal_planes == 1, rake, rake2)
        mechanisms = v_get_mechanism_based_on_rake(rakes)
    elif method == 'D':
        # Randomize strike and dip with a uniform distributions centered around the perscribed values
        strikes = np.random.uniform(strike - 30.0, strike + 30.0, total_sims)  # +/- 30 degrees
        strikes = np.where(strikes < 0, strikes + 360, strikes)                # strike >= 0 
        strikes = np.where(strikes >= 360.0, strikes - 360.0, strikes)         # strike < 360 
        dips = np.random.uniform(dip - 10.0, dip + 10.0, total_sims)           # +/- 10 degrees
        dips = np.where(dips < 10.0, 10.0, dips)                               # 10 <= dip
        dips = np.where(dips > 89.9999999, 89.9999999, dips) # dip <= 90 # cannot divide by zero
        rakes = np.repeat([rake], total_sims)
        mechanisms = v_get_mechanism_based_on_rake(rakes)
    elif method == 'E':
        # Randomize everything based on prescribed mechanism
        # if mechanism is not prescribed, then randomize everything
        if mechanism == None:
            mechanisms = np.random.choice(['SS','NM','RV'], total_sims)
        else:
            mechanisms = np.repeat(mechanism, total_sims)
        strikes = np.random.uniform(0.0, 360.0, total_sims)
        strikes = np.where(strikes >= 360.0, strikes - 360.0, strikes)    # 0 <= strike < 360
        rakes = np.where(mechanisms == 'SS', 0.0, mechanisms)
        rakes = np.where(rakes == 'NM', -90.0, rakes)
        rakes = np.where(rakes == 'RV', 90.0, rakes)
        dips = np.where(mechanisms == 'SS', 89.9999999, mechanisms) # cannot divide by zero
        dips = np.where(dips == 'NM', 55.0, dips)
        dips = np.where(dips == 'RV', 40.0, dips)
        dips = dips.astype(float)
    dips = np.where(dips > 89.9999999, 89.9999999, dips) # dip <= 90 # cannot divide by zero
    
    # Convert degrees to radians  -------------------------------------------------------------------------------------------
    strikes *= np.pi / 180.0
    dips *= np.pi / 180.0

    # Simulate hypocenter down-dip and along-strike position  ---------------------------------------------------------------
    down_dips = v_get_hyp_down_dip_position(np.repeat(eqType, total_sims), region)
    along_strikes = v_get_hyp_along_strike_position(np.repeat(eqType, total_sims))

    # Simulate rupture geometry using magnitude-scaling relationships  ------------------------------------------------------
    if nsims[0] > 0:
        nsim_WellsCoppersmith1994 = nsims[0]
        relationships_WellsCoppersmith1994 = v_WellsCoppersmith1994(np.repeat(magnitude, nsim_WellsCoppersmith1994), 
                                                                    np.repeat(eqType, nsim_WellsCoppersmith1994))
        areas_WellsCoppersmith1994 = relationships_WellsCoppersmith1994[0]
        aspect_ratios_WellsCoppersmith1994 = relationships_WellsCoppersmith1994[1]
        lengths_WellsCoppersmith1994 = relationships_WellsCoppersmith1994[2]
        widths_WellsCoppersmith1994 = relationships_WellsCoppersmith1994[3]
        models = np.repeat('WellsCoppersmith1994', nsim_WellsCoppersmith1994)
    else:
        areas_WellsCoppersmith1994, aspect_ratios_WellsCoppersmith1994 = np.array([]), np.array([])
        lengths_WellsCoppersmith1994, widths_WellsCoppersmith1994 = np.array([]), np.array([])
        models = np.array([])

    if nsims[1] > 0:
        nsim_Leonard2014 = nsims[1]
        mechanisms_Leonard2014 = mechanisms[np.sum(nsims[:1]): np.sum(nsims[:2])]
        relationships_Leonard2014 = v_Leonard2014(np.repeat(magnitude, nsim_Leonard2014), mechanisms_Leonard2014, np.repeat(eqType, nsim_Leonard2014))
        areas_Leonard2014 = relationships_Leonard2014[0]
        aspect_ratios_Leonard2014 = relationships_Leonard2014[1]
        lengths_Leonard2014 = relationships_Leonard2014[2]
        widths_Leonard2014 = relationships_Leonard2014[3]
        models = np.concatenate((models, np.repeat('Leonard2014', nsim_Leonard2014)))
    else:
        areas_Leonard2014, aspect_ratios_Leonard2014 = np.array([]), np.array([])
        lengths_Leonard2014, widths_Leonard2014 = np.array([]), np.array([])

    if nsims[2] > 0:
        nsim_ThingbaijamEtAl2017 = nsims[2]
        mechanisms_ThingbaijamEtAl2017 = mechanisms[np.sum(nsims[:2]): np.sum(nsims[:3])]
        relationships_ThingbaijamEtAl2017 = v_ThingbaijamEtAl2017(np.repeat(magnitude, nsim_ThingbaijamEtAl2017), 
                                                                  mechanisms_ThingbaijamEtAl2017, np.repeat(eqType, nsim_ThingbaijamEtAl2017))
        areas_ThingbaijamEtAl2017 = relationships_ThingbaijamEtAl2017[0]
        aspect_ratios_ThingbaijamEtAl2017 = relationships_ThingbaijamEtAl2017[1]
        lengths_ThingbaijamEtAl2017 = relationships_ThingbaijamEtAl2017[2]
        widths_ThingbaijamEtAl2017 = relationships_ThingbaijamEtAl2017[3]
        models = np.concatenate((models, np.repeat('ThingbaijamEtAl2017', nsim_ThingbaijamEtAl2017)))
    else:
        areas_ThingbaijamEtAl2017, aspect_ratios_ThingbaijamEtAl2017 = np.array([]), np.array([])
        lengths_ThingbaijamEtAl2017, widths_ThingbaijamEtAl2017 = np.array([]), np.array([])

    if nsims[3] > 0:
        nsim_ChiouYoungs2008_WellsCoppersmith1994 = nsims[3]
        mechanisms_ChiouYoungs2008_WellsCoppersmith1994 = mechanisms[np.sum(nsims[:3]): np.sum(nsims[:4])]
        relationships_ChiouYoungs2008_WellsCoppersmith1994 = v_ChiouYoungs2008(np.repeat(magnitude, nsim_ChiouYoungs2008_WellsCoppersmith1994), 
                                                                               mechanisms_ChiouYoungs2008_WellsCoppersmith1994,
                                                                               np.repeat(eqType, nsim_ChiouYoungs2008_WellsCoppersmith1994),
                                                                               'WellsCoppersmith1994')
        areas_ChiouYoungs2008_WellsCoppersmith1994 = relationships_ChiouYoungs2008_WellsCoppersmith1994[0]
        aspect_ratios_ChiouYoungs2008_WellsCoppersmith1994 = relationships_ChiouYoungs2008_WellsCoppersmith1994[1]
        lengths_ChiouYoungs2008_WellsCoppersmith1994 = relationships_ChiouYoungs2008_WellsCoppersmith1994[2]
        widths_ChiouYoungs2008_WellsCoppersmith1994 = relationships_ChiouYoungs2008_WellsCoppersmith1994[3]
        models = np.concatenate((models, np.repeat('ChiouYoungs2008_WellsCoppersmith1994', nsim_ChiouYoungs2008_WellsCoppersmith1994)))
    else:
        areas_ChiouYoungs2008_WellsCoppersmith1994, aspect_ratios_ChiouYoungs2008_WellsCoppersmith1994 = np.array([]), np.array([])
        lengths_ChiouYoungs2008_WellsCoppersmith1994, widths_ChiouYoungs2008_WellsCoppersmith1994 = np.array([]), np.array([])

    if nsims[4] > 0:
        nsim_ChiouYoungs2008_Leonard2014 = nsims[4]
        mechanisms_ChiouYoungs2008_Leonard2014 = mechanisms[np.sum(nsims[:4]): np.sum(nsims[:5])]
        relationships_ChiouYoungs2008_Leonard2014 = v_ChiouYoungs2008(np.repeat(magnitude, nsim_ChiouYoungs2008_Leonard2014), 
                                                                      mechanisms_ChiouYoungs2008_Leonard2014,
                                                                      np.repeat(eqType, nsim_ChiouYoungs2008_Leonard2014),
                                                                      'Leonard2014')
        areas_ChiouYoungs2008_Leonard2014 = relationships_ChiouYoungs2008_Leonard2014[0]
        aspect_ratios_ChiouYoungs2008_Leonard2014 = relationships_ChiouYoungs2008_Leonard2014[1]
        lengths_ChiouYoungs2008_Leonard2014 = relationships_ChiouYoungs2008_Leonard2014[2]
        widths_ChiouYoungs2008_Leonard2014 = relationships_ChiouYoungs2008_Leonard2014[3]
        models = np.concatenate((models, np.repeat('ChiouYoungs2008_Leonard2014', nsim_ChiouYoungs2008_Leonard2014)))
    else:
        areas_ChiouYoungs2008_Leonard2014, aspect_ratios_ChiouYoungs2008_Leonard2014 = np.array([]), np.array([])
        lengths_ChiouYoungs2008_Leonard2014, widths_ChiouYoungs2008_Leonard2014 = np.array([]), np.array([])

    if nsims[5] > 0:
        nsim_ChiouYoungs2008_ThingbaijamEtAl2017 = nsims[5]
        mechanisms_ChiouYoungs2008_ThingbaijamEtAl2017 = mechanisms[np.sum(nsims[:5]): np.sum(nsims[:6])]
        relationships_ChiouYoungs2008_ThingbaijamEtAl2017 = v_ChiouYoungs2008(np.repeat(magnitude, nsim_ChiouYoungs2008_ThingbaijamEtAl2017), 
                                                                              mechanisms_ChiouYoungs2008_ThingbaijamEtAl2017,
                                                                              np.repeat(eqType, nsim_ChiouYoungs2008_ThingbaijamEtAl2017),
                                                                              'ThingbaijamEtAl2017')
        areas_ChiouYoungs2008_ThingbaijamEtAl2017 = relationships_ChiouYoungs2008_ThingbaijamEtAl2017[0]
        aspect_ratios_ChiouYoungs2008_ThingbaijamEtAl2017 = relationships_ChiouYoungs2008_ThingbaijamEtAl2017[1]
        lengths_ChiouYoungs2008_ThingbaijamEtAl2017 = relationships_ChiouYoungs2008_ThingbaijamEtAl2017[2]
        widths_ChiouYoungs2008_ThingbaijamEtAl2017 = relationships_ChiouYoungs2008_ThingbaijamEtAl2017[3]
        models = np.concatenate((models, np.repeat('ChiouYoungs2008_ThingbaijamEtAl2017', nsim_ChiouYoungs2008_ThingbaijamEtAl2017)))
    else:
        areas_ChiouYoungs2008_ThingbaijamEtAl2017, aspect_ratios_ChiouYoungs2008_ThingbaijamEtAl2017 = np.array([]), np.array([])
        lengths_ChiouYoungs2008_ThingbaijamEtAl2017, widths_ChiouYoungs2008_ThingbaijamEtAl2017 = np.array([]), np.array([])
        
    if nsims[6] > 0:
        nsim_ContrerasEtAl2022 = nsims[6]
        mechanisms_ContrerasEtAl2022 = mechanisms[np.sum(nsims[:6]): np.sum(nsims[:7])]
        relationships_ContrerasEtAl2022 = v_ContrerasEtAl2022(np.repeat(magnitude, nsim_ContrerasEtAl2022), np.repeat(eqType, nsim_ContrerasEtAl2022))
        areas_ContrerasEtAl2022 = relationships_ContrerasEtAl2022[0]
        aspect_ratios_ContrerasEtAl2022 = relationships_ContrerasEtAl2022[1]
        lengths_ContrerasEtAl2022 = relationships_ContrerasEtAl2022[2]
        widths_ContrerasEtAl2022 = relationships_ContrerasEtAl2022[3]
        models = np.concatenate((models, np.repeat('ContrerasEtAl2022', nsim_ContrerasEtAl2022)))
    else:
        areas_ContrerasEtAl2022, aspect_ratios_ContrerasEtAl2022 = np.array([]), np.array([])
        lengths_ContrerasEtAl2022, widths_ContrerasEtAl2022 = np.array([]), np.array([])

    areas = np.concatenate((areas_WellsCoppersmith1994, areas_Leonard2014, areas_ThingbaijamEtAl2017, areas_ChiouYoungs2008_WellsCoppersmith1994,
                            areas_ChiouYoungs2008_Leonard2014, areas_ChiouYoungs2008_ThingbaijamEtAl2017,
                            areas_ContrerasEtAl2022))
    aspect_ratios = np.concatenate((aspect_ratios_WellsCoppersmith1994, aspect_ratios_Leonard2014, aspect_ratios_ThingbaijamEtAl2017, 
                                    aspect_ratios_ChiouYoungs2008_WellsCoppersmith1994, aspect_ratios_ChiouYoungs2008_Leonard2014, 
                                    aspect_ratios_ChiouYoungs2008_ThingbaijamEtAl2017,
                                    aspect_ratios_ContrerasEtAl2022))
    lengths = np.concatenate((lengths_WellsCoppersmith1994, lengths_Leonard2014, lengths_ThingbaijamEtAl2017, lengths_ChiouYoungs2008_WellsCoppersmith1994,
                              lengths_ChiouYoungs2008_Leonard2014, lengths_ChiouYoungs2008_ThingbaijamEtAl2017,
                              lengths_ContrerasEtAl2022))
    widths = np.concatenate((widths_WellsCoppersmith1994, widths_Leonard2014, widths_ThingbaijamEtAl2017, widths_ChiouYoungs2008_WellsCoppersmith1994,
                             widths_ChiouYoungs2008_Leonard2014, widths_ChiouYoungs2008_ThingbaijamEtAl2017,
                             widths_ContrerasEtAl2022))
    
    # Compute geometric parameters  -----------------------------------------------------------------------------------------
    xf = np.sin(strikes) * lengths * along_strikes                   
    yf = np.cos(strikes) * lengths * along_strikes                  
    xb = np.sin(strikes + np.pi) * lengths * (1.0 - along_strikes)
    yb = np.cos(strikes + np.pi) * lengths * (1.0 - along_strikes)
    rwh = widths * np.cos(dips)
    rwv = widths * np.sin(dips)
    top_depths = hypd - rwv * down_dips 
    down_dips = np.where(top_depths < 0.0, hypd / rwv, down_dips) # make sure hypocenter is in the ground
    top_depths = np.where(top_depths < 0.0, 0.0, top_depths)      # top_depth > 0.0
    bottom_depths = top_depths + rwv

    # Compute top points of the rupture surface 1 = ULC & 2 = URC  ----------------------------------------------------------
    rpx1s = xf + np.sin(strikes - np.arcsin(1.0)) * rwh * down_dips
    rpy1s = yf + np.cos(strikes - np.arcsin(1.0)) * rwh * down_dips
    rpz1s = top_depths

    rpx2s = xb + np.sin(strikes - np.arcsin(1.0)) * rwh * down_dips
    rpy2s = yb + np.cos(strikes - np.arcsin(1.0)) * rwh * down_dips
    rpz2s = top_depths

    # Compute bottom points of the rupture surface 3 = LRC & 4 = LLC  -------------------------------------------------------
    rpx3s = xf + np.sin(strikes + np.arcsin(1.0)) * rwh * (1.0 - down_dips)
    rpy3s = yf + np.cos(strikes + np.arcsin(1.0)) * rwh * (1.0 - down_dips)
    rpz3s = bottom_depths

    rpx4s = xb + np.sin(strikes + np.arcsin(1.0)) * rwh * (1.0 - down_dips)
    rpy4s = yb + np.cos(strikes + np.arcsin(1.0)) * rwh * (1.0 - down_dips)
    rpz4s = bottom_depths

    # Create a psudo-grid of stations  --------------------------------------------------------------------------------------
    r = np.array(list(np.arange(2, 20, 2)) + list(np.arange(25, 55, 5)) + list(np.arange(60, 110, 10)) + list(np.arange(125, 325, 25)))
    theta = np.pi*np.linspace(0, 360, 25)[:-1]/180
    psx = np.multiply(np.repeat(r.reshape(1,len(r)), len(theta), axis=0), np.sin(theta.reshape(len(theta),1))).flatten()
    psy = np.multiply(np.repeat(r.reshape(1,len(r)), len(theta), axis=0), np.cos(theta.reshape(len(theta),1))).flatten()
    
    # Compute distances to each psudo-station  ------------------------------------------------------------------------------
    SITES = np.array([[[x,y,0]] for x, y in zip(psx, psy)])
    TRI_123 = np.array([[[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]] for x1, x2, x3, y1, y2, y3, z1, z2, z3 in zip(rpx1s, rpx2s, rpx3s, rpy1s, rpy2s, rpy3s, rpz1s, rpz2s, rpz3s)])
    TRI_234 = np.array([[[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]] for x1, x2, x3, y1, y2, y3, z1, z2, z3 in zip(rpx2s, rpx4s, rpx3s, rpy2s, rpy4s, rpy3s, rpz2s, rpz4s, rpz3s)])
    rrups1 = pointTriangleDistance(TRI_123, SITES)
    rrups2 = pointTriangleDistance(TRI_234, SITES)
    rrups = np.minimum(rrups1, rrups2)

    # Locate the median Rrup for each site  ---------------------------------------------------------------------------------
    median_rrups = np.median(rrups, axis=0)
    median_rrup_indices = v_get_median_index(np.transpose(rrups))

    # Locate the simulation that minimizes the difference between median Rrup and Rrup across all sites  --------------------
    sum_of_difference_squared = np.sum((median_rrups - rrups)**2, axis=1)
    median_simulation = np.argmin(sum_of_difference_squared)

    # Convert simulated faults xy to ll using haversine formulation  --------------------------------------------------------
    lon1s, lat1s = v_XY2LL(elon, elat, rpx1s, rpy1s)
    lon2s, lat2s = v_XY2LL(elon, elat, rpx2s, rpy2s)
    lon3s, lat3s = v_XY2LL(elon, elat, rpx3s, rpy3s)
    lon4s, lat4s = v_XY2LL(elon, elat, rpx4s, rpy4s)

    # Convert radians to degrees  -------------------------------------------------------------------------------------------
    strikes = np.round( strikes * 180 / np.pi, 6) 
    dips = np.round( dips * 180 / np.pi, 6) 

    # Return everything  ----------------------------------------------------------------------------------------------------
    return(median_simulation,
           along_strikes, down_dips,
           lon1s, lat1s, rpz1s,
           lon2s, lat2s, rpz2s,
           lon3s, lat3s, rpz3s,
           lon4s, lat4s, rpz4s,
           strikes, dips, rakes,
           models, areas, aspect_ratios, lengths, widths, top_depths, bottom_depths)

def simulate_rupture_surface(eqid,
                             eqType, region,
                             elat, elon, hypd,
                             magnitude,
                             method, nsims,
                             mechanism = None,
                             strike = None, dip = None, rake = None,
                             strike2 = None, dip2 = None, rake2 = None):
    
    """
    Simulate earthquake rupture surface that minimizes the difference between the median distance of a pseudo-grid of sites 
    and a stochastic set of possible ruptures. Randomized fault attributes are dependent on the fault category. Fault scaling 
    models depend on earthquake type. Adapted from the CCLD routine originally coded in Frortran by Chiou and Youngs (2008).
    
    Input Arguments:
    - eqid = unique integer identifier for the event (to set the random seed)
    - eqType = type of earthquake:
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    - region = geographic region where the earthquake occured"
        "japan"
        "chile"
        "other"
    - elat = hypocenter latitude (degrees)
    - elon = hypocenter longiude (degrees)
    - hypd = hypocenter depth (km); positive into the ground
    - magnitude = earthquake moment magnitude (Mw)
    - method = code for rupture simulation constraints
        "A" = strike, dip, and rake for first nodal plane solution preferred 
              (optional "strike", "dip", and "rake" arguments are required)
              Warning: not recommended
        "B" = strike, dip, and rake for second nodal plane solution preferred
              (optional "strike2", "dip2", and "rake2" arguments are required)
              Warning: not recommended
        "C" = strike, dip, and rake are known for two nodal planes, and neither is preferred
              (optional "strike", "dip", "rake", "strike2", "dip2", and "rake2" arguments are required)
        "D" = One nodal plane solution for strike, dip, and rake; randomize the strike and dip
              (optional "strike", "dip", and "rake" arguments are required)
              Warning: not recommended
        "E" = No nodal plane solutions; randomize strike, dip, and rake
              (dip and rake are assigned based on faulting mechanism)
              (if optional "mechanism" argument is not speficied, simulations randomly assign one)
    - nsims = Number of simulations assigned to each M-scaling relationship. Total number of simulations should be odd.
        nsims[0] = Wells & Coppersmith (1994) - [recommended 334]
        nsims[1] = Leonard (2014) - [recommended 333]
        nsims[2] = Thingbaijam et al. (2017) [recommended 333]
        nsims[3] = Chiou & Youngs (2008) aspect ratio model with Wells & Coppersmith (1994) area relationship [recommended 111]
        nsims[4] = Chiou & Youngs (2008) aspect ratio model with Leonard (2014) area relationship [recommended 111]
        nsims[5] = Chiou & Youngs (2008) aspect ratio model with Thingbaijam et al. (2017) area relationship [recommended 111]
        nsims[6] = Contreras et al. (2022) - [recommended 333]
    - mechanim = known or preferred style-of-faulting [default None]
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    - strike = strike-angle (degrees) of the first nodal plane solution [default None]
    - dip = dip-angle (degrees) of the first nodal plane solution [default None]
    - rake = rake-angle (degrees) of the first nodal plane solution [default None]
    - strike2 = strike-angle (degrees) of the second nodal plane solution [default None]
    - dip2 = dip-angle (degrees) of the second nodal plane solution [default None]
    - rake2 = rake-angle (degrees) of the second nodal plane solution [default None]
    
    Output Objects:
    - SIMULATIONS = pandas DataFrame object containing all simulated rupture surfaces
    - SELECTED = pandas DataFrame object containing the selected rupture surface and statistics
    """
    
    (median_simulation,
     along_strikes, down_dips,
     lon1s, lat1s, rpz1s,
     lon2s, lat2s, rpz2s,
     lon3s, lat3s, rpz3s,
     lon4s, lat4s, rpz4s,
     strikes, dips, rakes,
     models, areas, aspect_ratios, lengths, widths, 
     top_depths, bottom_depths) = get_rupture_surface_simulations(eqid, eqType, region, elat, elon, hypd, magnitude, method, nsims,
                                                                  mechanism, strike, dip, rake, strike2, dip2, rake2)
    
    # Compute the statistics of the simulation results ----------------------------------------------------------------------
    areas_avg = 10**np.mean(np.log10(areas))
    areas_std = 10**np.std(np.log10(areas))
    aspect_ratios_avg = 10**np.mean(np.log10(aspect_ratios))
    aspect_ratios_std = 10**np.std(np.log10(aspect_ratios))
    lengths_avg = 10**np.mean(np.log10(lengths))
    lengths_std = 10**np.std(np.log10(lengths))
    widths_avg = 10**np.mean(np.log10(widths))
    widths_std = 10**np.std(np.log10(widths))

    max_top_depths = np.max(top_depths)
    min_top_depths = np.min(top_depths)
    max_bottom_depths = np.max(bottom_depths)
    min_bottom_depths = np.min(bottom_depths)
    
    # Construct pandas DataFrame objects to package the results -------------------------------------------------------------
    total_sims = np.sum(nsims)
    SIMULATIONS = pd.DataFrame({'Simulation':np.array(range(total_sims))+1,
                                'EQID':np.repeat(eqid, total_sims),'Magnitude':np.repeat(magnitude, total_sims),
                                'Hypocenter Longitude':np.repeat(elon, total_sims), 'Hypocenter Latitude':np.repeat(elat, total_sims),
                                'Hypocenter Depth (km)':np.repeat(hypd, total_sims),
                                'Hypocenter Along-Strike Position':along_strikes,'Hypocenter Down-Dip Position':down_dips,
                                'ULC Longitude':lon1s,'ULC Latitude':lat1s,'ULC Depth (km)':rpz1s,
                                'URC Longitude':lon2s,'URC Latitude':lat2s,'URC Depth (km)':rpz2s,
                                'LRC Longitude':lon4s,'LRC Latitude':lat4s,'LRC Depth (km)':rpz4s,
                                'LLC Longitude':lon3s,'LLC Latitude':lat3s,'LLC Depth (km)':rpz3s,
                                'Strike':strikes,'Dip':dips,'Rake':rakes,
                                'Scaling Relation':models,
                                'Area (km^2)':areas,'Aspect Ratio':aspect_ratios,
                                'Rupture Length (km)':lengths,'Rupture Width (km)':widths,
                                'Rupture Top Depth (km)':top_depths,'Rupture Bottom Depth (km)':bottom_depths})
    
    SELECTED = pd.DataFrame({'Median Simulation':[median_simulation+1], 'EQID':[eqid], 'Magnitude':[magnitude],
                             'Hypocenter Longitude':[elon], 'Hypocenter Latitude':[elat], 'Hypocenter Depth (km)':[hypd],
                             'Hypocenter Along-Strike Position':[along_strikes[median_simulation]],'Hypocenter Down-Dip Position':[down_dips[median_simulation]],
                             'ULC Longitude':[lon1s[median_simulation]],'ULC Latitude':[lat1s[median_simulation]],'ULC Depth (km)':[rpz1s[median_simulation]],
                             'URC Longitude':[lon2s[median_simulation]],'URC Latitude':[lat2s[median_simulation]],'URC Depth (km)':[rpz2s[median_simulation]],
                             'LRC Longitude':[lon4s[median_simulation]],'LRC Latitude':[lat4s[median_simulation]],'LRC Depth (km)':[rpz4s[median_simulation]],
                             'LLC Longitude':[lon3s[median_simulation]],'LLC Latitude':[lat3s[median_simulation]],'LLC Depth (km)':[rpz3s[median_simulation]],
                             'Strike':[strikes[median_simulation]],'Dip':[dips[median_simulation]],'Rake':[rakes[median_simulation]],
                             'Scaling Relation':[models[median_simulation]],
                             'Area (km^2)':[areas[median_simulation]],'Average Area (km^2)':[areas_avg],'Area Standard Deviation (km^2)':[areas_std],
                             'Aspect Ratio':[aspect_ratios[median_simulation]],'Average Aspect Ratio':[aspect_ratios_avg],'Aspect Ratio Standard Deviation':[aspect_ratios_std],
                             'Length (km)':[lengths[median_simulation]],'Average Length (km)':[lengths_avg],'Length Standard Deviation (km)':[lengths_std],
                             'Width (km)':[widths[median_simulation]],'Average Width (km)':[widths_avg],'Width Standard Deviation (km)':[widths_std],
                             'Rupture Top Depth (km)':[top_depths[median_simulation]],'Rupture Bottom Depth (km)':[bottom_depths[median_simulation]]})
    
    return(SIMULATIONS, SELECTED)
