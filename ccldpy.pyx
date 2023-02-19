import random
import numpy as np
import pandas as pd
cimport numpy as np
cimport cython
from cython cimport boundscheck, wraparound
from libc.math cimport sqrt, log, cos, sin, asin, tan, atan2, exp, log10

@boundscheck(False)
@wraparound(False)
def gasdev(IDUM):
    cdef float R=2.0
    cdef float V1, V2, FAC, GSET, GASDEV
    while R >= 1:
        V1 = 2.*random.random()-1.
        V2 = 2.*random.random()-1.
        R = V1**2 + V2**2
    FAC = sqrt(-2.*log(R)/R)
    GSET = V1*FAC
    GASDEV = V2*FAC
    return(random.choice((GSET, GASDEV)))

@boundscheck(False)
@wraparound(False)
def discrete(n,x,p,tp):
    cdef int i
    for i in range(10):
        if (tp >= p[i]) and (tp <= p[i+1]):
            tx = x[i] + (x[i+1] - x[i])*(tp - p[i])/(p[i+1] - p[i])
            return(tx)
    return

@boundscheck(False)
@wraparound(False)
def LL2XY(zlon, zlat, lon, lat):
    """
    
    Convert latitude/longitude to x/y coordinates on the basis of a spherical earth 
    (ignoring ellipsoidal effects) using the haversine formula.
    
    Reference: http://www.movable-type.co.uk/scripts/latlong.html
    
    Input Arguments:
    zlon -- centered longitude (degrees): corresponding to x/y coordinate of (0,0)
    zlat -- centered latitude (degrees): corresponding to x/y coordinate of (0,0)
    lon -- numpy array of longitudes (degrees)
    lat -- numpy array of latitudes (degrees)
    
    Output Arguments:
    x -- numpy array of x-coordinate
    y -- numpy array of y-coordinates
    
    """
    cdef int i, Rearth
    cdef int npts=len(lon)
    cdef float dtr, lam1, phi1, lam2, phi2, dhp, dlam, a, c, dist, theta
    cdef double[:] x = np.empty(npts, dtype='float')
    cdef double[:] y = np.empty(npts, dtype='float')
    Rearth = 6371
    dtr = asin(1.0)/90
    lam1 = (zlon+360)*dtr if zlon<0 else zlon*dtr
    phi1 = zlat*dtr
    for i in range(npts):
        lam2 = (lon[i]+360)*dtr if lon[i]<0 else lon[i]*dtr
        phi2 = lat[i]*dtr
        dphi = phi2 - phi1
        dlam = lam2 - lam1
        a = sin(dphi/2)**2 + cos(phi1)*cos(phi2)*sin(dlam/2)**2
        c = 2*atan2( sqrt(a), sqrt(1-a) )
        dist = Rearth*c
        if dist == 0:
            theta = 0
        else:
            theta = atan2( sin(dlam)*cos(phi2), 
                           cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(dlam) )
        x[i] = dist*sin(theta)
        y[i] = dist*cos(theta)
    return(x,y)

@boundscheck(False)
@wraparound(False)
def XY2LL(zlon, zlat, x, y):
    """
    
    Convert x/y coordinates to latitude/longitude on the basis of a spherical earth 
    (ignoring ellipsoidal effects) using the haversine formula.
    
    Reference: http://www.movable-type.co.uk/scripts/latlong.html
    
    Input Arguments:
    zlon -- centered longitude (degrees): corresponding to x/y coordinate of (0,0)
    zlat -- centered latitude (degrees): corresponding to x/y coordinate of (0,0)
    x -- numpy array of x-coordinate
    y -- numpy array of y-coordinates
    
    Output Arguments:
    lon -- numpy array of longitudes (degrees)
    lat -- numpy array of latitudes (degrees)
    
    """
    cdef int i, Rearth
    cdef int npts=len(x)
    cdef float dtr, lam1, phi1, d, delta, theta, phi2, lam2
    cdef double[:] lon = np.empty(npts, dtype='float')
    cdef double[:] lat = np.empty(npts, dtype='float')
    Rearth = 6371
    dtr = asin(1.0)/90
    lam1 = (zlon+360)*dtr if zlon<0 else zlon*dtr
    phi1 = zlat*dtr
    for i in range(npts):
        d = sqrt(x[i]**2 + y[i]**2)
        delta = d/Rearth
        theta = atan2(x[i], y[i]) 
        phi2 = asin( sin(phi1)*cos(delta) + cos(phi1)*sin(delta)*cos(theta) )
        lam2 = lam1 + atan2( sin(theta)*sin(delta)*cos(phi1),
                             cos(delta) - sin(phi1)*sin(phi2) )
        lat[i] = phi2/dtr
        lon[i] = lam2/dtr
        lon[i] = lon[i]-360 if lon[i]>180 else lon[i]
    return(lon, lat)     

@boundscheck(False)
@wraparound(False)
def DEMDT(P):
    """
    
    Compute minimum distance to trapizoid based on David Ebery (2008)
    originally coded in FORTRAN by Brian Chiou: modified to place site at (0,0)
    coded in PYTHON by Tristan Buckreis
    
    Input Argument:
    P -- numpy array representing trapezoid:
         P[0][n] -- x-coordinates
         P[1][n] -- y-coordinates
         P[2][n] -- z-coordinate
         where n = 0:3 (four points)
         
    Output Arguments:
    dclst  -- closest distance between point (0,0) and trapezoid
    P_clst -- numpy array representing closest point on the trapezoid 
              P_clst[0] -- x-coordinate
              P_clst[1] -- y-coordinate
              P_clst[2] -- z-coordinate
    
    """
    cdef double[:] P_clst = np.empty(3), B = np.empty(3)
    cdef double[:] E0 = np.empty(3), E1 = np.empty(3), 
    cdef np.ndarray[double, ndim=2] E_clst = np.empty((3,4))
    cdef int inside, iclst
    cdef float a1, b1, c1, d1, e11, detinv, s1, t1
    cdef float a2, b2, c2, d2, e2, s2, t2
    cdef float dclst, f, dis1, dis2, dis3, dis4
    
    inside = 0
    # Top triangle (P1P2P4) of trapezoid:
    B[0] = P[0][0]
    B[1] = P[1][0]
    B[2] = P[2][0]
    E0[0] = P[0][1] - B[0]
    E0[1] = P[1][1] - B[1]
    E0[2] = P[2][1] - B[2]
    E1[0] = P[0][3] - B[0]
    E1[1] = P[1][3] - B[1]
    E1[2] = P[2][3] - B[2]
    a1 = E0[0]*E0[0] + E0[1]*E0[1] + E0[2]*E0[2]
    b1 = E0[0]*E1[0] + E0[1]*E1[1] + E0[2]*E1[2]
    c1 = E1[0]*E1[0] + E1[1]*E1[1] + E1[2]*E1[2]
    d1 = E0[0]*B[0] + E0[1]*B[1] + E0[2]*B[2]
    e11 = E1[0]*B[0] + E1[1]*B[1] + E1[2]*B[2]
    detinv = 1.0/(a1*c1-b1*b1)
    s1 = (b1*e11-c1*d1)*detinv
    t1 = (b1*d1-a1*e11)*detinv

    if (s1 >= 0) and (s1 <= 1) and (t1 >= 0) and (t1 <= 1) and (s1 + t1 <= 1):
        inside = 1
        P_clst[0] = B[0]+s1*E0[0]+t1*E1[0]
        P_clst[1] = B[1]+s1*E0[1]+t1*E1[1]
        P_clst[2] = B[2]+s1*E0[2]+t1*E1[2]
        dclst = P_clst[0]**2+P_clst[1]**2+P_clst[2]**2
        dclst = sqrt(dclst)
        return(dclst, P_clst)

    # Bottom triangle (P3P2P4) of trapezoid:
    B[0] = P[0][2]
    B[1] = P[1][2]
    B[2] = P[2][2]
    E0[0] = P[0][1] - B[0]
    E0[1] = P[1][1] - B[1]
    E0[2] = P[2][1] - B[2]
    E1[0] = P[0][3] - B[0]
    E1[1] = P[1][3] - B[1]
    E1[2] = P[2][3] - B[2]
    a2 = E0[0]*E0[0] + E0[1]*E0[1] + E0[2]*E0[2]
    b2 = E0[0]*E1[0] + E0[1]*E1[1] + E0[2]*E1[2]
    c2 = E1[0]*E1[0] + E1[1]*E1[1] + E1[2]*E1[2]
    d2 = E0[0]*B[0] + E0[1]*B[1] + E0[2]*B[2]
    e2 = E1[0]*B[0] + E1[1]*B[1] + E1[2]*B[2]
    detinv = 1.0/(a2*c2-b2*b2)
    s2 = (b2*e2-c2*d2)*detinv
    t2 = (b2*d2-a2*e2)*detinv

    if (s2 >= 0) and (s2 <= 1) and (t2 >= 0) and (t2 <= 1) and (s2 + t2 <= 1):
        inside = 2
        P_clst[0] = B[0]+s2*E0[0]+t2*E1[0]
        P_clst[1] = B[1]+s2*E0[1]+t2*E1[1]
        P_clst[2] = B[2]+s2*E0[2]+t2*E1[2]
        dclst = P_clst[0]**2+P_clst[1]**2+P_clst[2]**2
        dclst = sqrt(dclst)
        return(dclst, P_clst)

    # Four Edges of trapezoid:
    # ---- Edge 1: top edge (P1P2):
    if -d1 < 0:
        f = 0
    elif -d1 > a1:
        f = 1
    else:
        f = -d1/a1
    E_clst[0][0] = P[0][0] + f * (P[0][1]-P[0][0])
    E_clst[1][0] = P[1][0] + f * (P[1][1]-P[1][0])
    E_clst[2][0] = P[2][0] + f * (P[2][1]-P[2][0])

    # ---- Edge 2: side edge (P3P2):
    if -d2 < 0:
        f = 0
    elif -d2 > a2:
        f = 1
    else:
        f = -d2/a2
    E_clst[0][1] = P[0][2] + f * (P[0][1]-P[0][2])
    E_clst[1][1] = P[1][2] + f * (P[1][1]-P[1][2])
    E_clst[2][1] = P[2][2] + f * (P[2][1]-P[2][2])

    # ---- Edge 3: bottom edge (P3P4):
    if -e2 < 0:
        f = 0
    elif -e2 > c2:
        f = 1
    else:
        f = -e2/c2
    E_clst[0][2] = P[0][2] + f * (P[0][3]-P[0][2])
    E_clst[1][2] = P[1][2] + f * (P[1][3]-P[1][2])
    E_clst[2][2] = P[2][2] + f * (P[2][3]-P[2][2])    

    # ---- Edge 4: side edge (P1P4):
    if -e11 < 0:
        f = 0
    elif -e11 > c1:
        f = 1
    else:
        f = -e11/c1
    E_clst[0][3] = P[0][0] + f * (P[0][3]-P[0][0])
    E_clst[1][3] = P[1][0] + f * (P[1][3]-P[1][0])
    E_clst[2][3] = P[2][0] + f * (P[2][3]-P[2][0])  

    # Find closest point on the edges of Trapezoid:
    dis1 = E_clst[0][0]**2 + E_clst[1][0]**2 + E_clst[2][0]**2
    dis2 = E_clst[0][1]**2 + E_clst[1][1]**2 + E_clst[2][1]**2
    dis3 = E_clst[0][2]**2 + E_clst[1][2]**2 + E_clst[2][2]**2
    dis4 = E_clst[0][3]**2 + E_clst[1][3]**2 + E_clst[2][3]**2
    iclst = 0
    dclst = dis1
    if dis2 < dclst:
        iclst = 1
        dclst = dis2
    if dis3 < dclst:
        iclst = 2
        dclst = dis3
    if dis4 < dclst:
        iclst = 3
        dclst = dis4  
    dclst = sqrt(dclst)
    P_clst[0] = E_clst[0][iclst]
    P_clst[1] = E_clst[1][iclst]
    P_clst[2] = E_clst[2][iclst]

    return(dclst, P_clst)

@boundscheck(False)
@wraparound(False)
def Simulate_Rupture_Surface(eqn, eqType, region, em, elon, elat, hypd, Category,
                             strike=None, dip=None, rake=None, 
                             strike2=None, dip2=None, rake2=None,
                             mech=None, saveto=None):
    """
    
    Simulate earthquake rupture surface that minimizes difference from median distance of a 
    pseudo-grid of sites. Randomized fault attributes are dependent on the fault category.
    Fault scaling model depends on earthquake type.
    
    Required Input Arguments:
    - eqn = unique integer identifier for the event
    - eqType = type of eathquake:
               "crustal" = shallow-crustal earthquake
               "intraslab" = intraslab-subduction earthquake
               "interface" = interface-subduction earthquake
               "stable" = stable continental earthquake
    - region = geographic region where the earthquake occured:
               0 = Japan
               1 = Chile
               2 = Other
    - em = earthquake moment magnitude (Mw)
    - elon = hypocenter longitude (degrees)
    - elat = hypocenter latitude (degrees)
    - hypd = hypocenter depth (km)
    - Category = code for rupture simulation constraints:
                 "A" = strike, dip, and rake for first nodal plane solution preferred 
                       (optional "strike", "dip", and "rake" arguments are required)
                 "B" = strike, dip, and rake for second nodal plane solution preferred
                       (optional "strike2", "dip2", and "rake2" arguments are required)
                 "C" = strike, dip, and rake are known for two nodal planes, and neither 
                       is preferred
                       (optional "strike", "dip", and "rake" arguments are required)
                       (optional "strike2", "dip2", and "rake2" arguments are required)
                 "D" = One nodal plane solution for strike, dip, and rake; randomize the 
                       strike and dip
                       (optional "strike", "dip", and "rake" arguments are required)
                 "E" = No nodal plane solutions; randomize strike, dip, and rake
                       (dip and rake are assigned based on faulting mechanism)
                       (if optional "mech" argument is not speficied, simulations randomly
                       assign one)
                       
    Optional Input Arguments:
    - strike = strike (degrees) of first nodal plane solution (default None)
               (required if Category "A", "C", or "D" is specified) 
    - dip = dip (degrees) of first nodal plane solution (default None)
            (required if Category "A", "C", or "D" is specified) 
    - rake = rake (degrees) of first nodal plane solution (default None)
             (required if Category "A", "C", or "D" is specified) 
    - strike2 = strike (degrees) of first nodal plane solution (default None)
                (required if Category "B" or "C" is specified) 
    - dip2 = dip (degrees) of first nodal plane solution (default None)
             (required if Category "B" or "C" is specified) 
    - rake2 = rake (degrees) of first nodal plane solution (default None)
              (required if Category "B" or "C" is specified) 
    - mech = faulting mechanism (default None)
             (if not specified, and eqType is "crustal" or "stable", simulations randomply
             assign one)
             "SS" = strike-slip
             "RV" = reverse-thrust
             "NM" = normal-thrust
    - saveto = directory where you would like to save output files (default None)
               (if not specified, results are not saved as files on local memory)
               
    Returns:
    - SIM_RESULTS = pandas DataFrame containing all simulated rupture surfaces and statistics
    - SELECTED_FAULT = pandas DataFrame containing the selected rupture surface and statistics
    
    """
    # Initialize DataTypes -----------------------------------------------------------------
    cdef int idum = -eqn
    cdef int mxd, nxf, nsim, nmed, d, ia, nps, i_s, ksim, nffpts
    cdef float seisd, sumold, deg2rad, rad2deg, pi, pio2, ang, a, b, tp
    cdef float snm, srv, xf, yf, xb, yb, rwh, rwv, hh, SUM
    cdef float area_ave, area_sig, ar_ave, ar_sig, rl_ave, rl_sig, rw_ave, rw_sig
    cdef float topd_mx, topd_mn, botd_mx, botd_mn, fflon, fflat
    cdef double[:] xdf, ssd, nmd, rvd, ied0, ied1, ied2, iad, hypx, hypxe, hypxa
    cdef double[:] psx, psy, area, fd, fl, ar, frake, fstrike, fdip, rw, rl, topd, botd
    cdef double[:] ckseisd, clsd, jbd, cld, ffx, ffy, ffz, p_clst
    cdef double[:] medrd, medjbd, msim, jbsim, sigrd, averd, fflon2, fflat2
    cdef np.ndarray[double, ndim=2] pt, dminr, dminjbr
    cdef np.ndarray[double, ndim=3] rpx, rpy, rpz
    #cdef np.ndarray[float, ndim=2] E_clst = np.empty((3,4))

    # Set random seed based on eqn ---------------------------------------------------------
    random.seed(idum)
    
    # Define Constants ---------------------------------------------------------------------
    mxd = 1000; nxf = 11; nsim = 101; nmed = 51; seisd=3.0; sumold = 1.0e20
    xdf = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ssd = np.array([0, 0.025, 0.05, 0.1, 0.175, 0.275, 0.4, 0.55, 0.7, 0.85, 1])
    nmd = np.array([0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.25, 0.5, 1])
    rvd = np.array([0, 0.03, 0.07, 0.12, 0.2, 0.3, 0.42, 0.55, 0.69, 0.84, 1])
    ied0 = np.array([0.000, 0.024, 0.085, 0.206, 0.389, 0.599, 0.783, 0.906, 0.969, 0.993, 1.00])
    ied1 = np.array([0.0, 0.002, 0.012, 0.044, 0.121, 0.262, 0.460, 0.671, 0.843, 0.950, 1.00])
    ied2 = np.array([0.00, 0.013, 0.053, 0.143, 0.297, 0.500, 0.703, 0.857, 0.947, 0.987, 1.00])
    iad = np.array([0.00, 0.012, 0.051, 0.139, 0.294, 0.500, 0.706, 0.861, 0.949, 0.988, 1.00]) 
    hypx = np.array([0, 0.05, 0.125, 0.225, 0.35, 0.5, 0.65, 0.775, 0.875, 0.95, 1])
    hypxe = np.array([0.000, 0.007, 0.034, 0.112, 0.272, 0.500, 0.728, 0.888, 0.966, 0.993, 1.00])
    hypxa = np.array([0.00, 0.015, 0.057, 0.148, 0.301, 0.500, 0.699, 0.852, 0.943, 0.985, 1.00])
    deg2rad=asin(1.0)/90; rad2deg=90./asin(1.0); 
    pi=2.0*asin(1.0); pio2=asin(1.0)
    
    # Create a pseudo-grid -----------------------------------------------------------------
    psx = np.empty(mxd); psy = np.empty(mxd); 
    psx[0] = 0.0; psy[0] = 0.0; nps = 0
    for d in np.arange(2,20,2):
        for ia in range(24):
            ang = pi * 2.0 * ia / 24
            nps += 1
            psx[nps] = d * sin(ang)
            psy[nps] = d * cos(ang)     
    for d in np.arange(25,55,5):
        for ia in range(24):
            ang = pi * 2.0 * ia / 24
            nps += 1
            psx[nps] = d * sin(ang)
            psy[nps] = d * cos(ang)      
    for d in np.arange(60,110,10):
        for ia in range(24):
            ang = pi * 2.0 * ia / 24
            nps += 1
            psx[nps] = d * sin(ang)
            psy[nps] = d * cos(ang)    
    for d in np.arange(125,325,25):
        for ia in range(24):
            ang = pi * 2.0 * ia / 24
            nps += 1
            psx[nps] = d * sin(ang)
            psy[nps] = d * cos(ang)     
    psx = psx[:nps]
    psy = psy[:nps]
    
    # Initilize variable arrays and pandas DataFrames --------------------------------------
    area = np.empty(nsim); 
    fd = np.empty(nsim); 
    fl = np.empty(nsim); 
    ar = np.empty(nsim); 
    
    frake = np.empty(nsim); 
    fstrike = np.empty(nsim); 
    fdip = np.empty(nsim); 
    rw = np.empty(nsim); 
    
    rl = np.empty(nsim); 
    topd = np.empty(nsim); 
    botd = np.empty(nsim); 
    
    rpx = np.empty((nsim,4,2)); 
    rpy = np.empty((nsim,4,2)); 
    rpz = np.empty((nsim,4,2)); 
    
    ckseisd = np.empty(nsim); 
    clsd = np.empty(nsim); 
    jbd = np.empty(nsim);  
    cld = np.empty(nsim);
    
    ffx = np.empty(5); 
    ffy = np.empty(5); 
    ffz = np.empty(5); 
    
    #fflon = np.empty(5); 
    #fflat = np.empty(5); 
    p_clst = np.empty(3)
    
    dminr = np.empty((nps,nsim)); 
    dminjbr = np.empty((nps,nsim)); 
    pt = np.empty((3,4)); 
    
    medrd = np.empty(nps); 
    medjbd = np.empty(nps); 
    msim = np.empty(nps); 
    jbsim = np.empty(nps)
    
    sigrd = np.zeros(nps); 
    averd = np.zeros(nps) # keep as zeros
    SIM_RESULTS = pd.DataFrame(columns=['SimN','FL','FD','PX1','PY1','PZ1','PX2','PY2','PZ2',
                                        'PX3','PY3','PZ3','PX4','PY4','PZ4','Strk','Dip','Rake',
                                        'Mag','AR','Area','RL','RW','TopD','BotD'])
    SELECTED_FAULT = pd.DataFrame(columns = ['EQID','EQNAME','Ksim','ELat','ELon','Hypd',
                                             'FL','FD','PX1','PY1','PZ1','PX2','PY2','PZ2',
                                             'PX3','PY3','PZ3','PX4','PY4','PZ4',
                                             'Strike','Dip','Rake','Mag','AspR','AveAR',
                                             'SigAR','Area','AveArea','SigArea','RL','AveRL',
                                             'SigRL','RW','AveRW','SigRW'])
    
    # Simulate the fault rupture surface nsim times ----------------------------------------
    for isim in range(nsim):
        
        # Simulate fault area (area) .......................................................
        if eqType == 'crustal':
            # Wells & Coppersmith (1994):
            a = -3.49; b = 0.91
            area[isim] = 10.0**(a + b*em + 0.24*gasdev(idum))
            # Leonard (2010):
            
        elif eqType == 'interface':
            a = -3.8290; b = 1.0
            area[isim] = 10.0**(a + b*em + 0.270*gasdev(idum)) 
        elif eqType == 'intraslab':
            a = -3.251; b = 0.890
            area[isim] = 10.0**(a + b*em + 0.184*gasdev(idum))
        elif eqType == 'stable':
            a = 1.5; b = 6.38
            area[isim] = 10.0**(((1.5*em+16.05)/10**(7)-b)/a+0.8*gasdev(idum))
        else:
            print('ERROR: eqType "' + eqType + '" not found.')

        # Simulate fault mechanism - if unknown (mech) .....................................
        if mech == None:
            tp = random.random()
            if tp <= 1/3:
                mech = 'SS'
            elif tp <= 2/3:
                mech = 'RV'
            else:
                mech = 'NM'

        # Simulate hypocenter's down-dip position on the fault (fd) ........................
        tp = random.random()
        if (eqType == 'crustal') or (eqType == 'stable'):
            if mech == 'SS':
                snm = 0.0; srv = 0.0
                fd[isim] = discrete(nxf,xdf,ssd,tp)
                if rake == None:
                    rake = 0.0
                if dip == None:
                    dip = 89.99
            elif mech == 'RV':
                snm = 0.0; srv = 1.0
                fd[isim] = discrete(nxf,xdf,ssd,tp)
                if rake == None:
                    rake = 90.0
                if dip == None:
                    dip = 40.0
            elif mech == 'NM':
                snm = 1.0; srv = 0.0
                fd[isim] = discrete(nxf,xdf,ssd,tp)
                if rake == None:
                    rake = -90.0
                if dip == None:
                    dip = 55.0
            else:
                print('ERROR: mech "' + mech + '" not found.')

        # ----- simulate intraslab hypocenter along dip. snm and srv are NOT used for intraslab event, dummies:
        elif eqType == 'intraslab':
            snm = 0.0; srv = 0.0
            fd[isim] = discrete(nxf,xdf,iad,tp)
        # ----- simulate interface hypocenter along dip. snm and srv are NOT used for interface event, dummies:
        elif eqType == 'interface':
            snm = 0.0; srv = 0.0
            if region == 0: # Japan
                fd[isim] = discrete(nxf,xdf,ied0,tp)
            elif region == 1: # Chile
                fd[isim] = discrete(nxf,xdf,ied1,tp)
            elif region == 2: # Other
                fd[isim] = discrete(nxf,xdf,ied2,tp)
        
        # Simulate aspect ratio (ar) .......................................................
        if eqType == 'crustal':
            if em > 4.0:
                ar[isim] = 10.0**( (0.017518 - 0.004717*snm - 0.010992*srv) * (em - 4.0)**3.096977 + 0.16*gasdev(idum) )
            else:
                ar[isim] = 10.0**( 0.0 + 0.16*gasdev(idum) )
        elif eqType == 'interface':
            if em > 7.25:
                ar[isim] = 10.0**( 0.2759*(em - 7.25) + 0.192*gasdev(idum) )    # std = 0.441 in ln unit
            else:
                ar[isim] = 10.0**( 0.0 + 0.0717*gasdev(idum) )                  # std = 0.165 in ln unit
        elif eqType == 'intraslab':
            if em > 6.5:
                ar[isim] = 10.0**( 0.0938*(em - 6.5) + 0.164*gasdev(idum) )     # std = 0.378 in ln unit
            else:
                ar[isim] = 10.0**( 0.0 + 0.104*gasdev(idum) )                   # std = 0.239 in ln unit
        elif eqType == 'stable':
            ar[isim] = 1*gasdev(idum)
            
        
        # Simulate fault strike and dip ....................................................
        frake[isim] = rake

        tp = random.random()
        if Category == 'E':
            fstrike[isim] = tp*360.0
            fdip[isim] = dip
        elif Category == 'D':
            fstrike[isim] = (tp - 0.5)*2.0*30.0 + strike     # +/- 30 degrees uncertainty
            tp = random.random()
            fdip[isim] = (tp - 0.5)*2.0*10.0 + dip           # +/- 10 degrees uncertainty
            if fdip[isim] < min(10, dip):
                fdip[isim] = min(10, dip)                    # floor on dip for category D
        elif Category == 'C':
            if tp < 0.5:     # Use first nodal plane
                fstrike[isim] = strike
                fdip[isim] = dip
                frake[isim] = rake
            else:            # Use second nodal plane
                fstrike[isim] = strike2
                fdip[isim] = dip2
                frake[isim] = rake2
        elif Category == 'B':
            fstrike[isim] = strike2
            fdip[isim] = dip2
            frake[isim] = rake2
        elif Category == 'A':
            fstrike[isim] = strike
            fdip[isim] = dip
            frake[isim] = rake
        else:
            print('ERROR: Category "' + Category + '" not found.')

        fstrike[isim] = fstrike[isim] * deg2rad

        if fdip[isim] >= 90:
            fdip[isim] = 89.99 # required to avoid asymptotic solution (i.e., cannot divide by zero)
        fdip[isim] = fdip[isim] * deg2rad
        
        # Simulate hypocenter's position along-strike (fl) .................................
        tp = random.random()
        if eqType == 'crustal':
            fl[isim] = discrete(nxf,xdf,hypx,tp)
        elif eqType == 'intraslab':
            fl[isim] = discrete(nxf,xdf,hypxa,tp)
        elif eqType == 'interface':
            fl[isim] = discrete(nxf,xdf,hypxe,tp)
        
        # Compute other parameters .........................................................
        rw[isim] = sqrt(area[isim]/ar[isim])                   # rupture width
        rl[isim] = area[isim]/rw[isim]                         # rupture length
        xf = sin(fstrike[isim])*rl[isim]*fl[isim]              # x-coord of epicenter
        yf = cos(fstrike[isim])*rl[isim]*fl[isim]              # y-coord of epicenter 
        xb = sin(fstrike[isim]+pi)*rl[isim]*(1.0-fl[isim])
        yb = cos(fstrike[isim]+pi)*rl[isim]*(1.0-fl[isim])
        rwh = rw[isim]*cos(fdip[isim])
        rwv = rw[isim]*sin(fdip[isim])
        topd[isim] =hypd-rwv*fd[isim] 
        if topd[isim] < 0:      # Limit topd to 0.0
            topd[isim] = 0.0
            fd[isim] = hypd/rwv
        botd[isim] = topd[isim] + rwv
        
        # Compute top points ...............................................................
        rpx[isim][0][0] = xf + sin(fstrike[isim] - pio2)*rwh*fd[isim]
        rpx[isim][0][1] = xb + sin(fstrike[isim] - pio2)*rwh*fd[isim]
        rpy[isim][0][0] = yf + cos(fstrike[isim] - pio2)*rwh*fd[isim]
        rpy[isim][0][1] = yb + cos(fstrike[isim] - pio2)*rwh*fd[isim]
        rpz[isim][0][0] = topd[isim]
        rpz[isim][0][1] = topd[isim]
        
        # Compute bottom points ............................................................
        rpx[isim][1][0] = xf + sin(fstrike[isim] + pio2)*rwh*(1.0 - fd[isim])
        rpx[isim][1][1] = xb + sin(fstrike[isim] + pio2)*rwh*(1.0 - fd[isim])
        rpy[isim][1][0] = yf + cos(fstrike[isim] + pio2)*rwh*(1.0 - fd[isim])
        rpy[isim][1][1] = yb + cos(fstrike[isim] + pio2)*rwh*(1.0 - fd[isim])
        rpz[isim][1][0] = botd[isim]
        rpz[isim][1][1] = botd[isim]
        
        # Locate top of seisd ..............................................................
        if topd[isim] < seisd:
            ckseisd[isim] = True
            hh = tan(pio2 - fdip[isim])*(hypd - seisd)
            rpx[isim][2][0] = xf + sin(fstrike[isim] - pio2)*hh
            rpx[isim][2][1] = xb + sin(fstrike[isim] - pio2)*hh
            rpy[isim][2][0] = yf + cos(fstrike[isim] - pio2)*hh
            rpy[isim][2][1] = yb + cos(fstrike[isim] - pio2)*hh
            rpz[isim][2][0] = seisd
            rpz[isim][2][1] = seisd
        else:
            ckseisd[isim] = False
        
        # Place line through hypocenter ....................................................
        rpx[isim][3][0] = xf
        rpx[isim][3][1] = xb
        rpy[isim][3][0] = yf
        rpy[isim][3][1] = yb
        rpz[isim][3][0] = hypd
        rpz[isim][3][1] = hypd
        
        # Compute closest distance for each pseudo-site ....................................
        for i_s in range(nps):

            # load data into pt array for rupture distances:
            pt[0][0] = rpx[isim][0][1] - psx[i_s]
            pt[1][0] = rpy[isim][0][1] - psy[i_s]
            pt[2][0] = rpz[isim][0][1] 

            pt[0][1] = rpx[isim][0][0] - psx[i_s]
            pt[1][1] = rpy[isim][0][0] - psy[i_s]
            pt[2][1] = rpz[isim][0][0] 

            pt[0][2] = rpx[isim][1][0] - psx[i_s]
            pt[1][2] = rpy[isim][1][0] - psy[i_s]
            pt[2][2] = rpz[isim][1][0] 

            pt[0][3] = rpx[isim][1][1] - psx[i_s]
            pt[1][3] = rpy[isim][1][1] - psy[i_s]
            pt[2][3] = rpz[isim][1][1] 

            dminr[i_s][isim], p_clst_r = DEMDT(pt)

            # load data into pt array for JB distance:
            pt[2][0] = 0.0
            pt[2][1] = 0.0
            pt[2][2] = 0.0
            pt[2][3] = 0.0

            dminjbr[i_s][isim], p_clst_jb = DEMDT(pt)
        
        # Save Data to a pandas dataframe to be returned ...................................
        SIM_RESULTS.loc[len(SIM_RESULTS.index)+1,:] = [isim+1, fl[isim], fd[isim],
                                                       pt[0][0]+psx[i_s], pt[1][0]+psy[i_s], rpz[isim][0][1],
                                                       pt[0][1]+psx[i_s], pt[1][1]+psy[i_s], rpz[isim][0][0],
                                                       pt[0][2]+psx[i_s], pt[1][2]+psy[i_s], rpz[isim][1][0],
                                                       pt[0][3]+psx[i_s], pt[1][3]+psy[i_s], rpz[isim][1][1],
                                                       fstrike[isim]*rad2deg, fdip[isim]*rad2deg, frake[isim],
                                                       em, ar[isim], area[isim], 
                                                       rl[isim], rw[isim], topd[isim], botd[isim] ]
    
    # Locate median value of dminr for each site -------------------------------------------
    for i_s in range(nps):
        for isim in range(nsim):
            cld[isim] = dminr[i_s][isim]
            jbd[isim] = dminjbr[i_s][isim]
        medrd[i_s] = np.median(cld)
        medjbd[i_s] = np.median(jbd)
        for isim in range(nsim):
            if medrd[i_s] == dminr[i_s][isim]:
                msim[i_s] = isim
            if medjbd[i_s] == dminjbr[i_s][isim]:
                jbsim[i_s] = isim
    
    # Find simulation that minimizes difference between med rupd and rupd across all sites -
    for isim in range(nsim):
        SUM = 0
        for i_s in range(nps):
            SUM += ( medrd[i_s] - dminr[i_s][isim] )**2
        if SUM < sumold:
            ksim = isim
            sumold = SUM
    
    # Compute log mean value and log sigma of dminr for each site --------------------------
    for i_s in range(nps):
        for isim in range(nsim):
            cld[isim] = dminr[i_s][isim]
            averd[i_s] += log(cld[isim])
            sigrd[i_s] += log(cld[isim])**2
        averd[i_s] = averd[i_s]/nsim
        sigrd[i_s] = ( (sigrd[i_s]/nsim - averd[i_s]**2) * nsim/(nsim-1) )
        if sigrd[i_s] < 0:
            sigrd[i_s] = 0
        else:
            sigrd[i_s] = sqrt(sigrd[i_s])
        averd[i_s] = exp(averd[i_s])
    
    # Compute the statistics of simulated faults -------------------------------------------
    area_ave, area_sig, ar_ave, ar_sig, rl_ave, rl_sig, rw_ave, rw_sig = 0,0,0,0,0,0,0,0
    for isim in range(nsim):
        area_ave += log10(area[isim])
        area_sig += log10(area[isim])**2
        ar_ave += log10(ar[isim])
        ar_sig += log10(ar[isim])**2
        rl_ave += log10(rl[isim])
        rl_sig += log10(rl[isim])**2
        rw_ave += log10(rw[isim])
        rw_sig += log10(rw[isim])**2
    area_ave = area_ave/nsim
    area_sig = sqrt( (area_sig/nsim - area_ave**2) * nsim/(nsim-1) )
    area_ave = 10**area_ave
    ar_ave = ar_ave/nsim
    ar_sig = sqrt( (ar_sig/nsim - ar_ave**2) * nsim/(nsim-1) )
    ar_ave = 10**ar_ave
    rl_ave = rl_ave/nsim
    rl_sig = sqrt( (rl_sig/nsim - rl_ave**2) * nsim/(nsim-1) )
    rl_ave = 10**rl_ave
    rw_ave = rw_ave/nsim
    rw_sig = sqrt( (rw_sig/nsim - rw_ave**2) * nsim/(nsim-1) )
    rw_ave = 10**rw_ave

    topd_mx = max(topd)
    topd_mn = min(topd)
    botd_mx = max(botd)
    botd_mn = min(botd)
   
    # Convert selected FF xy to LL using haversine formulation -----------------------------
    nffpts=5

    ffx[0]=rpx[ksim][0][1]
    ffy[0]=rpy[ksim][0][1]
    ffz[0]=rpz[ksim][0][1]

    ffx[1]=rpx[ksim][0][0]
    ffy[1]=rpy[ksim][0][0]
    ffz[1]=rpz[ksim][0][0]

    ffx[2]=rpx[ksim][1][0]
    ffy[2]=rpy[ksim][1][0]
    ffz[2]=rpz[ksim][1][0]

    ffx[3]=rpx[ksim][1][1]
    ffy[3]=rpy[ksim][1][1]
    ffz[3]=rpz[ksim][1][1]

    ffx[4]=rpx[ksim][0][1]
    ffy[4]=rpy[ksim][0][1]
    ffz[4]=rpz[ksim][0][1]

    fflon2, fflat2 = XY2LL(elon, elat, ffx, ffy)

    # Save Data to a pandas dataframe to be returned ---------------------------------------
    SELECTED_FAULT.loc[0,:] = [eqn, eqn, ksim, elat, elon, hypd, fl[ksim], fd[ksim],
                               fflon2[0], fflat2[0], ffz[0],
                               fflon2[1], fflat2[1], ffz[1],
                               fflon2[2], fflat2[2], ffz[2],
                               fflon2[3], fflat2[3], ffz[3],
                               fstrike[ksim]*rad2deg, fdip[ksim]*rad2deg, frake[ksim], em, 
                               ar[ksim], ar_ave, ar_sig, 
                               area[ksim], area_ave, area_sig,
                               rl[ksim], rl_ave, rl_sig, 
                               rw[ksim], rw_ave, rw_sig ] 
    
    return(SIM_RESULTS, SELECTED_FAULT)