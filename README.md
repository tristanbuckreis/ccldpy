# ccldpy

Python package for simulating earthquake rupture surface representation.


## Function 
```python
Simulate_Rupture_Surface(eqn, eqType, region, em, elon, elat, hypd, Category,
                         strike=None, dip=None, rake=None, 
                         strike2=None, dip2=None, rake2=None,
                         mech=None, saveto=None)
```

### Input Parameters:

#### Required Keys:

    eqn = unique integer identifier for the event
    eqType = type of eathquake:
             "crustal" = shallow-crustal earthquake
             "intraslab" = intraslab-subduction earthquake
             "interface" = interface-subduction earthquake
             "stable" = stable continental earthquake
    region = geographic region where the earthquake occured:
             0 = Japan
             1 = Chile
             2 = Other
    em = earthquake moment magnitude
    elon = hypocenter longitude (degrees)
    elat = hypocenter latitude (degrees)
    hypd = hypocenter depth (km)
    Category = code for rutpure simulation constraints:
               "A" = strike, dip, and rake are known for two nodal planes, however nodal plane 1 is preferrd.
               "B" = strike, dip, and rake are known for two nodal planes, however nodal plane 2 is preferrd.
               "C" = strike, dip, and rake are known for two nodal planes, and neither is preferred.
               "D" = One nodal plane solution for strike, dip, and rake. Randomize the strike and dip.
               "E" = No nodal plane solutions. Randomize strike. 
                            
#### Optional Keys:
    strike = strike (degrees) of nodal plane solution 1. Required if dip and/or rake are speficied. (default None)
    dip = dip (degrees) of nodal plane solution 1. Required if strike and/or rake are speficied. (default None)
    rake = rake (degrees) of nodal plane solution 1, Required if strike and/or dip are speficied. (default None)
    strike2 = strike (degrees) of nodal plane solution 2. Required if dip2 and/or rake2 are speficied. (default None)
    dip2 = dip (degrees) of nodal plane solution 2. Required if strike2 and/or rake2 are speficied. (default None)
    rake2 = rake (degrees) of nodal plane solution 2, Required if strike2 and/or dip2 are speficied. (default None)
    mech = faulting mechanism: (default None)
           "SS" = strike slip
           "RV" = reverse thrust
           "NM" = normal thrust
    saveto = directory where you would like to save output files. 
             (default None; the function will not save files, results are only returned as pandas DataFrame objects)
   
#### Returns:
    SIM_RESULTS = pandas DataFrame object with 101 simulated rupture surfaces.
    SELECTED_FAULT = pandas DataFrame object with selected representative rupture surface.
    
### Subroutines:

