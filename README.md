# ccldpy

Python package for simulating earthquake rupture surface representation.

The CCLD program was originially coded in Fortran by Brian Chiou and Robert Youngs (Chiou and Youngs 2008; Appendix B), then later updated by Brian Chiou during the NGA-West2 and NGA-Subduction projects (Chiou et al. 2008; Contreras et al. 2020). This version of CCLD replicates the originial simulation code in Python. Key changes to the program include the following:
1. Updated magnitude-area scaling relationship for shallow-<em>crustal</em> type events (Leonard 2014)
2. Added magnitude-area scaling relationship for <em>stable</em>-continental type events (Leonard 2014). Shallow-crustal type relations for the position of hypocenter along the rupture surface are assumed.
3. Flexibility to select different magnitude-based scaling relationships.

The current version of ccldpy (0.0.1) does not compute distances for real seismic stations, which can be perfomed using the P4CF program (Chiou, B.S-J. 2021).

# Installation/Usage

Cython code is pre-compiled into the <em>ccldpy.cp39-win_amd64.pyd</em> file, which can be imported into Python environments. 
### Simplified instructions:

1. Download <em>ccldpy.cp39-win_amd64.pyd</em> into your {<em>desired path</em>}.
2. Import <em>ccldpy</em> into your code:<br>
```python
import sys
sys.path.append({desired path})
import ccldpy
```
# 
### Simulation Methods:

The current version of ccldpy (0.0.1) supports five methods of simulation, which are specificed using the <em>Category</em> indicators described below. For each method, the program runs 101 simulations, and identifies the rupture surface of the simulation that minimizes the median distance computed for a pseudo-grid of locations.

A = One or two nodal plane solutions (strike, dip, and rake) are known, however the first solution is preferred. Only the area, aspect-ratio, and position of hypocenter on the rupture surface are randomized between simulations assuming a Gaussian Distribution centered about the median values.

B = Two nodal plane solutions are known, however the second solution is preferred. Only the area, aspect-ratio, and position of hypocenter on the rupture surface are randomized between simulations.

C = Two nodal plane solutions are known, and neither is preferred over the other. Each simulation will randomly select which nodal plane solution to use, and the area, aspect-ratio, and position of hypocenter on the rupture surface are also randomized.

D = Only one nodal plane solution is known or assumed, but with some uncertainty. The given rake angle is used to assign rupture mechanism (<u>S</u>trike-<u>S</u>lip, <u>N</u>or<u>M</u>al dip-slip, or <u>R</u>e<u>V</u>erse dip-slip - for shallow-crustal or stable-continetal type events) and the strike and dip are randomized with &#177; 30<sup>o</sup> and &#177; 10<sup>o</sup>, respectively. The area, aspect-ratio, and position of hypocenter on the rupture surface are also randomized between simulations.

E = No nodal plane solutions are known or assumed. Rake (faulting mechanism), strike [0<sup>o</sup> - 360<sup>o</sup>), and dip [0<sup>o</sup> - 90<sup>o</sup>] are randomly assigned with equal probability during each simulation.

![Category Illustration](https://user-images.githubusercontent.com/71461454/220185818-708986c3-28ff-4dfa-b54b-e225dfe261f3.png)
<b>Figure 1:</b> Schematic illustration of simulation results for each <em>Category</em>.


# Main Function
```python
Simulate_Rupture_Surface(eqn, eqType, region, em, elon, elat, hypd, Category,
                         strike=None, dip=None, rake=None, 
                         strike2=None, dip2=None, rake2=None,
                         mech=None, saveto=None, model=None)
```

### Input Parameters:

#### Required Keys:

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
                 "C" = strike, dip, and rake are known for two nodal planes, and neither is preferred
                       (optional "strike", "dip", and "rake" arguments are required)
                       (optional "strike2", "dip2", and "rake2" arguments are required)
                 "D" = One nodal plane solution for strike, dip, and rake; randomize the strike and dip
                       (optional "strike", "dip", and "rake" arguments are required)
                 "E" = No nodal plane solutions; randomize strike, dip, and rake
                       (dip and rake are assigned based on faulting mechanism)
                       (if optional "mech" argument is not speficied, simulations randomly assign one)
                            
#### Optional Keys:
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
             (if not specified, and eqType is "crustal" or "stable", simulations randomply assign one)
             "SS" = strike-slip
             "RV" = reverse-thrust
             "NM" = normal-thrust
    - saveto = directory where you would like to save output files (default None)
               (if not specified, results are not saved as files on local memory)
    - model = Magnitude-based scaling relationship model to simulate area and aspect ratio (width/length)
             "wells_and_coppersmith_1994" = Wells, D.L. and K.J. Coppersmith (1994) New empirical relationships 
                              among magnitude, rupture length, rupture width, rupture area, and surface 
                              displacement, Bull. Seismol. Soc. Am. 84 (4), 974 - 1002.
                              Can be used for "crustal" and "stable" eqTypes
             "leonard_2014" = Leonard, M. (2014) Self-consistent earthquake fault-scaling relations: Update
                              and extension to stable continental strike-slip faults, Bull. Seismol. Soc. 
                              Am. 104 (6), 2953 - 2965.
                              Can be used for "custal" and "stable" eqTypes (default)
             "contreras_et_al_2022" = Contreras, V., J.P. Stewart, T. Kishida, R.B. Darragh, B.-S.J. Chiou,
                              S. Mazzoni, R.R. Youngs, N.M. Kuehn, S.K. Ahdi, K. Wooddell, R. Boroschek, F.
                              Rojas, and J. Ordenes (2022) NGA-Sub source and path database, Earthq. Spectra
                              38 (2), 799 - 840.
                              Can be used for "intraslab" and "interface" eqTypes (default)
   
#### Returns:
    - SIM_RESULTS = pandas DataFrame containing all simulated rupture surfaces and statistics
    - SELECTED_FAULT = pandas DataFrame containing the selected rupture surface and statistics

# References:

Chiou, B. S.-J. (2021). P4CF [https://github.com/bc88bc/P4CF]

Chiou ,B. S.‐J., and Youngs R. R. (2008). <em>NGA Model for Average Horizontal Component of Peak Ground Motion and Response Spectra</em>, PEER Rept. 2008/09, Pacific Earthquake Engineering Research Center, Berkeley, California.

Contreras V. Stewart J. P. Kishida T. Darragh R. B. Chiou B. S.‐J. Mazzoni S. Kuehn N. Ahdi S. K. Wooddell K., and Youngs R. R., et al. (2020). Source and path database, in <em>Data Resources for NGA‐Subduction Project</em>, Stewart J. P. (Editor), Chapter 4, PEER Rept. 2020/02, Pacific Earthquake Engineering Research Center, UC Berkeley, Berkeley, California.

Contreras, V., J.P. Stewart, T. Kishida, R.B. Darragh, B.-S.J. Chiou,S. Mazzoni, R.R. Youngs, N.M. Kuehn, S.K. Ahdi, K. Wooddell, R. Boroschek, F. Rojas, and J. Ordenes (2022) NGA-Sub source and path database, <em>Earthq. Spectra</em> <b>38</b> (2), 799 - 840.

Leonard, M. (2014). Self-consistent earthquake fault-scaling relations: Update and extension to stable continental strike-slip faults. <em>Bull. Seismol. Soc. Am.</em> <b>104</b> (6), 2953 - 2965.

Veness, C. (n.d.). <em>Movable type scripts</em>. Calculate distance and bearing between two Latitude/Longitude

Wells, D.L. and K.J. Coppersmith (1994) New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement, <em>Bull. Seismol. Soc. Am.</em> <b>84</b> (4), 974 - 1002.
