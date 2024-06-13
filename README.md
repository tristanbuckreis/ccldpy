# ccldpy

Python package for simulating earthquake rupture surface representation for the purpose of computing finite-fault distance metrics (e.g., closest distance to the rupture surface, $R_{RUP}$; Joyner-Boore distance, $R_{JB}$).

The CCLD program was originally coded in Fortran by Brian Chiou and Robert Youngs (Chiou and Youngs 2008; Appendix B), then later updated by Brian Chiou during the NGA-West2 and NGA-Subduction projects (Chiou and Youngs 2008; Contreras et al. 2022). This version of CCLD replicates the original simulation code in Python. Key changes to the program include the following:
1. Added magnitude-scaling relationships for shallow-<em>crustal</em> type events (Wells & Coppersmith 1994; Leonard 2014; and Thingbaijam et al. 2017)
2. Added magnitude-scaling relationship for <em>stable</em>-continental type events (Leonard 2014). Shallow-crustal type relations for the position of hypocenter along the rupture surface are assumed (Chiou and Youngs 2008).
3. Added magnitude-scaling relationship for <em>interface</em> subduction type events (Thingbaijam et al. 2017). Subduction-type relations for the position of the hypocenter on the rupture surface are assumed (Contreras et al. 2022)
4. Flexibility to specify the number of simulations to use for each magnitude-scaling relation (e.g., use <em>N</em> simulations from one scaling-relation or partition <em>N</em> simulations into <em>M</em> appropriate scaling-relations.
   
The current version of ccldpy (2.0.0) does not compute distances for real seismic stations, which can be perfomed using the P4CF program (Chiou, B.S-J. 2021).

# Magnitude-Scaling Relationships:

Scaling-relationships are defined for active shallow-<em>crustal</em>, <em>stable</em>-continental, <em>interface</em>-subduction, and <em>intra-slab</em> subduction-type earthquakes. These relationships can only be used in their respective tectonic regimes, as summarized herein. 

### Active Shallow-Crustal Earthquakes:

Four scaling-relationships for active shallow-crustal type earthquakes are implemented within <em>ccldpy</em>:
1. Wells and Coppersmith (1994): relationships for area, length, and width; one set of coefficients for all styles-of-faulting (i.e., mechanism)
2. Chiou and Youngs (2008): only an aspect ratio relationship; separate coefficients for strike-slip, normal, and reverse styles-of-faulting
3. Lenoard (2014): relationships for area, length, and width; separate coefficients for strike-slip and dip-slip (normal and reverse) styles-of-faulting
4. Thingbaijam et al. (2017): relationships for area, length, and width; separate coefficients for strike-slip, normal, and reverse styles of faulting

Wells and Coppersmith (1994) present separate coefficients for different styles-of-faulting, however the authors do not recommend using these coefficients and instead recommend the "all" set of coefficients. Chiou and Youngs (2008), Leonard (2014), and Thingbaijam et al. (2017) recommend using different sets of coefficients for different styles-of-faulting. Since Chiou and Youngs (2008) only present an aspect ratio ($AR$) relationship, it must be used in combination with the other models to develop a rupture surface (details discussed in the <em>Generating Rupture Geometries</em> section) of the documentation.

![ccldpy - crustal scaling relations](https://github.com/tristanbuckreis/ccldpy/assets/71461454/bf358c1c-c250-45a1-b03d-6dc3f2e4e162)
<b>Figure 1:</b> Median magnitude-scaling relationships for active shallow-crustal events.

All of these models were develop using data from events with magnitudes generally greater than 6.0, meaning the relationships must be extrapolated for use with smaller-magnitude events. The $L$, $W$, and $AR$ (= $L/W$) relationships of Wells and Coppersmith (1994), Leonard (2014), and Thingbaijam et al. (2017) for small-magnitudes (generally < 5.5 - 6.0) are not reflective of the as-published relations. This is because extrapolation of these relationships to small-magnitudes leads to unrealistic ruptures (e.g., $AR$ < 0.1, which means a fault 1 km long would have a rupture width/depth of 10 km). As a means of addressing these issues, <em>ccldpy</em> will compute $AR$ from the as-published $L$ and $W$ relationships. If $AR \leq$ 1.0, <em>ccldpy</em> will constrain $AR$ = 1.0, in effect modeling the rupture as a square which approximates a circular rupture. It follows that $L$ and $W$ can be derived from the published area ($A$)-relation and the constrained $AR$-relation as:

$$L = \sqrt{A * AR}$$

and

$$W = \sqrt{A / AR}$$ 

### Stable-Continental Shallow-Crustal Earthquakes:

Only one scaling-relationship for stable-continental shallow-crustal type earthquakes is implemented within <em>ccldpy</em>:
1. Leonard (2014): relationships for area, length, and width; separate coefficients for strike-slip, normal, and reverse styles-of-faulting

![ccldpy - stable scaling relations](https://github.com/tristanbuckreis/ccldpy/assets/71461454/18319da1-861d-43bd-8b06-0f0a52248ac6)
<b>Figure 2:</b> Median magnitude-scaling relationships for stable-continental shallow-crustal events.

### Subduction (Interface and Intra-Slab) Earthquakes:

Two scaling-relationships are implemented within <em>ccldy</em>:
1. Thingbaijam et al. (2017): relationships for area, length, and width; only interface earthquakes
2. Contreras et al. (2022): relationships for area and aspect ratio; separate coefficients for interface and intra-slab type events

![ccldpy - subduction scaling relations](https://github.com/tristanbuckreis/ccldpy/assets/71461454/4afcbef6-9967-45e0-8d92-fbc6ba879760)
<b>Figure 3:</b> Median magnitude-scaling relationships for subduction (interface and intra-slab) events.

# Generating Rupture Geometries:

Most scaling relationships provide estimates for rupture area ($A$), length ($L$), and width ($W$), each of which coming from a different model with quantified uncertainty ($\sigma$), however Chiou and Youngs (2008) and Contreras et al. (2022) provide relationships for aspect ratio ($AR$) to be used to derive the rupture geometry ($L$ and $W$). <em>ccldy</em> implements a consistent approach, regardless of which types of information are provided by the published scaling relations. This approach begins by computing $A$, which, for example, can be expressed as:

$$A = 10^{a_1 + b_1 M + \epsilon\sigma_1}$$  

where $M$ is the moment magnitude, $a_1$ and $b_1$ are model coefficients, and $\epsilon$ is a standard normal variate (mean = 0 and standard deviation = 1). For relationships that provide $AR$, it is simply computed using the published relationship, otherwise $L$ and $W$ are computed:

$$L = 10^{a_2 + b_2 M + \epsilon\sigma_2}$$ 

and 

$$W = 10^{a_3 + b_3 M + \epsilon\sigma_3}$$

If $L$ and $W$ were directly computed, $AR$ is checked for reasonableness. If $L/W$ > 1.0, then these dimensions are used for this particular realization of the rupture surface. Otherwise, $AR$ is constrained to be 1.0 + $\epsilon \sigma_{AR,CY08}$, where $\sigma_{AR,CY08}$ = 0.16 (Chiou and Youngs 2008). $L$ and $W$ are then derived from $A$ and $AR$.

The position of the hypocenter on the rupture surface is also randomized for each realization. Shallow-crustal and stable-continental type events use relative down-dip ($f_d$) and along-strike ($f_l$) distributions which range from 0.0 to 1.0 (inclusive) developed by Chiou and Youngs (2008), whereas subduction-type events (interface and intra-slab) use distributions developed by Contreras et al. (2022). The exact location of the hypocenter (latitude, longitude, and depth) is not randomized in this approach, only the relative position of the hypocenter on the rupture plane is. 

### Specifying Scaling-Relationships:

The scaling-relationships proposed by Wells & Coppersmith (1994), Leonard (2014), Thingbaijam et al. (2017), and Contreras et al. (2022) are self-consistent and provide all the necessary information needed to generate a rupture surface. As such, these models are implemented as separate branches within <em>ccldpy</em>. Chiou and Youngs (2008), which is only applicable for shallow-crustal type events, only provides $AR$, which is insufficient by itself to define a rupture geometry. Therefore, the Chiou and Youngs (2008) $AR$ relationship can be used with the $A$-relationships proposed by Wells and Coppersmith (1994), Leonard (2014), and Thingbaijam et al. (2017) within <em>ccldpy</em> as three separate branches. Table 1 summarizes all current branches implemented in <em>ccldpy</em> for each type of earthquake.

<b>Table 1</b>: Summary of scaling-relationship branches currently implemented in <em>ccldpy</em>.
<table>
    <thead>
        <tr>
            <th>Earthquake Type</th>
            <th>Model</th>
            <th>$A$ Relationship</th>
            <th>$L$ & $W$ or $AR$ Relationship(s)</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan=6>crustal</td>
            <td>WellsCoppersmith1994</td>
            <td>Wells & Coppersmith (1994)</td>
            <td>Wells & Coppersmith (1994)</td>
        </tr>
        <tr>
            <td>Leonard2014</td>
            <td>Leonard (2014)</td>
            <td>Leonard (2014)</td>
        </tr>
        <tr>
            <td>ThingbaijamEtAl2017</td>
            <td>Thingbaijam et al. (2017)</td>
            <td>Thingbaijam et al. (2017)</td>
        </tr>
        <tr>
            <td>ChiouYoungs2008_WellsCoppersmith1994</td>
            <td>Wells & Coppersmith (1994)</td>
            <td>Chiou & Youngs (2008)</td>
        </tr>
        <tr>
            <td>ChiouYoungs2008_Leonard2014</td>
            <td>Leonard (2014)</td>
            <td>Chiou & Youngs (2008)</td>
        </tr>
        <tr>
            <td>ChiouYoungs2008_ThingbaijamEtAl2017</td>
            <td>Thingbaijam et al. (2017)</td>
            <td>Chiou & Youngs (2008)</td>
        </tr>
        <tr>
            <td>stable</td>
            <td>Leonard2014</td>
            <td>Leonard (2014)</td>
            <td>Leonard (2014)</td>
        </tr>
        <tr>
            <td rowspan=2>interface</td>
            <td>ThingbaijamEtAl2017</td>
            <td>Thingbaijam et al. (2017)</td>
            <td>Thingbaijam et al. (2017)</td>
        </tr>
        <tr>
            <td>ContrerasEtAl2022</td>
            <td>Contreras et al. (2022)</td>
            <td>Contreras et al. (2022)</td>
        </tr>
        <tr>
            <td>intraslab</td>
            <td>ContrerasEtAl2022</td>
            <td>Contreras et al. (2022)</td>
            <td>Contreras et al. (2022)</td>
        </tr> 
    </tbody>
</table>

# Selection of the Preferred Rupture Geometry:
<em>ccldpy</em> generates a stochastic set of possible rupture surfaces given the available source metadata as outlined above. The objective is to select the most probable surface which does not result in atypically short or long finite-fault distances (e.g., $R_{RUP}$) for any given site, when published finite-fault models from the literature are not available. This is done by computing rupture distances ($R_{RUP}$) for a grid of pseudo-stations for each rupture realization (Figure 4). The most-probable rupture surface, which is ultimately selected to compute real distance metrics, is that which minimizes the squared difference between $R_{RUP}$ and the median $R_{RUP}$ at each pseudo-station. In other words, the optimal rupture surface for the purpose of computing reasonable finite-fault distances is that which minimizes the following expression:

$$\sum_{r = 1}^{N_r} \sum_{s = 1}^{N_S} (R_{RUP,median,s} - R_{RUP,r,s})^2$$

where $N_r$ and $N_s$ represent the total number of simulated rupture surfaces and pseudo-stations, respectively; $R_{RUP,r,s}$ is the rupture distance between simulated rupture $r$ and pseudo-station $s$; and $R_{RUP,median,s}$ is the median rupture distance at pseudo-station $s$ from all simulated rupture surfaces.

![ccldpy - grid of pseudo-stations](https://github.com/tristanbuckreis/ccldpy/assets/71461454/db6ad04c-1b57-41a1-a3aa-f77069abc471)
<b>Figure 4:</b> Plan view of grid of pseudo-stations distributed around a realization of a simulated rupture surface.

# Installation/Usage

The <em>ccldpy.py</em> file can be imported into Python environments. 
### Simplified instructions:

1. Download <em>ccldpy.py</em> into your {<em>desired path</em>}.
2. Import <em>ccldpy</em> into your code:<br>
```python
import sys
sys.path.append({desired path})
import ccldpy
```
# 
### Simulation Methods:

The current version of ccldpy (2.0.0) supports five methods of simulation, which are specified using the <em>method</em> indicators described below. 

A = One or two nodal plane solutions (strike, dip, and rake) are known, however the first solution is preferred. Only the area, aspect-ratio, and position of hypocenter on the rupture surface are randomized between simulation. This method is not recommended because in reality we cannot be certain of a preferred orientation.

B = Two nodal plane solutions are known, however the second solution is preferred. Only the area, aspect-ratio, and position of hypocenter on the rupture surface are randomized between simulations. This method is not recommended because in reality we cannot be certain of a preferred orientation.

C = Two nodal plane solutions are known, and neither is preferred over the other. Each simulation will randomly select which nodal plane solution to use, and the area, aspect-ratio, and position of hypocenter on the rupture surface are also randomized. This is the recommended method when two nodal plane solutions are known.

D = Only one nodal plane solution is known or assumed, but with some uncertainty. The given rake angle is used to assign rupture mechanism (<u>S</u>trike-<u>S</u>lip, <u>N</u>or<u>M</u>al dip-slip, or <u>R</u>e<u>V</u>erse dip-slip - for shallow-crustal or stable-continental type events) and the strike and dip are randomized with &#177; 30<sup>o</sup> and &#177; 10<sup>o</sup>, respectively. The area, aspect-ratio, and position of hypocenter on the rupture surface are also randomized between simulations. This method is not recommended for general applications, exceptions are when there is outstanding evidence that supports a known or preferred nodal plane solution.

E = No nodal plane solutions are known or assumed. Rake (faulting mechanism), strike [0<sup>o</sup> - 360<sup>o</sup>), and dip [0<sup>o</sup> - 90<sup>o</sup>] are randomly assigned with equal probability during each simulation. This is the recommended method when there is missing nodal plane information.

![Category Illustration](https://user-images.githubusercontent.com/71461454/220185818-708986c3-28ff-4dfa-b54b-e225dfe261f3.png)
<b>Figure 5:</b> Schematic illustration of simulation results for each simulation <em>method</em>.


# Main Function
```python
simulate_rupture_surface(eqid, eqType, region, elat, elon, hypd, magnitude, method, nsims,
                         mechanism=None,
                         strike=None, dip=None, rake=None, 
                         strike2=None, dip2=None, rake2=None)
```

### Input Parameters:

#### Required Keys:

    - eqid = unique integer identifier for the event (to set the random seed)
    - eqType = type of earthquake:
        "crustal" = shallow-crustal in active tectonic regimes
        "intraslab" = intraslab-type in subduction regimes
        "interface" = interface-type in subduction regimes
        "stable" = shallow-crustal in stable-continental regimes
    - region = geographic region where the earthquake occurred"
        "japan"
        "chile"
        "other"
    - elat = hypocenter latitude (degrees)
    - elon = hypocenter longitude (degrees)
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
              (if optional "mechanism" argument is not specified, simulations randomly assign one)
    - nsims = Number of simulations assigned to each M-scaling relationship. Total number of simulations should be odd.
        nsims[0] = Wells & Coppersmith (1994) - [recommended 334]
        nsims[1] = Leonard (2014) - [recommended 333]
        nsims[2] = Thingbaijam et al. (2017) [recommended 333]
        nsims[3] = Chiou & Youngs (2008) aspect ratio model with Wells & Coppersmith (1994) area relationship [recommended 111]
        nsims[4] = Chiou & Youngs (2008) aspect ratio model with Leonard (2014) area relationship [recommended 111]
        nsims[5] = Chiou & Youngs (2008) aspect ratio model with Thingbaijam et al. (2017) area relationship [recommended 111]
        nsims[6] = Contreras et al. (2022) - [recommended 333]
                            
#### Optional Keys:
    - mechanism = known or preferred style-of-faulting [default None]
        "SS" = strike-slip: (-180 < rake < -150) or (-30 < rake < 30) or (150 < rake < 180)
        "NM" = normal: -150 < rake < -30
        "RV" = reverse: 30 < rake < 150
    - strike = strike-angle (degrees) of the first nodal plane solution [default None]
    - dip = dip-angle (degrees) of the first nodal plane solution [default None]
    - rake = rake-angle (degrees) of the first nodal plane solution [default None]
    - strike2 = strike-angle (degrees) of the second nodal plane solution [default None]
    - dip2 = dip-angle (degrees) of the second nodal plane solution [default None]
    - rake2 = rake-angle (degrees) of the second nodal plane solution [default None]
   
#### Returns:
    - SIMULATIONS = pandas DataFrame object containing all simulated rupture surfaces
    - SELECTED = pandas DataFrame object containing the selected rupture surface and statistics

# References:

Chiou, B. S.-J. (2021). P4CF [https://github.com/bc88bc/P4CF]

Chiou ,B. S.‐J. and Youngs R. R. (2008) <em>NGA Model for Average Horizontal Component of Peak Ground Motion and Response Spectra</em>, PEER Rept. 2008/09, Pacific Earthquake Engineering Research Center, Berkeley, California.

Contreras V., Stewart J.P., Kishida T., Darragh R.B., Chiou B. S.‐J., Mazzoni S., Kuehn N., Ahdi S.K., Wooddell K., and Youngs R.R., et al. (2020) Source and path database, in <em>Data Resources for NGA‐Subduction Project</em>, Stewart J. P. (Editor), Chapter 4, PEER Rept. 2020/02, Pacific Earthquake Engineering Research Center, UC Berkeley, Berkeley, California.

Contreras, V., J.P. Stewart, T. Kishida, R.B. Darragh, B.-S.J. Chiou,S. Mazzoni, R.R. Youngs, N.M. Kuehn, S.K. Ahdi, K. Wooddell, R. Boroschek, F. Rojas, and J. Ordenes (2022) NGA-Sub source and path database. <em>Earthquake Spectra</em> 38(2): 799 - 840.

Leonard, M. (2014). Self-consistent earthquake fault-scaling relations: Update and extension to stable continental strike-slip faults. <em>Bulletin of the Seismological Society of America</em> 104(6): 2953 - 2965.

Thingbaijam K.K.S., Mai P.M., and Goda K. (2017) New empirical earthquake source-scaling laws. <em>Bulletin of the Seismological Society of America</em> 107(5): 2225 - 2246.

Veness, C. (n.d.). <em>Movable type scripts</em>. Calculate distance and bearing between two Latitude/Longitude

Wells, D.L. and K.J. Coppersmith (1994) New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. <em>Bulletin of the Seismological Society of America</em> 84(4): 974 - 1002.
