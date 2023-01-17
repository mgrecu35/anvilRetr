### Synergistic retrievals of ice in high clouds from lidar, Ku-band radar and submillimeter wave radiometer observations
**_Mircea Grecu and John Yorks_**

### 1.Introduction ###
The future NASA Atmospheric Observing System (AOS) mission (Braun 2022) is expected to feature new combinations of observations that may be used to quantify the amounts of ice in high clouds and characterize the microphysical properties of ice particles. These observations include lidar backscatter, Ku-band radar reflectivity and submillimeter wave radiometer brightness temperature measurements.  While not optimal for cloud ice estimation, but for the characterization of a broader spectrum of cloud and precipitation processes, these observations are nevertheless synergistic from the cloud ice perspective. That is, despite the fact that lidar observations attenuate quickly in thick ice clouds and the Ku-band radar will not be able to detect clouds characterized by an echo weaker than 8.0 dBZ, the active observations are expected to provide context that may be incorporate into the radiometer retrievals. In this study, we investigate the impact of incorporating the lidar and radar observations into the radiometer retrieval of ice clouds. Because the existent amount of coincident backscatter lidar, Ku-band radar, and submillimeter wave radiometer observations is rather insufficient to derive conclusive results, we employ accurate physical models to simulate the lidar, radar and radiometer observations and use a cross-validation methodology to characterize the retrieval accuracy. As retrievals from passive instrument observations strongly depend on the "a priori" information (Rodgers 2000), for the results to be relevant in real applications it is necessary to base them on realistic vertical distributions of ice properties.  Such distributions may be derived from cloud-resolving-model (CRM) simulations (Pfreundschuh et al. 2020) or directly from observations.  In this study, we employ the latter approach, as CRMs may still be deficient in properly reproducing the vertical distribution of ice clouds and their associated microphysical properties.

### 2. Methdology ###
As previously mentioned, we use CloudSat (CS) observations to derive the vertical distributions of ice properties needed in the investigation (Stephens et 2002).  Although research quality CS cloud ice products exist, to maximize the physical consistency of the approach, we do not use them and derive ice amounts and associated properties from CS reflectivity observations.  This ensures the consistency between the particle distribution assumptions and the electromagnetic scattering properties used in the CS reflectivity processing and those used the simulation of the lidar, Ku-band radar and radiometer observations.  

#### 2.1. Assumptions and forward models
To quantify the number of ice particles in an elementary atmospheric volume as a function of their size, we use normalized gamma functions (Bringi et al. 2003).  The benefit of normalized gamma functions is that encapsulate the variability of Ice Water Content (IWC) - reflectivity relationship into a single parameter, i.e. the normalized Particle Size Distribution (PSD) intercept (Testud et al. 2001; Bringi et al. 2003). The normalized PSD intercept is defined as $N_w=\frac {4^4} {\pi \rho_w} \frac {IWC} {D_m^4}$, where $IWC$ is the ice water content associated with the PSD, and $D_m$ is the mass weighted mean diameter.  Testud et al. (2001) showed that the variability in IWC reflectivity (Z) relationships may be fully explained by variability in $N_w$, and that formulae of the type 

$IWC=N_w^{1-b}aZ^b$

perfectly explain the relationships between IWC and Z calculated from observed PSDs. Equation (1) is not sufficient to derive accurate, unbiased estimates of ice water contents, because $N_w$ varies considerably in time and space. Nevertheless, multiple studies showed that it is beneficial to parameterize $N_w$ as a function of temperature (e.g. Hogan et al. 2006; Delanoe and Hogan 2008; Deng et al. 2010).  In this study we parameterize $N_w$ as a function of temperature based on the CloudSat 2C-ICE product (Deng et al. 2010; Deng et al. 2013).  Specifically, we cluster, based on their similarity, a large set 2C-ICE profiles into 36 classes using a k-Means procedure. The mean IWC profiles associated with the 36 classes are shown in continuous lines in Fig. 1.  Alternative estimates, based on the PSD assumptions and electromagnetic scattering calculations that enable accurate and physically consistent simulations of radar observations at Ku-band and radiometer observations of submillimeter-wave frequencies are also shown in Fig. 2 (symbolized as SSRG). Details regarding the estimations process are provided in the subsequent paragraphs.  As apparent in Fig. 1, the CS and SSRG estimates are in good agreement.  Some discrepancies due to discrepancies in the SSRG $N_w$ parameterization and the CS 2C-ICE "a priori assumptions" are also apparent, but they are not deemed critical in this study, whose objective is the investigation of synergistic lidar, Ku-band radar and submillimeter-wave radiometer retrievals, because the outcome is not
likely to be sensitive to such details.

![Cloud and SSRG-based IWC retrievals clustered by distribution.](./iwcClasses.png)
Fig.1 Conditional average IWC profiles for each of the 36 classes determined by the k-Means algorithm. The CS 2C-ICE profiles (Deng et al. 2013) are shown in continuous line, while the profiles derived in this are shown using the "*" symbol and referred to the SSRG estimates due to the use the self-similarity Rayleigh-Gans approximation (SSRGA) theory in the estimation process (see the text for details).  The relative height is defined relative to the freezing level.

 One may notice that the average IWC profiles in Fig.1 are characterized by different peak values and locations.  This facilitates a simple way to reverse-engineer to (some extent) the "a priori" assumptions used in the CS 2C-ICE product and use them in formulation of the type described in Eq. (1).  Specifically, the conditional derivation of relationships of the type $IWC=a_i Z^{b_i}$ for every class i and the representation of $a_i$ as a function of the height of class IWC peak. Shown in Fig. 2 is a relative height - $a_i$ scatter plot.  As apparent in the figure, and as expected, the $a_i$ exhibits a strong variation with the relative height. Coefficient $b_i$ exhibits a height dependency as well, but the range of variation is significantly smaller, almost zero relative to the mean value of $b$. Given that any deviation of the multiplicative coefficient in an IWC-Z relation from an average is equivalent to a deviation of the associated $N_w$ from its mean value (Testud et al. 2001), the relative height $a$ relationship may be converted into a relative-height $N_w$ relationship.  We, therefore, use the data in Fig. 2 to parameterize $N_w$ as a function of height.

 ![](./iwcCoeff.png)
Fig. 2.

 For the determination of the reference $a$ and $b$ values, we assume that PSD are normalized gamma distributions with $N_w=0.08cm^{-4}$ and $mu=2$. The mass weighted mean diameter is equidistantly sampled to span the entire range of IWC values in the CS 2C-ICE dataset. The assumed mass-size relation is that of Brown and Francis (1995) because it works well with the self-similarity Rayleigh-Gans approximation (SSRGA) of Hogan et al. (2017) upon the electromagnetic scattering property calculations are based (Heymsfield et al. 2022).  The open source software scatter-1.1 of Hogan (2019) is used to provide the actual scattering properties. The SSRGA theory was developed for millimeter-wave calculations and may not be applicable at lidar's wavelength.  Therefore, for lidar calculations we are using the Mie solution of scatter-1.1. Although more accurate calculations based on more realistic ice particle shapes exist, they are rather incomplete and not readily available.  Moreover, Wagner and Deleny (2022) compared lidar backscatter observations with backscatter calculations based on coincident PSD observations and the Mie solution and found good agreement, which suggests that electromagnetic properties derived from Mie calculations are adequate for practical applications. The lidar molecular backscatter and extinction are calculated using the lidar module of COSP stands for CFMIP Observation Simulator Package (COSP; Bodas-Salcedo et al. 2011).  To account for multiple-scattering in the lidar observations, we are using the multiscatter-1.2.11 model (Hogan 2015) of Hogan and Battaglia (2008).

### 2.2 Estimation and evaluation

## 3. Results


### _References_ ###

* Bodas-Salcedo, A., Webb, M. J., Bony, S., Chepfer, H., Dufresne, J., Klein, S. A., Zhang, Y., Marchand, R., Haynes, J. M., Pincus, R., & John, V. O. (2011). COSP: Satellite simulation software for model assessment, Bulletin of the American Meteorological Society, 92(8), 1023-1043. Retrieved Jul 3, 2022, from https://journals.ametsoc.org/view/journals/bams/92/8/2011bams2856_1.xml

* Braun and co-authors, 2022. Aerosol, Cloud, Convection, and Precipitation (ACCP) Science and Applications https://aos.gsfc.nasa.gov/docs/ACCP_Science_Narrative-(March2022).pdf

* Bringi VN, Chandrasekar V, Hubbert J, Gorgucci E, Randeu WL, Schoenhuber M. 2003.  Raindrop size distribution in different climatic regimes from disdrometer and dual-polarized radar analysis. Journal of the atmospheric sciences. 2003;60(2):354-65.

* Brogniez, H., Roca, R., Auguste, F., Chaboureau, J.P., Haddad, Z., Munchak, S.J., Li, X., Bouniol, D., Dépée, A., Fiolleau, T. and Kollias, P., 2022. Time-delayed tandem microwave observations of tropical deep convection: Overview of the C2OMODO mission. Frontiers in Remote Sensing.

* Brown, P. R. A., and P. N. Francis, 1995: Improved measurements of ice water content in cirrus
using a total-water probe. J. Atmos. Oceanic Technol., 12, 410–414.

* Delanoë, J., and R. J. Hogan (2008), A variational scheme for retrieving icecloud properties from combined radar, lidar, and infrared radiometer,J. Geophys. Res.,113, D07204, doi:10.1029/2007JD009000.

* Deng, M., G. G. Mace, Z. Wang, and E. Berry, 2015: CloudSat 2C-ICE product update with a new Ze parameterization in lidar-only region, J. Geophys. Res. Atmos., 120, 12,198-12,208, doi:10.1002/2015JD023600.

* Deng M, Gerald G. Mace, Zhien Wang, and R. Paul Lawson, 2013: Evaluation of Several A-Train Ice Cloud Retrieval Products with In Situ Measurements Collected during the SPARTICUS Campaign. J. Appl. Meteor. Climatol., 52, 1014–1030.

* Deng, M., G. G. Mace, Z. Wang, and H. Okamoto, 2010: Tropical Composition, Cloud and Climate Coupling Experiment validation for cirrus cloud profiling retrieval using CloudSat radar and CALIPSO lidar, J. Geophys. Res., 115, D00J15, doi:10.1029/2009JD013104 (6) 

* Grecu, M., Tian, L., Heymsfield, G. M., Tokay, A., Olson, W. S., Heymsfield, A. J., & Bansemer, A. (2018). Nonparametric Methodology to Estimate Precipitating Ice from Multiple-Frequency Radar Reflectivity Observations, Journal of Applied Meteorology and Climatology, 57(11), 2605-2622. Retrieved Jan 16, 2023, from https://journals.ametsoc.org/view/journals/apme/57/11/jamc-d-18-0036.1.xml

* Heymsfield, A., Bansemer, A., Heymsfield, G., Noone, D., Grecu, M., & Toohey, D. (2022). Relationship of Multiwavelength Radar Measurements to Ice Microphysics from the IMPACTS Field Program, Journal of Applied Meteorology and Climatology (published online ahead of print 2022). Retrieved Jan 17, 2023, from https://journals.ametsoc.org/view/journals/apme/aop/JAMC-D-22-0057.1/JAMC-D-22-0057.1.xml

* Hogan, R. J., Mittermaier, M. P., & Illingworth, A. J. (2006). The Retrieval of Ice Water Content from Radar Reflectivity Factor and Temperature and Its Use in Evaluating a Mesoscale Model, Journal of Applied Meteorology and Climatology, 45(2), 301-317. Retrieved Jan 16, 2023, from https://journals.ametsoc.org/view/journals/apme/45/2/jam2340.1.xml

* Hogan, R.J., Honeyager, R., Tyynelä, J. and Kneifel, S., 2017. Calculating the millimetre‐wave scattering phase function of snowflakes using the self‐similar Rayleigh–Gans Approximation. Quarterly Journal of the Royal Meteorological Society, 143(703), pp.834-844.

* Hogan, R.,J., 2019. scatter-1.1. Retrieved from http://www.met.reading.ac.uk/clouds/ssrga/scatter-1.1.tar.gz

* Pfreundschuh, S., Eriksson, P., Buehler, S.A., Brath, M., Duncan, D., Larsson, R. and Ekelund, R., 2020. Synergistic radar and radiometer retrievals of ice hydrometeors. Atmospheric Measurement Techniques, 13(8), pp.4219-4245.

* Rodgers, C.D., 2000. Inverse methods for atmospheric sounding: theory and practice (Vol. 2). World scientific.

* Stephens, G.L., Vane, D.G., Boain, R.J., Mace, G.G., Sassen, K., Wang, Z., Illingworth, A.J., O'connor, E.J., Rossow, W.B., Durden, S.L. and Miller, S.D., 2002. The CloudSat mission and the A-Train: A new dimension of space-based observations of clouds and precipitation. Bulletin of the American Meteorological Society, 83(12), pp.1771-1790.

* Testud J, Oury, S, Black RA, Amayenc P,  Dou X. 2001.  The Concept of “Normalized” Distribution to Describe Raindrop Spectra: A Tool for Cloud Physics and Cloud Remote Sensing. 
Journal of Applied Meteorology. 2001; 40(6): 1118--1140.
