Configuration
=============

Main config
-----------

``requests`` (bool):
    Whether to use Web-based API to retrieve data.

``apiKey`` (str):
    API key for IllustrisTNG, required if ``requests`` is True.

``simulation`` (str):
    Simulation name. TNG50-1, TNG100-1, TNG300-1, etc.

``hydrodynamicSolver`` (str):
    hydrodynamic solver for gas particles, can be VoronoiMesh (TNG) or smoothParticle (EAGLE)  

``TNGPath`` (str):
    Path to simulation data in local. 

``postprocessingPath`` (str):
    Path to postprocessing data in local. Offset for each snapshot is required to read data from TNG locally.

``workingDir`` (str):
    Working directory for running SKIRT code. 

``simulationMode`` (str):
    | Simulation mode. NoMedium, ExtinctionOnly or DustEmission.
    | NoMedium is an ideal case with no medium;
    | ExtinctionOnly is employed when extinction and scattering for dust are considered;
    | DustEmission also includes secondary emission from dust.

``includeDust`` (bool):
    Whether to include dust particles. Be True if ``simulationMode`` is ExtinctionOnly or DustEmission

``dustEmissionType`` (str):
    Dust emission type. Equilibrium or Stochastic.

``dustModel`` (str):
    Dust model. ZubkoDustMix, DraineLiDustMix or ThemisDustMix.

``includeVelocity`` (bool):
    Whether to include velocities for particle data.

``minWavelength``, ``maxWavelength`` (float):
    Wavelength range in rest frame, in micron.

``boxLengthScale`` (float):
    Particles are retrieved from the box with length ``halfStellarMassRadius * boxLengthScale``.

``maxBoxLength`` (float):
    Maximum box length in kpc.

``fieldOfView`` (float)
    Field of view in arcsec, equals to box length if 0.

``wavelengthGrid`` (str):
    Wavelength grid type. Linear or Log.

``numWavelengths`` (int):
    Number of wavelengths bins.

``minLevel``, ``maxLevel`` (int):
    Octree min/max level refinement for dust calculation of SKIRT.

``numPackets`` (float):
    Number of photon packets. Main parameter to affect the SNR of SKIRT simulation.

``SEDFamily`` (str):
    SED family for quenched star particles. BC03 or FSPS.

``initialMassFunction`` (str):
    | Initial mass function for quenched star particles.  
    | Chabrier or Salpeter for BC03 SED family. Chabrier, Kroupa or Salpeter for FSPS family.

``starformingSEDFamily`` (str):
    SED family for star-forming regions. MAPPINGS (`Groves et al. 2008 <https://iopscience.iop.org/article/10.1086/528711>`_) or TODDLERS (`Kapoor et al. 2023 <https://academic.oup.com/mnras/article/526/3/3871/7287615>`_, `2024 <https://www.aanda.org/articles/aa/full_html/2024/12/aa51207-24/aa51207-24.html>`_).

``minStellarMass``, ``maxStellarMass`` (float):
    Stellar mass range for subhalos in Msun, inf for infinite.

``faceAndEdge`` (bool):
    Whether to use face-on and edge-on angles derived by angular momentum.

``numViews`` (int):
    Number of instrument views (observing directions).

``randomViews`` (bool):
    Whether to use random viewing angles.

``inclinations``, ``azimuths`` (list):
    Inclinations and azimuths for instrument views.

``postProcessing`` (bool):
    Whether to perform postprocessing.

``spatialResol`` (float):
    Base spatial resolution in pc, must be provided if postprocessing is False.

``imageUnit`` (str):
    Unit type of bandpass image output. electron or flux.

``surveys`` (list[str]):
    Surveys considered. Euclid, JWST, HST, CSST, HSC are built-in.

``displaySED`` (bool):
    Whether to display SED.

``displaySEDxlogscale`` (bool):
    Whether to display SED in logscale for wavelength.

``outputSEDOnly`` (bool):
    If only output SED.

``snapNum`` (int):
    Snapshot ID.

``snapRedshift`` (float):
    Snapshot redshift of Snapshot ID.
    
``viewRedshift`` (float):
    | viewing redshift, should be close to snapRedshift.  
    | For generating galaxies at continous redshift with a small offset on snapRedshift.  
    | For example, viewRedshift can be 0.105 for snapshot-91 at 0.1.

``inLocal`` (bool):
    if use local cosmology.

``viewDistance`` (float):
    viewing distance, in Mpc, must be provided if inLocal is True.

``numThreads`` (int):
    Number of threads. No speedup for threads larger than 24.

``recordComponents`` (bool):
    | Whether to record individual components including transparent, primary direct, primary scattered  
    | and secondarytransparent, secondarydirect, secondaryscattered, apart from total components.  
    | Memory cosumptions will be 7 times if True.

``ageThreshold`` (float):
    Age threshold for discriminating star-forming regions and quenched star particles, in Myr.

``logCompactnessMean``, ``logCompactnessStd`` (float):
    | Only used when starformingSEDFamily is MAPPINGS.  
    | Logarithmic mean and standard deviation of compactness for star-forming star particles, see `Kapoor et al. 2021 <https://academic.oup.com/mnras/article/506/4/5703/6324023>`_.

``logPressure`` (float):
    | Only used when starformingSEDFamily is MAPPINGS.  
    | Logarithmic pressure for star-forming star particles. log10[(Pressure/k_B)/cm^-3 K] = logPressure

``constantCoveringFactor`` (bool):
    | Only used when starformingSEDFamily is MAPPINGS.  
    | Whether to use constant covering factor.

``coveringFactor`` (float):
    | Only used when starformingSEDFamily is MAPPINGS.  
    | Constant covering factor, if constantCoveringFactor is True, see `Groves et al. 2008 <https://iopscience.iop.org/article/10.1086/528711>`_.

``PDRClearingTimescale`` (float):
    | Only used when starformingSEDFamily is MAPPINGS.  
    | PDR clearing timescale, in Myr, if constantCoveringFactor is False, see `Baes et al. 2024 <https://www.aanda.org/articles/aa/full_html/2024/03/aa48418-23/aa48418-23.html>`_.

``temperatureThreshold`` (float):
    Temperature threshold for creating dusts from gas particles, in K.

``pahfraction`` (float):
    | Only used when starformingSEDFamily is TODDLERS.  
    | Fraction of PAH mass in dust.

``sedMode`` (str):
    | Only used when starformingSEDFamily is TODDLERS.  
    | SED mode, Cloud or SFRNormalized.

``stellarTemplate`` (str):
    | Only used when starformingSEDFamily is TODDLERS.  
    | Stellar template, SB99Kroupa100Sin, BPASSChab100Bin or BPASSChab300Bin.

``sfrPeriod`` (int):
    | Only used in SFRNormalized mode.  
    | SFR period in Myr, can be 10 or 30.

``starFormationEfficiency`` (float):
    | Only used when starformingSEDFamily is TODDLERS.  
    | star formation efficiency.

``cloudNumberDensity`` (float):
    | Only used when starformingSEDFamily is TODDLERS.  
    | Cloud number density, in cm^-3.

``alpha`` (float):
    | Only used when starformingSEDFamily is TODDLERS.  
    | Alpha for power-law distribution of cloud mass.

``scaling`` (float):
    | Only used when starformingSEDFamily is TODDLERS.  
    | Scaling factor.

``massFraction`` (float):
    Fraction of the metallic gas locked up in dust.

``DISMModel`` (str):
    Dust ISM recipe. `Camps_2016 <https://academic.oup.com/mnras/article/462/1/1057/2589990>`_ or `Torrey_2012 <https://academic.oup.com/mnras/article/427/3/2224/1099996>`_.

``numSilicateSizes``, ``numGraphiteSizes``, ``numPAHSizes``, ``numHydrocarbonSizes`` (int):
    Number of bins for dust grains.


Survey config
-------------

``filters`` (list[str]):
    Considered filters for survey.

``resolFromPix`` (bool):
    Whether to use resolution derived from pixel scale.

``resolution`` (float):
    Spatial resolution, in pc, must be provided if resolFromPix is False, override spatialResol in config.ini.

``pixelScales`` (float, list[float]):
    Pixel scale for considered filters, in arcsec.

``numExposure`` (float, list[float]):
    Number of exposure for considered filters.

``exposureTime`` (float):
    Exposure time, in second.

``aperture`` (float):
    Aperture size for instrument, in meter.

``includePSF`` (bool):
    Whether to include PSF effects.

``PSFFromFile`` (bool):
    Whether to use PSF from file.

``PSFFWHM`` (float, list[float]):
    FWHM of PSF, in arcsec.

``includeBkg`` (bool):
    Whether to include background.

``noiseType`` (str):
    | By what parameters noise levels are calculated, can be instrument or limitingMagnitude.  
    | For noiseType == instrument, skyBkg, darkCurrent, readOut are required.  
    | For noiseType == limitingMagnitude, limitMag, limitSNR, limitAperture and zeroPoint are required.

``skyBkg`` (list[float]):
    Background level for considered filters. Please refer to notebook `calc_sky_bkg.ipynb <https://github.com/xczhou-astro/galaxyGenius/blob/main/Bkg/calc_sky_bkg.ipynb>`_ for calculating skyBkg.

``darkCurrent`` (float, list[float]):
    Dark current.

``readOut`` (float, list[float]):
    Readout noise.

``limitMag`` (float, list[float]):
    limiting magnitudes.

``limitSNR`` (float, list[float]):
    SNR for the limiting magnitude.

``limitAperture`` (float, list[float]):
    The aperture used to measure the limiting magnitude, in arcsec. 

``zeroPoint`` (float): 
    Zero point for magnitude system.

``imgDisplay`` (bool):
    Whether to display image.

``RGBImg`` (bool):
    | Whether to create and display RGB image. 
    | RGB image is created from `astropy.visualization.make_rgb`, a new feature added in version 7.0.0.
    | The RGB image may not be as one expected, feel free to edit the `convert_to_rgb` function in `utils.py`.

``RGBFilters`` (list[str, str, str]):
    Considered three filters for RGB image.

``displayFilter`` (str):
    Filter for displaying image if RGBImg is False.