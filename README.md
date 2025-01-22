# galaxyGenius

A Mock Galaxy Image Generator for Various Telescopes from Hydrodynamical Simulations

## Results  
### SEDs
![SEDs](assets/fiducial_SED.png)  
### HSC
![HSC](assets/HSC_combined.png)  
### HST WFC3
![HST](assets/HST_combined.png)  
### JWST NIRCam
![JWST](assets/JWST_combined.png)  
### CSST
![CSST](assets/CSST_combined.png)  

## Dependence
Python verison:  
`python==3.11.10`  

Python packages:
`tomlkit==0.13.2`  
`numpy==1.24.0` ⚠️  
`scipy==1.13.1`  
`matplotlib==3.9.2`  
`astropy==7.0.0` ⚠️  
`scikit-image==0.24.0`  
`joblib==1.4.2`  
`matplotlib_scalebar==0.8.1`  
`h5py==3.9.0` ⚠️  
`termcolor==2.4.0`  
`photutils==2.0.2`  

⚠️: Package versions should strictly match, otherwise there may be compatibility issues.

Other packages:  
[`illustris_python`](https://github.com/illustristng/illustris_python)  
[`SKIRT9`](https://skirt.ugent.be/root/_installation_guide.html)  

A new Python environment is recommended  

```Python
conda create -n galaxyGenius python=3.11.10  

conda activate galaxyGenius  

pip install .  
```  
`Data` directory should be set in environment, named `GALAXYGENIUS_DATA_DIR` or use
```Bash
export GALAXYGENIUS_DATA_DIR=/path/to/Data
```

## Recommended folder tree
```Bash
├── galaxyGenius/  
│   ├── Data/  
│   ├── galaxyGenius/  
│   ├── Notebooks/  
│   └── (workspace)/  
├── TNG100/  
│   ├── groups_N/  
│   │   ├── fof_subhalo_tab_N.n.hdf5  
│   │   └── ...  
│   └── snapdir_N/  
│       ├── snap_N.n.hdf5  
│       └── ...  
└── postprocessing/  
    └── offsets/  
        ├── offsets_N.hdf5  
        └── ...  
```
We also implement [Web-based API](https://www.tng-project.org/data/docs/api/) in `galaxyGenius`, therefore, users do not even need to download the large snapshots and groups files. This feature can be activated by setting `requests=true` and provide `apiKey` in config. However, the generation speed will be slower and depends on the network speed.  

## Filters and PSFs
Filter throughputs and PSFs of specific surveys are saved in `Data/filters/(survey)` and `Data/PSFs/(survey)` directories.  
The format of filters are recommended to be plain text file with extension `.fil`, including two columns: wavelength (in angstrom) and throughputs.  
PSFs are recommended to be in `numpy.array`, and can be opened by `numpy.loadtxt` or `numpy.load`.  
Please make sure that the filters and PSFs exist in correct directories and formats before running galaxyGenius.

## Sky backgrounds  
Sky background noises are calculated based on the throughput of each filter and sky emission curve.  
We provide a helper notebook for calculation of instrumental noises in [Notebooks/calc_sky_bkg.ipynb](https://github.com/xczhou-astro/galaxyGenius/blob/main/Notebooks/calc_sky_bkg.ipynb).  

## Usage
### Initialization
```Python
python init.py --workspace=workspace --surveys="CSST,HST"
```
`config.ini`, `config_CSST.ini` and `config_HST.ini` will be created in `workspace` directory. Please freely edit them as wish.  

if `surveys` are not specified, `postPostprocessing` will be set as False.

Currently, throughputs and PSFs of filters for CSST, HST, JWST, Roman and HSC are all uploaded. For Euclid, we only upload the filters, since the PSF files cannot be accessed.  

### Run
Enter workspace, and create a python file named `run.py`
```Python
# run.py

from galaxyGenius.config import Configuration
from galaxyGenius.generation import DataGeneration
from galaxyGenius.preprocess import PreProcess
from galaxyGenius.postprocess import PostProcess

config = Configuration() # initialize Configuration class
conf = config.get_config() # read config from current directory

preprocess = PreProcess(conf) # initialize PreProcess class
subhalos = preprocess.get_subhalos() # get subhalos 
subhaloIDs = subhalos['subhaloIDs'] # get subhaloIDs

for i, ID in enumerate(subhaloIDs):

    # Dynamically set parameters
    data = {}
    data['numPackets'] = 1e6

    if i < 100:
        data['inclinations'] = [0]
        data['azimuths'] = [0]
        data['faceAndEdge'] = False
    else:
        data['faceAndEdge'] = True


    preprocess.subhalo(subhaloID=ID) # get properties of subhalo
    preprocess.prepare(data) # prepare for simulation, including retrieving particles and creating .ski file

    dataGeneration = DataGeneration(config=conf) # initialize DataGeneration class
    dataGeneration.runSKIRT() # run SKIRT

    postprocess = PostProcess(subhaloID=ID, config=conf) # initialize PostProcess class
    postprocess.runPostprocess() # run postprocessing
```  
then `python run.py`. 

Or you can interactively run in jupyter as illustrated in [Notebooks/tutorial.ipynb](https://github.com/xczhou-astro/galaxyGenius/blob/main/Notebooks/tutorial.ipynb).  

### For hydrodynamical simulations except TNG
Necessary particle data and properties need to be created and saved in `workingDir (run)` in advance, and call `preprocess.inputs(data)` and subsequent methods to perform the galaxy generation. An example for EAGLE simulation is provided in [EAGLE/eagle.ipynb](https://github.com/xczhou-astro/galaxyGenius/blob/main/EAGLE/eagle.ipynb).  
data: `dict`  
`snapRedshift`, `cosmology`, `stellarMass`, `subhaloID`, `boxLength` keys are required. List of `inclinations` and `azimuths` are optional.


### Outputs
`./dataCubes` includes the IFU-like raw dataCubes and SEDs of **SubhaloID** generated by SKIRT and auxiliary files.  
`./mock_CSST` includes the bandpass images and SEDs of **SubhaloID** for CSST.  
`./mock_HST` includes the bandpass images and SEDs of **SubhaloID** for HST.  

Folder tree:
```bash
dataCubes/
└── Subhalo_(ID)/
    ├── config.json # recorded config for dynamical modification
    ├── dusts.txt # dust particles from gas particles
    ├── properties.json # properties of the subhalo and surveys
    ├── skirt_log.txt # log file of SKIRT
    ├── skirt_parameters.xml # SKIRT execution file 
    ├── skirt_view_(ViewID)_sed.dat # Generated SEDs
    ├── skirt_view_(ViewID)_total.fits # Generated dataCubes
    ├── starforming_regions.txt # starforming regions
    └── stars.txt # star particles
```

```bash
mock_(survey)/
└── Subhalo_(ID)/
    ├── galaxy_images.fits # bandpass images 
    ├── galaxy_SED_(ViewID).png # Plot of SEDs
    ├── galaxy_SEDs.fits # SEDs
    └── galaxy_view_(ViewID).png # Plot of bandpass images
```

Bandpass images are saved in pages of fits file, and each page includes images in different views.

## Classes
### Configuration
```Python
config = Configuration(surveys=None)
```
surveys: `str (N,)` or `list`, considered surveys. If None, configurations will be read from current directory, or configuration files will be created from templates.  
```Python
conf = config.get_config()
```  
return configurations.  
```
config.add_survey(surveys)
```
surveys: `str (N,)` or `list`, new surveys to be added. Call `get_config()` to update configurations.  
```
config.remove_survey(surveys)
```
surveys: `str (N,)` or `list`, surveys will be removed. Call `get_config()` to update configurations.  
```
config.save_config()
```
manually save configurations to `config.ini`.  
```
config.check_config()
```
manually check configurations.  

### PreProcess
```Python
preprocess = PreProcess(config)
```  
```Python
preprocess.get_subhalos()
```  
Subhalos in stellar mass range in `config.ini` will be read.  
```Python
subhaloIDs = preprocess.get_subhaloIDs()
```  
return subhaloIDs for subhalos obtained in `get_subhalos()`.  
```Python
stellarMasses =preprocess.get_stellarMasses()
```
return stellar masses for subhalos obtained in `get_subhalos()`.
```Python
preprocess.subhalo(subhaloID)
```
`subhaloID`: `int`: Intialize simulation for subhalo with subhaloID.  
```Python
preprocess.prepare()
```
Preparation for simulation, including retrieving particles and creating .ski file.  
```Python
properties = preprocess.get_properties()
```
Return properties used for postprocessing.  
```Python
preprocess.runSKIRT()
```
Run SKIRT.  

### PostProcess
```Python
postprocess = PostProcess(config)
```  
`config`: `dict`: configurations.  
```Python
postprocess.runPostProcess(showImages=False)
```  
`showImages`: `bool`: if show galaxy images and SEDs, should be used in Notebooks. Run postprocessing for galaxy simulation.  

## Config.ini
`dataDir`:  
`str`, Data directory of galaxyGenius.  

`filePath`:  
`str`, Directory for TNG simulation.  

`workingDir`:  
`str`, Directory for execution of SKIRT.  

`simulationMode`:  
`str`, `ExtinctionOnly` or `DustEmission`, Simulation Mode, `DustEmission` should be used when mid-infrared or far-infrared bands are considered.  

`includeDust`:  
`bool`, If include dusts generated from Gas particles; must be True if `simulationMode=DustEmission`  

`dustEmissionType`:  
`str`, `Equilibrium` or `Stochastic`, Dust emission type.  

`dustModel`:  
`str`, `ZubkoDustMix`, `DraineLiDustMix` or `ThemisDustMix`, Dust model.  

`minWavelength`, `maxWavelength`:  
`float`, Considered wavelength range in rest frame, in micron; `maxWavelength` should be higher than maximum wavelength of filters considered.  

`smoothLengthFromNeigbors`:  
`bool`, If smoothing length are derived from neigbors as illustrated in [Baes et al. 2024](https://www.aanda.org/articles/aa/full_html/2024/03/aa48418-23/aa48418-23.html).  

`NthNeigbor`:  
`int`, Nth nearest neighbor for smoothing length, if smoothLengthFromNeighbors is True.  

`maxSmoothingLength`:  
`float`, maximum smoothing length, if smoothLengthFromNeighbors is True, in kpc.  

`boxLengthScale`:  
`float`, Determine the boxsize to retrieve particles;  

`maxBoxLength`:  
`float`, Maximum boxsize, in kpc; `boxisze = min(halfStellarMassRadius * boxLengthScale, maxBoxLength)`.  

`wavelengthGrid`:  
`str`, `Linear` or `log`, Wavelength grid type for SKIRT.  

`numWavelengths`:  
`int`, Wavelength bins for SKIRT.  

`minLevel`, `maxLevel`:  
`int`: Octree min/max level refinement for dust calculation.  

`numPackets`:  
`float`: number of photon packets launched during simulation; determine the level of signal to noise in the results.  

`SEDFamily`:  
`str`, `BC03` or `FSPS`, SED family for quenched star particles.  

`initialMassFunction`:  
`str`, `Chabrier` or `Salpeter` for `BC03`; `Chabrier`, `Salpeter` or `Kroupa` for `FSPS`, Intial mass function for SED family.  

`minStellarMass`, `maxStellarMass`:  
`float` or `inf`, Stellar mass range for subhalos, in 1e10 Msun; `inf` for infinite.  

`numViews`:  
`int`, Number of instrument views for datacube generation by SKIRT.  

`randomViews`:  
`bool`, If specify views by randoms.  

`inclinations`, `azimuths`:  
`float (numViews,)`, Views; must be provided if `randomViews=False`; inclinations: 0 ~ 180, azimuths: -360 ~ 360.  

`FoVboxLengthRatio`:  
`float`, Ratio for field of view, `Fov = Boxsize * FoVboxLengthRatio`.  

`postProcessing`:  
`bool`, If run postprocessing.  

`saveDataCube`:  
`bool`, If save data cubes; must be True if `postProcessing=False`.  

`spatialResol`:  
`float`, Spatial resolution for data cube, in pc; must be provided if `postProcessing=False`.  

`imageUnit`:  
`str`, `electron`, `flux` or `magnitude`, unit type of image output.  

`surveys`:  
`str (N,)`, Considered surveys; must be provided if `postProcessing=True`.  

`displaySED`:  
`bool`, If display SED.  

`displaySEDxlogscale`:  
`bool`, If display SED in logscale for wavelength.  

`snapNum`:  
`int`,  Snapshot ID of TNG simulation.  

`fixedRedshift`:  
`float`, Redshift of snapshot ID; should be larger than 0 when `snapNum=99` to avoid error.  

`numThreads`:  
`int`, Number of Threads to run SKIRT. The perforamnce will not be improved if `numThreads > 24`.  

`recordComponents`:  
`bool`, If including transparent, primary direct, primary scattered, secondarytransparent, secondarydirect, secondaryscattered components for data cube and SED; Consumption of memory increases to 7 times if True.  

`ageThreshold`:  
`float`, Age threshold for discriminating star-forming and quenched star particles, in Myr.  

`logCompactnessMean`, `logCompactnessStd`:  
`float`,  logCompactness for star-forming particles, sampled from normal distribution, see [Kapoor et al. 2021](https://academic.oup.com/mnras/article/506/4/5703/6324023).  

`logPressure`:  
`float`, `log10[(Pressure/k_B)/cm^-3 K] = logPressure` for star-forming particles.  

`constantCoveringFactor`:  
`bool`, If use constant covering factor.  

`coveringFactor`:  
`float`, constant covering factor, if `constantCoveringFactor=True`, see [Groves et al. 2008](https://iopscience.iop.org/article/10.1086/528711).  

`PDRClearingTimescale`:  
`float`, covering factor is `f = e^(-t / PDRClearingTimescale)`, where t is the age, if `constantCoveringFactor=False`, see [Baes et al. 2024](https://www.aanda.org/articles/aa/full_html/2024/03/aa48418-23/aa48418-23.html).  

`DISMModel`:  
`str`, `Camps_2016` or `Torrey_2012`, dust-containing ISM (DISM) model, see [Camps et al. 2016](https://academic.oup.com/mnras/article/462/1/1057/2589990) and [Torrey et al. 2012](https://academic.oup.com/mnras/article/427/3/2224/1099996).  

`temperatureThreshold`:  
`float`, Gas particles lower than `temperatureThreshold` will be considered as dusts, for `DISMModel=Camps_2016`.  

`massFraction`:  
`float`, Fraction of the metallic gas locked up in dust.  

`numSilicateSizes`, `numGraphiteSizes`, `numPAHSizes`, `numHydrocarbonSizes`:  
`int`, Number of bins for dust grains.  

## Config_\[survey\].ini
Config_\[survey\].ini is generated if `postProcessing=True` and `surveys` are provided.  

`filters`:  
`str (N,)`, Considered filters for survey.  

`resolFromPix`:  
`bool`, If use resolution derived from pixel scales; instrumental effects can only be added if `resolFromPix=True`.  

`resolution`:  
`float`, Spatial resolution, in pc; override `spatialResol` in config.ini.  

`pixelScales`:  
`float (N,)` or `float (1,)`, Pixel scales for considered filters, in arcsec; provide one value if homegenous, otherwise must be consistent with filters.  

`numExposure`:  
`int (N,)` or `int (1,)` Number exposures or number of filters.  

`exposureTime`:  
`float (N,)` or `int (1,)`, Exposure time, in second.  

`aperture`:  
`float`, Aperture size for instrument, in meter.  

`includePSF`:  
`bool`, If include PSF effects; can be True if `resolFromPix=True`, otherwise must be False.  

`PSFFromFile`:  
`bool`, If PSFs are from files; pixel scales of PSF array must be consistent with `pixelScales`.  

`PSFFWHM`:  
`float (N,)` or `float (1,)`, FWHM to create Gaussian PSF kernels, in arcsec.  

`includeBkg`:  
`bool`, If include background; can be True if `resolFromPix=True`, otherwise must be False.  

`gaussianNoise`:  
`bool`, if add instrumental noise as Gaussian distribution.  

`skyBkg`, `darkCurrent`, `readOut`:  
`float (N,)` or `float (1,)`: sky backgrounds, dark currents, and read-out noise, required for calculating instrumental noise.  
Please refer to notebook [Notebooks/calc_sky_bkg.ipynb](https://github.com/xczhou-astro/galaxyGenius/blob/main/Notebooks/calc_sky_bkg.ipynb) for calculating `skyBkg`.  

`imgDisplay`:  
`bool`, If display galaxy images.  

`RGBImg`:  
`bool`, If display rgb image created by 3 bands.  
RGB image is created from `astropy.visualization.make_rgb`, a new feature added in version 7.0.0.  
The RGB image may not be as one expected, feel free to edit the `convert_to_rgb` function in `utils.py`.  

`RGBFilters`:  
`str (3,)`, 3 bands to create rgb image.  

`displayFilter`:  
`str`, Band image to display if `RGBImg=False`.  
