import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck15
import sys
import os
import shutil
from astropy.visualization import ManualInterval, LogStretch, make_rgb
from typing import Union
from types import NoneType
import logging
import inspect
import json
import re

class Units:
    _instance = None
    _initialized = False
    
    def __new__(cls, cosmology=None, snapRedshift=None, reinitialize=False):
        if cls._instance is None:
            cls._instance = super(Units, cls).__new__(cls)
        return cls._instance
        
    def __init__(self, cosmology=None, snapRedshift=None, reinitialize=False):
        if not self._initialized or reinitialize:
            
            self.cosmology = cosmology if cosmology is not None else Planck15
            self.snapRedshift = snapRedshift if snapRedshift is not None else 0.
            
            self.a = 1 / (1 + self.snapRedshift)
            self.h = self.cosmology.h
            
            self.distance = u.kpc * self.a / self.h
            self.position = u.kpc * self.a / self.h
            self.density = (10**10 * u.Msun / self.h) / (u.kpc * self.a / self.h)**3
            self.mass = 10**10 * u.Msun / self.h
            self.sfr = u.Msun / u.yr
            self.velocity = u.km / u.s * np.sqrt(self.a)
            self.potential = (u.km / u.s)**2 / self.a
            self.temperature = u.K
            self.energy = (u.km / u.s)**2
            self.dimless = u.dimensionless_unscaled
            
            # mostly used: ckpc
            self.ckpc = u.kpc * self.a
            self.cMpc = u.Mpc * self.a
            self.pkpc = u.kpc
            self.pMpc = u.Mpc
            
            self.ckpc_over_h = u.kpc * self.a / self.h
            self.cMpc_over_h = u.Mpc * self.a / self.h
            self.pc_over_h = u.pc * self.a / self.h
            
            self._initialized = True
    
    def get_cosmology(self):
        return self.cosmology
    
    def get_snapRedshift(self):
        return self.snapRedshift
    
    def unit_convention(self):
        print(f'Current unit convention is in redshift {self.snapRedshift} for cosmology {self.cosmology.name}')
                
    def explain(self):
        """
        Default unit conventions are used by TNG simulation:
        
        - distance: ckpc / h
        - density: (10**10 * Msun) / (ckpc / h)^3
        - mass: 10**10 * Msun
        - sfr: Msun / yr
        - velocity: km * sqrt(a) / s
        - potential: (km / s)^2 / a
        - temperature: K
        - energy: (km / s)^2
        
        comoving quantities can be converted to physical ones by multiplying for the 
        appropriate power of the scale factor a. For instance, to convert a length 
        in physical units it is sufficient to multiply it by a, volumes need a factor a^3,
        densities a^-3 and so on. 
        
        When using other simulations, the unit convention can be re-defined. 
        """
        print(self.explain.__doc__)
        
def setup_logging(log_file="galaxygenius.log", log_level=logging.INFO, force_reconfigure=False):
    """Set up logging configuration that can be used across multiple modules
    
    Args:
        log_file: Path to the log file
        log_level: Logging level (default: logging.INFO)
        force_reconfigure: If True, forces reconfiguration even if handlers exist.
                          Use this after moving log files to update file handlers.
    """
    
    root_logger = logging.getLogger()
    
    # Check if we need to reconfigure (different log file or no handlers)
    needs_reconfigure = force_reconfigure
    if not needs_reconfigure:
        if not root_logger.handlers:
            needs_reconfigure = True
        elif log_file:
            # Check if any existing FileHandler points to a different file
            existing_log_files = [h.baseFilename for h in root_logger.handlers 
                                 if isinstance(h, logging.FileHandler)]
            if not existing_log_files or os.path.abspath(log_file) not in existing_log_files:
                needs_reconfigure = True
    
    if needs_reconfigure:
        # Remove all existing handlers to prevent duplicate logging
        for handler in root_logger.handlers[:]:
            handler.close()
            root_logger.removeHandler(handler)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        )
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)
        
        # File handler
        if log_file:
            # Ensure log directory exists
            log_dir = os.path.dirname(log_file)
            if log_dir and not os.path.exists(log_dir):
                os.makedirs(log_dir)
                
            file_handler = logging.FileHandler(log_file, mode='a')  # Append mode
            file_handler.setFormatter(formatter)
            root_logger.addHandler(file_handler)
        
        # Set log level
        root_logger.setLevel(log_level)
    
    # Return logger for the calling module
    caller_module = inspect.currentframe().f_back.f_globals['__name__']
    return logging.getLogger(caller_module)

def galaxygenius_data_dir():
    
    if os.environ.get('GALAXYGENIUS_DATA_DIR') is not None:
        dataDir = os.environ.get('GALAXYGENIUS_DATA_DIR').split(':')[0]
        if not os.path.exists(dataDir):
            raise FileNotFoundError('Data directory not found. Please set GALAXYGENIUS_DATA_DIR environment variable.')
    else:
        print('GALAXYGENIUS_DATA_DIR not set in environment variables. ' + 'Data directory falling to default path: ../Data')
        dataDir = '../Data'
        if not os.path.exists(dataDir):
            raise FileNotFoundError('Data directory not found. Please set GALAXYGENIUS_DATA_DIR environment variable.')
    
    return dataDir

def u2temp(u_energy: Union[float, u.Quantity], x_e: float) -> float:
    '''
    u_energy: InternelEnergy
    x_e: ElectronAbundance
    
    return:
    T: temperatrure in K
    '''
    
    if isinstance(u_energy, u.Quantity):
        u_energy = u_energy.value
    
    X_H = 0.76
    u_energy = u_energy * (u.km/u.s)**2
    mu = 4 / (1 + 3 * X_H + 4 * X_H * x_e) * const.m_p
    T = ((5/3 - 1) * u_energy/const.k_B * mu).to(u.K)
    T = T.value
    return T

def convert_to_rgb(bandpassImage: Union[np.ndarray, list], idx: list=[2, 3, 5]) -> np.ndarray:
    
    '''
    Convert the bandpass image to RGB image
    
    Args:
        bandpassImage: bandpass image
        idx: index of the bandpass image used to create the RGB image
        
    Returns:
        rgb: RGB image
    '''

    img_red = bandpassImage[idx[2]]
    img_green = bandpassImage[idx[1]]
    img_blue = bandpassImage[idx[0]]
    
    pctl = 99.5
    maxv = 0
    for img in [img_red, img_green, img_blue]:
        val = np.percentile(img, pctl)
        if val > maxv:
            maxv = val

    rgb = make_rgb(img_red, img_green, img_blue, interval=ManualInterval(vmin=0, vmax=maxv),
                   stretch=LogStretch(a=1000)) 
    
    return rgb
        
def split(string: str, castType: Union[type, NoneType]=None) -> list:
    
    '''
    Split the string into a list
    
    Args:
        string: string to be split
        castType: type of the elements in the list
        
    Returns:
        splits: list of the split string
    '''
    
    splits = [inc for inc in "".join(string.split()).split(',')]
    if castType is not None:
        splits = [castType(sp) for sp in splits]
        
    return splits

def get_wavelength_scale(filename: str) -> float:
    
    '''
    Get the wavelength scale of the filter
    
    Args:
        filename: filename of the filter
        
    Returns:
        wavelength_scale: wavelength scale of the filter
    '''
    
    with open(filename) as file:
        header = file.readline()
        
    if 'angstrom' in header or 'AA' in header:
        wavelength_scale = 1
    elif 'nm' in header:
        wavelength_scale = 10
    elif 'um' in header or 'micron' in header:
        wavelength_scale = 10 * 10**3
    else:
        # default consider as AA
        wavelength_scale = 1
    return wavelength_scale

def get_wavelength_unit(filename: str) -> u.Unit:
    
    with open(filename) as file:
        header = file.readline()
        
    if 'angstrom' in header or 'AA' in header:
        wave_unit = u.angstrom
    elif 'nm' in header:
        wave_unit = u.nm
    elif 'um' in header or 'micron' in header:
        wave_unit = u.um
    elif 'cm' in header:
        wave_unit = u.cm
    elif 'm' in header:
        wave_unit = u.m
    else:
        # default consider as AA
        wave_unit = u.angstrom
    return wave_unit
    

def calc_pivot(dataDir: str, survey: str, filter: str) -> float:
    
    '''
    Calculate the pivot wavelength of the filter
    
    Args:
        dataDir: directory of the data
        survey: survey name
        filter: filter name
        
    Returns:
        pivot: pivot wavelength of the filter
    '''
    
    filterDir = f'{dataDir}/filters/{survey}'
    filterLs = os.listdir(filterDir)
    filterNames = [name.split('.')[0] for name in filterLs]
    filename = filterLs[filterNames.index(filter)]
    filterName = os.path.join(filterDir, filename)
        
    wavelength_scale = get_wavelength_scale(filterName)
    
    try:
        transmission = np.loadtxt(filterName)
    except:
        transmission = np.load(filterName)
    
    transmission[:, 0] = transmission[:, 0] * wavelength_scale
    
    numerator = np.trapz(transmission[:, 1], transmission[:, 0])
    denomerator = np.trapz(transmission[:, 1] * transmission[:, 0]**-2,
                               transmission[:, 0])
    pivot = np.sqrt(numerator/denomerator)
    
    return pivot

def copyfile(src: str, tar: str):
    if os.path.exists(src) and not os.path.exists(tar):
        shutil.copyfile(src, tar)

def extend(values: Union[int, float, list, u.Quantity], nums: int) -> list:
    
    '''
    Extend the values to list with size consistent with nums
    
    Args:
        values: values to be extended
        nums: number of values to be extended
        
    Returns:
        values: extended values
    '''
    
    # if the values is a scaler
    if isinstance(values, (int, float, np.int32, np.int64, np.float32, np.float64)):
        values = u.Quantity([values] * nums)
    # elif the values is a list
    elif isinstance(values, list):
        if len(values) == nums:
            try:
                values = u.Quantity(values)
            except Exception as e:
                raise e
        else:
            raise ValueError(f'length of values is {len(values)}, which is not consistent with nums {nums}.')
    elif isinstance(values, u.Quantity):
        if values.size == 1:
            values = u.Quantity([values] * nums)
        elif values.size == nums:
            values = u.Quantity(values)
        else:
            raise ValueError(f'length of values is {len(values)}, which is not consistent with nums {nums}.')
    return values

def custom_serialier(obj):
    
    '''
    Custom serializer for json dump
    
    Args:
        obj: object to be serialized
        
    Returns:
        obj: serialized object
    '''
    
    if isinstance(obj, (np.float32, np.float64)):
        return float(obj)
    elif isinstance(obj, (np.int32, np.int64)):
        return int(obj)
    elif type(obj) == u.Quantity:
        if obj.size == 1:
            return {
                'value': obj.value.tolist(),
                'unit': str(obj.unit)
            }
        else:
            return {
                'value': obj.value.tolist(),
                'unit': str(obj.unit)
            }
    elif type(obj) == np.ndarray:
        return obj.tolist()
    
def _assign(value, unit):
    if u.Quantity(value).unit == unit:
        value = u.Quantity(value)
    else:
        value = u.Quantity(value) * unit
    return value
        
    
def check_exist_and_assign_unit(conf, key, unit):
    if key in conf:
        conf[key] = _assign(conf[key], unit)
        
    return conf
        
def assign_unit(conf):
    
    conf['minWavelength'] = _assign(conf['minWavelength'], u.um)
    conf['maxWavelength'] = _assign(conf['maxWavelength'], u.um)
    conf['maxBoxLength'] = _assign(conf['maxBoxLength'], u.kpc)
    conf['fieldOfView'] = _assign(conf['fieldOfView'], u.arcsec)
    conf['minStellarMass'] = _assign(conf['minStellarMass'], u.Msun)
    conf['maxStellarMass'] = _assign(conf['maxStellarMass'], u.Msun)
    conf['spatialResol'] = _assign(conf['spatialResol'], u.pc)
    conf['viewDistance'] = _assign(conf['viewDistance'], u.Mpc)
    conf['ageThreshold'] = _assign(conf['ageThreshold'], u.Myr)
    conf['PDRClearingTimescale'] = _assign(conf['PDRClearingTimescale'], u.Myr)
    conf['temperatureThreshold'] = _assign(conf['temperatureThreshold'], u.K)
    for survey in conf['surveys']:
        check_exist_and_assign_unit(conf, f'resolution_{survey}', u.pc)
        check_exist_and_assign_unit(conf, f'pixelScales_{survey}', u.arcsec)
        check_exist_and_assign_unit(conf, f'exposureTime_{survey}', u.s)
        check_exist_and_assign_unit(conf, f'aperture_{survey}', u.m)
        if f'skyBkg_{survey}' in conf:
            if isinstance(conf[f'skyBkg_{survey}'], list):
                conf[f'skyBkg_{survey}'] = _assign(conf[f'skyBkg_{survey}'], u.s**-1)
            else:
                match = re.match(r"\[([^\]]+)\]\s*(.+)", conf[f'skyBkg_{survey}'])
                if match:
                    values_str, unit_str = match.groups()
                    values = np.fromstring(values_str, sep=" ")
                    unit = u.Unit(unit_str)
                conf[f'skyBkg_{survey}'] = _assign(values, unit)
                
                
        check_exist_and_assign_unit(conf, f'darkCurrent_{survey}', u.s**-1)
        check_exist_and_assign_unit(conf, f'readOut_{survey}', u.dimensionless_unscaled)
        check_exist_and_assign_unit(conf, f'PSFFWHM_{survey}', u.arcsec)
        check_exist_and_assign_unit(conf, f'limitAperture_{survey}', u.arcsec)
    
    return conf

def to_quantity_list(value):
    match = re.match(r"\[([^\]]+)\]\s*(.+)", value)    
    if match:
        values_str, unit_str = match.groups()
    else:
        raise Exception("Input format not recognized")

    # Parse numbers (split by spaces)
    values = np.fromstring(values_str, sep=" ")
    # return a list
    return values.tolist()

def read_properties(workingDir: str) -> dict:
        with open(workingDir + '/properties.json', 'r') as file:
            properties = json.load(file)
            
        for key in properties.keys():
            if isinstance(properties[key], dict):
                subkeys = properties[key].keys()
                if 'value' in subkeys and 'unit' in subkeys:
                    properties[key] = u.Quantity(properties[key]['value']).astype(np.float32) * u.Unit(properties[key]['unit'])
                
                
        return properties
    
def read_config(directory: str) -> dict:
    with open(directory + '/config.json', 'r') as file:
        configs = json.load(file)
    
    for key in configs.keys():
            if isinstance(configs[key], dict):
                subkeys = configs[key].keys()
                if 'value' in subkeys and 'unit' in subkeys:
                    configs[key] = u.Quantity(configs[key]['value']).astype(np.float32) * u.Unit(configs[key]['unit'])
                
                
    return configs