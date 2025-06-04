import numpy as np
import astropy.units as u
import astropy.constants as const
import sys
import os
import shutil
from astropy.visualization import ManualInterval, LogStretch, make_rgb
from typing import Union
from types import NoneType

def u2temp(u_energy: float, x_e: float) -> float:
    '''
    u_energy: InternelEnergy
    x_e: ElectronAbundance
    
    return:
    T: temperatrure in K
    '''
    
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
        print('unrecognized unit in wavelength!')
        sys.exit(0)
    return wavelength_scale
        

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

def extend(values: Union[int, float, list], nums: int) -> list:
    
    '''
    Extend the values to list with size consistent with nums
    
    Args:
        values: values to be extended
        nums: number of values to be extended
        
    Returns:
        values: extended values
    '''
    
    if isinstance(values, (int, float, np.int32, np.int64, np.float32, np.float64)):
        values = [values] * nums
    else:
        values = values
    
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
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, object):
        return str(obj)
    else:
        return obj