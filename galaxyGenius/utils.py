import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck15
from PIL import Image
import sys
import os
import shutil
from astropy.visualization import ManualInterval, LogStretch, make_rgb

def u2temp(u_energy, x_e):
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

def fage(z):
    
    '''
    z: redshift
    
    return:
    time: age in Myr
    '''
    
    time = Planck15.age(z).to(u.Myr)
    time = time.value
    return time

def convert_to_rgb(bandpassImage, idx=[2, 3, 5]):
    
    # RGB -- long -> short wavelength

    bandpassImage = np.array(bandpassImage)

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
        
def split(string, castType=None):
    splits = [inc for inc in "".join(string.split()).split(',')]
    if castType is not None:
        splits = [castType(sp) for sp in splits]
        
    return splits

def get_wavelength_scale(filename):
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
        

def calc_pivot(survey, filter):
    
    filterDir = f'../Data/filters/{survey}'
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

def copyfile(src, tar):
    if os.path.exists(src) and not os.path.exists(tar):
        shutil.copyfile(src, tar)

def extend(values, nums):
    if len(values) == 1:
        values = [values][0] * nums
    else:
        values
    
    return values

