import numpy as np
import h5py
import os
import illustris_python as ill
from astropy.cosmology import Planck15, WMAP9
import astropy.units as u
import astropy.constants as const
from .utils import *
import requests
import json
import tempfile
from scipy.interpolate import interp1d
from typing import Union

class Galaxy:
    def __init__(self, position: np.ndarray, smoothLength: np.ndarray, 
                 mass: np.ndarray, subhaloPos: np.ndarray, a: float, h: float):
        self.pos = position * a / h - subhaloPos
        self.smoothLength = smoothLength * a / h
        self.mass = mass * 10**10 / h

class PreProcess:
    
    def __init__(self, config: dict):
        self.config = config
        
        self.dataDir = self.config['dataDir']
        
        self.snapnum = np.int32(self.config['snapNum'])
        self.cosmology = Planck15
        self.snapz = self.config['snapRedshift']
        self.fage = self.__fage()
        
        # to avoid potential error in SKIRT execution
        if np.isclose(self.snapz, 0, atol=0.005):
            self.snapz = 0.005
            
        self.h = self.cosmology.h
        self.a = 1 / (1 + self.snapz)
        self.workingDir = self.config['workingDir']
        self.dataDir = self.config['dataDir']
        
        # TNG100-1 is the main high-resolution IllustrisTNG100 run including the full TNG physics model.
        self.simulation = self.config['simulation']
        
        if 'TNG' in self.config['simulation']:
            name = 'tng'
            self.cosmology = Planck15 # TNG uses Planck15 cosmology
        else:
            name = None
            print('Simulation name is unrecognized. Please manually input data by calling inputs().')
        
        if self.config['requests']:
            self.base_url = f"https://www.{name}-project.org/api/"
            self.headers = {'api-key': self.config['apiKey']}
            self.results = self.__get_subhalos_remote()
            
    def __get_snapz(self) -> float:
        
        if not self.config['requests']:
            
            snap = h5py.File(os.path.join(self.config['TNGPath'],
                                      f'snapdir_{self.snapnum:03d}/snap_{self.snapnum:03d}.0.hdf5'), 'r')
            snapz = dict(snap['Header'].attrs.items())['Redshift']
        else:
            url = f'{self.base_url}{self.simulation}/snapshots/{self.snapnum}'
            response = requests.get(url, headers=self.headers)
            data = response.json()
            snapz = data['redshift']
        
        return snapz
        

    def __read_subhalos(self) -> dict:
        snap_subhalos = ill.groupcat.loadSubhalos(self.config['TNGPath'], self.snapnum)
        return snap_subhalos
    
    
    def get_subhalos(self) -> dict:
        
        minStellarMass = np.float32(self.config['minStellarMass'])
        maxStellarMass = np.float32(self.config['maxStellarMass'])
        
        minStellarMass_in_10 = np.around(np.log10(minStellarMass), 2)
        maxStellarMass_in_10 = np.around(np.log10(maxStellarMass), 2)
        
        subhalos = {}
        
        if not self.config['requests']:
        
            minStellarMass = np.float32(self.config['minStellarMass'])
            maxStellarMass = np.float32(self.config['maxStellarMass'])
            
            snap_subhalos = self.__read_subhalos()
            stellarMass = snap_subhalos['SubhaloMassType'][:, 4] * 10**10 / self.h # Msun
            
            subhalo_indices = np.where((stellarMass > minStellarMass) \
                                        & (stellarMass < maxStellarMass) \
                                        & (snap_subhalos['SubhaloFlag'] == 1))[0]
            self.subhaloIDs = subhalo_indices
            self.stellarMasses = stellarMass[self.subhaloIDs] # Msun
            halfMassRad = snap_subhalos['SubhaloHalfmassRadType'][:, 4] * self.a / self.h # kpc
            self.halfMassRadii = halfMassRad[self.subhaloIDs]
            self.subhaloNum = self.subhaloIDs.shape[0]
            subhaloPos = snap_subhalos['SubhaloPos'] * self.a / self.h # in kpc
            self.subhaloPos = subhaloPos[self.subhaloIDs]
            self.subhaloSFR = snap_subhalos['SubhaloSFR'][self.subhaloIDs] # Msun/yr
            
            subhalos['subhaloNum'] = self.subhaloNum
            subhalos['subhaloIDs'] = self.subhaloIDs
            # subhalos['stellarMasses'] = self.stellarMasses
            subhalos['subhaloSFR'] = self.subhaloSFR
            
            subhalos['units'] = ['1', '1', 'Msun/yr']
            
            if maxStellarMass == np.inf:
                print_info = f'{self.subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass higher than 10^{minStellarMass_in_10:.2f} [M_sun]'
            else:
                print_info = f'{self.subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass from 10^{minStellarMass_in_10:.2f} to 10^{maxStellarMass_in_10:.2f} [M_sun]'
        else:
            
            subhaloNum = len(self.results)
            subhalosIDs = [result['id'] for result in self.results]
            # stellarMasses = [10 ** result['mass_log_stars'] for result in self.results]
            subhaloSFRs = [result['sfr'] for result in self.results]
            
            subhalos['subhaloNum'] = subhaloNum
            subhalos['subhaloIDs'] = subhalosIDs
            # subhalos['stellarMasses'] = stellarMasses
            subhalos['subhaloSFR'] = subhaloSFRs
            
            subhalos['units'] = ['1', '1', 'Msun', 'Msun/yr']
            
            if maxStellarMass == np.inf:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass higher than 10^{minStellarMass_in_10:.2f} [M_sun]'
            else:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass from 10^{minStellarMass_in_10:.2f} to 10^{maxStellarMass_in_10:.2f} [M_sun]'
        
        print(print_info)
        
        return subhalos
    
    def __get_subhalos_remote(self) -> list:
        
        cache_dir = os.path.join('.', 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        
        minStellarMass = np.float32(self.config['minStellarMass']) / 10**10 * self.h # in Msun / h
        if self.config['maxStellarMass'] == np.inf:
            maxStellarMass = np.inf
            url = f'{self.base_url}{self.simulation}/snapshots/{self.snapnum}' \
                + f'/subhalos/?mass_stars__gt={minStellarMass}&subhaloflag=1&limit=1000000'
            filename = f'subhalos_snap_{self.snapnum}_mass_gt_{minStellarMass}_1e10_Msun_over_h_subhaloflag1.json'
            cache_file = os.path.join(cache_dir, filename)
        else:
            maxStellarMass = np.float32(self.config['maxStellarMass']) / 10**10 * self.h
            url = f'{self.base_url}{self.simulation}/snapshots/{self.snapnum}' \
                + f'/subhalos/?mass_stars__gt={minStellarMass}&mass_stars' \
                + f'__lt={maxStellarMass}&subhaloflag=1&limit=1000000'
                
            filename = f'subhalos_snap_{self.snapnum}_mass_gt_{minStellarMass}' \
                + f'_lt_{maxStellarMass}_1e10_Msun_over_h_subhaloflag1.json'
            cache_file = os.path.join(cache_dir, filename)
        
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as f:
                results = json.load(f)
        else:
            response = requests.get(url, headers=self.headers)
            data = response.json()
            results = data.get('results', [])
            with open(cache_file, 'w') as f:
                json.dump(results, f, indent=4)
        
        return results

    
    def subhalo(self, subhaloID: int):
        
        if not self.config['requests']:
            idx = list(self.subhaloIDs).index(subhaloID)
            self.id = self.subhaloIDs[idx]
            self.mass = self.stellarMasses[idx]
            self.radius = self.halfMassRadii[idx]
            self.pos = self.subhaloPos[idx]

            mass_in_10 = np.log10(self.mass)
            print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [M_sun].')
        else:
            subhaloIDs = [result['id'] for result in self.results]
            idx = list(subhaloIDs).index(subhaloID)
            self.id = subhaloIDs[idx]
            subhalo = self.results[idx]
            subhalo_url = subhalo['url']
            print('Subhalo URL: ', subhalo_url)
            subhalo_response = requests.get(subhalo_url, headers=self.headers)
            self.data = subhalo_response.json()
            self.mass = self.data['mass_stars'] * 10**10 / self.h # Msun
            self.radius = self.data['halfmassrad_stars'] * self.a / self.h # kpc
            self.pos_x = self.data['pos_x'] * self.a / self.h # kpc
            self.pos_y = self.data['pos_y'] * self.a / self.h # kpc
            self.pos_z = self.data['pos_z'] * self.a / self.h # kpc
            self.pos = np.array([self.pos_x, self.pos_y, self.pos_z]) # kpc
            
            mass_in_10 = np.log10(self.mass)
            print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [M_sun].')
            
    def __fage(self) -> interp1d:
        z = np.linspace(0, 4, 1000)
        t = self.cosmology.age(z).to(u.Myr).value
        fage = interp1d(z, t, kind='cubic', fill_value='extrapolate')
        return fage
    
    def __calculate_angular_momentum_and_angles(self) -> tuple:
        
        
        print('------Calculating face-on and edge-on viewing angles------')
        
        stars = np.loadtxt(self.workingDir + '/stars.txt')
        
        positions = stars[:, :3] # positions
        velocities = stars[:, 7:10] # velocities
        masses = stars[:, 10] # current mass
        
        # Initialize total angular momentum vector
        total_angular_momentum = np.zeros(3)
        
        mask = np.where((np.abs(positions[:, 0]) < 30) \
                    & (np.abs(positions[:, 1]) < 30) \
                    & (np.abs(positions[:, 2]) < 30))[0]
        
        positions = positions[mask]
        velocities = velocities[mask]
        masses = masses[mask]

        # Calculate angular momentum for each particle
        for i in range(len(positions)):
            r = positions[i]
            v = velocities[i]
            m = masses[i]
            
            # Calculate angular momentum for the current particle
            angular_momentum = np.cross(r, m * v)
            
            # Add to total angular momentum
            total_angular_momentum += angular_momentum

        # Normalize the angular momentum vector to get the direction
        norm = np.linalg.norm(total_angular_momentum)
        angular_momentum_direction = total_angular_momentum / norm if norm != 0 else total_angular_momentum
        
        # Calculate inclination and azimuth angles
        x, y, z = angular_momentum_direction
        inclination = np.arccos(z / np.linalg.norm(angular_momentum_direction)) if norm != 0 else 0
        azimuth = np.arctan2(y, x)
        
        # Calculate perpendicular direction
        # Choose an arbitrary vector not parallel to angular_momentum_direction
        arbitrary_vector = np.array([1, 0, 0]) if x != 0 else np.array([0, 1, 0])
        perpendicular_direction = np.cross(angular_momentum_direction, arbitrary_vector)
        perpendicular_direction /= np.linalg.norm(perpendicular_direction)  # Normalize

        # Calculate perpendicular inclination and azimuth
        x_perp, y_perp, z_perp = perpendicular_direction
        perp_inclination = np.arccos(z_perp)  # Inclination with respect to z-axis
        perp_azimuth = np.arctan2(y_perp, x_perp)  # Azimuth in the xy-plane
        
        incs = [inclination, perp_inclination] # in radians
        azis = [azimuth, perp_azimuth] # in radians
        
        incs = [np.rad2deg(inc) for inc in incs]
        azis = [np.rad2deg(azi) for azi in azis]
        
        print('Face-on angle:', (incs[0], azis[0]))
        print('Edge-on angle:', (incs[1], azis[1]))
        
        return incs, azis
        
    def __get_value(self, key, castType=None):
        
        if hasattr(self, key):
            value = getattr(self, key)
        else:
            value = self.config[key]
        
        if castType is not None:
            value = castType(value)
        
        return value
        
        
    def __get_particles(self):
        
        print('Retrieving Stellar and Gas particles.')
        
        boxLengthScale = self.__get_value('boxLengthScale', np.float32)
        maxBoxLength = self.__get_value('maxBoxLength', np.float32)

        boxLength = self.radius * boxLengthScale
        boxLength = np.min([boxLength, maxBoxLength])
        self.boxLength = boxLength
        
        region = boxLength / 2.

        fields = ['GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime',
                'Coordinates', 'Velocities', 'Masses']
        
        if not self.config['requests']:
            starPart = ill.snapshot.loadSubhalo(self.config['TNGPath'], self.snapnum,
                                                self.id, 'star', fields)
        else:
            cutout_url = f"{self.base_url}{self.simulation}/snapshots/{self.snapnum}/subhalos/{self.id}/cutout.hdf5"
            
            params = {'stars': ','.join(fields)}
            response = requests.get(cutout_url, headers=self.headers, params=params)
            
            with tempfile.NamedTemporaryFile(suffix='.hdf5') as tmp:
                tmp.write(response.content)
                tmp.flush()
                
                with h5py.File(tmp.name, 'r') as f:
                    
                    starPart = {key: f[f'PartType4/{key}'][:] for key in fields}
            
            starPart['count'] = starPart['Coordinates'].shape[0]
            
        if 'TNG50' in self.simulation:
            stellar_smoothLength = 288
        elif 'TNG100' in self.simulation:
            stellar_smoothLength = 740
        elif 'TNG300' in self.simulation:
            stellar_smoothLength = 1480

        if self.snapz <= 1:
            # z = 0 ~ 1, smoothing length are changing
            stellar_smoothLength = stellar_smoothLength * 10**-3 / self.a * self.h # in ckpc/h, comoving scale
        else:
            # z > 1, smoothing length are constant as z = 1 value
            stellar_smoothLength = stellar_smoothLength * 10**-3 / 0.5 * self.h # in ckpc/h, comoving scale
        
        starPart['smoothLength'] = np.array([stellar_smoothLength] * starPart['count'])
        
        mask = np.where(starPart['GFM_StellarFormationTime'] > 0)[0]
        for key in starPart.keys():
            if key == 'count':
                pass
            else:
                starPart[key] = starPart[key][mask] # select stellar particles instead of wind
        
        snapshot_age = self.fage(self.snapz)
        starPart['age'] = snapshot_age - self.fage(1/starPart['GFM_StellarFormationTime'] - 1) # in Myr
        g = Galaxy(starPart['Coordinates'], starPart['smoothLength'], starPart['GFM_InitialMass'],
                   self.pos, a=self.a, h=self.h)

        ageThreshold = self.__get_value('ageThreshold', np.float32)
            

        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                    & (np.abs(g.pos[:, 1]) < region)\
                    & (np.abs(g.pos[:, 2]) < region)\
                    & (starPart['age'] < ageThreshold))[0]
        size = mask.shape[0]
        print('Starforming regions: ', size)

        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask]
        part['z'] = g.pos[:, 2][mask]
        part['smoothLength'] = g.smoothLength[mask] # smoothing length in kpc
        part['sfr'] = g.mass[mask] / (ageThreshold * 10**6) # constant SFR in Msun/yr
        part['Z'] = starPart['GFM_Metallicity'][mask]
        logCompactnessMean = self.__get_value('logCompactnessMean', np.float32)
        logCompactnessStd = self.__get_value('logCompactnessStd', np.float32)
        part['compactness'] = np.random.normal(loc=logCompactnessMean, 
                                               scale=logCompactnessStd, size=size) # from Kapoor et al. 2021
        logPressure = self.__get_value('logPressure', np.float32)
        pressure = (10**logPressure * const.k_B) * (u.J * u.cm**-3).to(u.J * u.m**-3) # J * m**3 == Pa
        pressure = pressure.value
        part['pressure'] = np.array([pressure] * size)
        age = starPart['age'][mask] # in Myr
        constantCoveringFactor = self.__get_value('constantCoveringFactor', bool)
        if constantCoveringFactor:
            coveringFactor = self.__get_value('coveringFactor', np.float32)
            part['covering'] = np.array([coveringFactor] * size)
        else:
            PDRClearingTimescale = self.__get_value('PDRClearingTimescale', np.float32)
            part['covering'] = np.array(np.exp(-age / PDRClearingTimescale)) # from Baes, M., et al. 2024
            
        part['velocities'] = starPart['Velocities'][mask] * np.sqrt(self.a) # in km/s
        part['masses'] = starPart['Masses'][mask] * 10**10 / self.h # Msun

        header = 'starforming regions\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n' \
                + 'Column 2: y-coordinate (kpc)\n' \
                + 'Column 3: z-coordinate (kpc)\n' \
                + 'Column 4: smoothing length (kpc)\n' \
                + 'Column 5: star formation rate (Msun/yr)\n' \
                + 'Column 6: metallicity (1)\n' \
                + 'Column 7: compactness (1)\n' \
                + 'Column 8: pressure (J/m3)\n' \
                + 'Column 9: coveringFactor (1)\n' \
                + 'Column 10: x-velocity (km/s)\n' \
                + 'Column 11: y-velocity (km/s)\n' \
                + 'Column 12: z-velocity (km/s)\n' \
                + 'Column 13: Mass (Msun)'

        info = np.column_stack((part['x'], part['y'], part['z'],
                            part['smoothLength'], part['sfr'],
                            part['Z'], part['compactness'],
                            part['pressure'], part['covering'],
                            part['velocities'], part['masses'])) # 13 params

        np.savetxt(self.workingDir + '/starforming_regions.txt', info, header=header)
        

        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                        & (np.abs(g.pos[:, 1]) < region)\
                        & (np.abs(g.pos[:, 2]) < region)\
                        & (starPart['age'] > ageThreshold))[0]
        size = mask.shape[0]
        print('Stars: ', size)
        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask] # in kpc
        part['z'] = g.pos[:, 2][mask] # in kpc
        part['smoothLength'] = g.smoothLength[mask] # smoothing length in kpc
        part['initialMass'] = g.mass[mask]
        part['Z'] = starPart['GFM_Metallicity'][mask]
        part['age'] = starPart['age'][mask]
        part['velocities'] = starPart['Velocities'][mask] * np.sqrt(self.a) # in km/s
        part['masses'] = starPart['Masses'][mask] * 10**10 / self.h # Msun
        
        header = 'Stars\n' \
                    + '\n' \
                    + 'Column 1: x-coordinate (kpc)\n'\
                    + 'Column 2: y-coordinate (kpc)\n'\
                    + 'Column 3: z-coordinate (kpc)\n'\
                    + 'Column 4: smoothing length (kpc)\n'\
                    + 'Column 5: initial mass (Msun)\n'\
                    + 'Column 6: metallicity (1)\n'\
                    + 'Column 7: age (Myr)\n' \
                    + 'Column 8: x-velocity (km/s)\n' \
                    + 'Column 9: y-velocity (km/s)\n' \
                    + 'Column 10: z-velocity (km/s)\n' \
                    + 'Column 11: Mass (Msun)'
        
        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['smoothLength'], part['initialMass'],
                                part['Z'], part['age'],
                                part['velocities'], part['masses']))
        np.savetxt(self.workingDir + '/stars.txt', info, header=header)
        
        part = {}
        includeDust = self.__get_value('includeDust', bool)
        if includeDust:
            try:
                fields = ['GFM_Metallicity', 'Coordinates', 'Masses',
                            'InternalEnergy', 'StarFormationRate', 'ElectronAbundance', 
                            'Density', 'Velocities']
                if not self.config['requests']:
                    gasPart = ill.snapshot.loadSubhalo(self.config['TNGPath'], 
                            self.snapnum, self.id, 'gas', fields)
                else:
                    cutout_url = f"{self.base_url}{self.simulation}/snapshots/{self.snapnum}/subhalos/{self.id}/cutout.hdf5"
                    
                    params = {'gas': ','.join(fields)}
                    response = requests.get(cutout_url, headers=self.headers, params=params)
                    
                    with tempfile.NamedTemporaryFile(suffix='.hdf5') as tmp:
                        tmp.write(response.content)
                        tmp.flush()
                        
                        with h5py.File(tmp.name, 'r') as f:
                            gasPart = {key: f[f'PartType0/{key}'][:] for key in fields}
                            
                    gasPart['count'] = gasPart['Coordinates'].shape[0]
                            
                volumes = gasPart['Masses'] / gasPart['Density']
                cellRadius = (3 * volumes / 4 / np.pi)**(1/3)
                gasPart['smoothLength'] = 2.5 * cellRadius # in ckpc/h
                
                gasPart['Temperature'] = u2temp(gasPart['InternalEnergy'], 
                                                gasPart['ElectronAbundance'])
                gasPart['velocities'] = gasPart['Velocities'] * np.sqrt(self.a) # in km/s
                
                gas = Galaxy(gasPart['Coordinates'], gasPart['smoothLength'], gasPart['Masses'], 
                            self.pos, a=self.a, h=self.h)
                
                spatialmask = np.where((np.abs(gas.pos[:, 0]) < region)\
                            & (np.abs(gas.pos[:, 1]) < region)\
                            & (np.abs(gas.pos[:, 2]) < region))[0]
                
                DISMModel = self.__get_value('DISMModel', str)
                if DISMModel == 'Camps_2016':
                    temperatureThreshold = self.__get_value('temperatureThreshold', np.float32)
                    othermask = np.where((gasPart['StarFormationRate'] > 0) \
                                    | (gasPart['Temperature'] < temperatureThreshold))[0]
                    mask = np.intersect1d(spatialmask, othermask)
                    size = mask.shape[0]
                elif DISMModel == 'Torrey_2012':
                    # Density from (10^10 Msun/h) / (ckpc/h)^3 to 10^10 h^2 Msun / kpc^3
                    density = gasPart['Density'] * self.a**-3
                    othermask = np.where(np.log10(gasPart['Temperature']) < (6 + 0.25 * np.log10(density)))[0]
                    mask = np.intersect1d(spatialmask, othermask)
                    size = mask.shape[0]
                                         
                print('Dusts: ', size)
                part['x'] = gas.pos[:, 0][mask] # in kpc
                part['y'] = gas.pos[:, 1][mask]
                part['z'] = gas.pos[:, 2][mask]
                part['smoothLength'] = gas.smoothLength[mask]
                part['mass'] = gas.mass[mask]
                part['Z'] = gasPart['GFM_Metallicity'][mask]
                part['Temperature'] = gasPart['Temperature'][mask]
                part['velocities'] = gasPart['Velocities'][mask]
            except:
                part['x'] = np.array([])
                part['y'] = np.array([])
                part['z'] = np.array([])
                part['smoothLength'] = np.array([])
                part['mass'] = np.array([])
                part['Z'] = np.array([])
                part['Temperature'] = np.array([])
                part['velocities'] = np.array([])
        
        else:
            part['x'] = np.array([])
            part['y'] = np.array([])
            part['z'] = np.array([])
            part['smoothLength'] = np.array([])
            part['mass'] = np.array([])
            part['Z'] = np.array([])
            part['Temperature'] = np.array([])
            part['velocities'] = np.array([])
        header = 'dusts\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n' \
                + 'Column 2: y-coordinate (kpc)\n' \
                + 'Column 3: z-coordinate (kpc)\n' \
                + 'Column 4: smoothing length (kpc)\n' \
                + 'Column 5: mass (Msun)\n' \
                + 'Column 6: metallicity (1)\n' \
                + 'Column 7: temperature (K)\n' \
                + 'Column 8: x-velocity (km/s)\n' \
                + 'Column 9: y-velocity (km/s)\n' \
                + 'Column 10: z-velocity (km/s)'
                
        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['smoothLength'], part['mass'],
                                part['Z'], part['Temperature'],
                                part['velocities']))

        np.savetxt(self.workingDir + '/dusts.txt', info, header=header)
        
    def __get_properties_survey(self, distance: float, properties: dict) -> dict:
        
        postProcessing = self.__get_value('postProcessing', bool)
        if postProcessing:
            surveys = self.__get_value('surveys', list)
            spatialResols = []
            for survey in surveys:
                filters = self.__get_value(f'filters_{survey}', list)
                numfilters = len(filters)
                if self.__get_value(f'resolFromPix_{survey}', bool):
                    pixelScales = self.__get_value(f'pixelScales_{survey}')
                    pixelScales = extend(pixelScales, numfilters)
                    resol = [distance * 10**6 * np.deg2rad(ang / 3600) for ang in pixelScales]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                else:
                    resol = self.__get_value(f'resolution_{survey}')
                    resol = extend(resol, numfilters)
                    pixelScales = [np.rad2deg(res / (distance * 10**6)) * 3600 for res in resol]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                
                properties[f'resolution_{survey}'] = resol
                properties[f'angleRes_{survey}'] = pixelScales
                properties[f'ps_in_sr_{survey}'] = ps_in_sr

                spatialResols += resol

                numExposure = self.__get_value(f'numExposure_{survey}')
                numExposure = extend(numExposure, numfilters)
                properties[f'numExposure_{survey}'] = numExposure
                
                exposureTime = self.__get_value(f'exposureTime_{survey}')
                exposureTime = extend(exposureTime, numfilters)
                properties[f'exposureTime_{survey}'] = exposureTime
                
                aperture = self.__get_value(f'aperture_{survey}', np.float32)
                properties[f'aperture_{survey}'] = aperture
                
                properties[f'numfilters_{survey}'] = numfilters
                properties[f'filters_{survey}'] = filters
            
            spatialResol = np.min(spatialResols)
            
            for survey in surveys:
                ratios = [spatialResol / res for res in properties[f'resolution_{survey}']]
                properties[f'ratios_{survey}'] = ratios
            
        else:
            spatialResol = self.__get_value('spatialResol', np.float32)
        
        properties['baseRes'] = spatialResol
        
        return properties
    
    def __get_properties(self) -> dict:
        
        properties = {}
        properties['subhaloID'] = self.id
        properties['stellarMass'] = np.log10(self.mass)
        properties['redshift'] = self.snapz
        distance = self.cosmology.luminosity_distance(self.snapz).value
        properties['lumiDis'] = distance
        
        properties = self.__get_properties_survey(distance, properties)
        
        faceAndEdge = self.__get_value('faceAndEdge', bool)
        
        if hasattr(self, 'inclinations') and hasattr(self, 'azimuths'):
            
            inclinations = self.inclinations
            azimuths = self.azimuths
            numViews = len(self.inclinations)
            self.numViews = numViews
            
        if faceAndEdge:
            inclinations, azimuths = self.__calculate_angular_momentum_and_angles()
            numViews = 2
            self.numViews = numViews
            
        if not faceAndEdge and not hasattr(self, 'inclinations') and not hasattr(self, 'azimuths'):
            
            numViews = self.__get_value('numViews', np.int32)
            randomViews = self.__get_value('randomViews', bool)
            
            if randomViews:
                inclinations = np.random.uniform(0, 180, numViews).tolist()
                azimuths = np.random.uniform(-360, 360, numViews).tolist()
            else:
                inclinations = self.config['inclinations']
                azimuths = self.config['azimuths']
            
        properties['numViews'] = numViews
        properties['inclinations'] = inclinations
        properties['azimuths'] = azimuths
        fov = self.boxLength * 10**3
        properties['FoV'] = fov
        
        return properties
    
    def __create_ski(self):
        
        print('Creating .ski file.')
        
        mode = self.__get_value('simulationMode', str)
        ski_file = os.path.join(self.dataDir,  f'ski_templates/{mode}_template.ski')
        
        with open(ski_file, 'r') as file:
            data = file.read()
        
        data = data.replace('redshift="0.008"', f'redshift="{self.snapz}"')
        
        self.properties = self.__get_properties()
        
        SEDFamily = {'BC03': 'BruzualCharlotSEDFamily',
                     'FSPS': 'FSPSSEDFamily'}
        SEDFamily = SEDFamily[self.config['SEDFamily']]
        
        if SEDFamily == 'FSPSSEDFamily':
            data = data.replace('resolution="High"', '')
        
        if mode == 'DustEmission':
            dustEmissionType = self.config['dustEmissionType']
            data = data.replace('dustEmissionType="Equilibrium"',
                                f'dustEmissionType="{dustEmissionType}"')
        
        initialMassFunction = self.config['initialMassFunction']
        
        data = data.replace('BruzualCharlotSEDFamily', SEDFamily)
        data = data.replace('Chabrier', initialMassFunction)
        
        numPackets = self.__get_value('numPackets', np.float32)
        data = data.replace('numPackets="1e7"', f'numPackets="{numPackets}"')
        
        minWavelength = self.__get_value('minWavelength', np.float32) * (1 + self.snapz)
        maxWavelength = self.__get_value('maxWavelength', np.float32) * (1 + self.snapz)
        
        data = data.replace('minWavelength="0.01 micron"', f'minWavelength="{minWavelength} micron"')
        data = data.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength} micron"')
        
        self.properties['minWavelength'] = minWavelength # micron
        self.properties['maxWavelength'] = maxWavelength # micron    

        massFraction = self.__get_value('massFraction', np.float32)
        data = data.replace('massFraction="0.3"', f'massFraction="{massFraction}"')
        dustConfig = '<ZubkoDustMix numSilicateSizes="15" numGraphiteSizes="15" numPAHSizes="15"/>'
        dustModel = self.__get_value('dustModel', str)
        numSilicateSizes = self.__get_value('numSilicateSizes', np.int32)
        numGraphiteSizes = self.__get_value('numGraphiteSizes', np.int32)
        numPAHSizes = self.__get_value('numPAHSizes', np.int32)
        numHydrocarbonSizes = self.__get_value('numHydrocarbonSizes', np.int32)
        if dustModel == 'ThemisDustMix':
            dustConfigNew = f'<{dustModel} numHydrocarbonSizes="{numHydrocarbonSizes}" numSilicateSizes="{numSilicateSizes}"/>'
        else:
            dustConfigNew = f'<{dustModel} numSilicateSizes="{numSilicateSizes}" numGraphiteSizes="{numGraphiteSizes}" numPAHSizes="{numPAHSizes}"/>'
        data = data.replace(dustConfig, dustConfigNew)

        minLevel = self.__get_value('minLevel', np.int32)
        maxLevel = self.__get_value('maxLevel', np.int32)
        data = data.replace('minLevel="8"', f'minLevel="{minLevel}"')
        data = data.replace('maxLevel="12"', f'maxLevel="{maxLevel}"')
        
        spatialRange = self.boxLength / 2. * 10**3 # in pc
        data = data.replace('minX="-5e4 pc"', f'minX="{-spatialRange} pc"')
        data = data.replace('maxX="5e4 pc"', f'maxX="{spatialRange} pc"')
        data = data.replace('minY="-5e4 pc"', f'minY="{-spatialRange} pc"')
        data = data.replace('maxY="5e4 pc"', f'maxY="{spatialRange} pc"')
        data = data.replace('minZ="-5e4 pc"', f'minZ="{-spatialRange} pc"')
        data = data.replace('maxZ="5e4 pc"', f'maxZ="{spatialRange} pc"')
        
        wavelengthGrid = self.__get_value('wavelengthGrid', str)
        grid_type = {'Linear': 'LinWavelengthGrid',
                'Log': 'LogWavelengthGrid'}
        grid_type = grid_type[wavelengthGrid]
        data = data.replace('LinWavelengthGrid', grid_type)
        
        numWavelengths = self.__get_value('numWavelengths', np.int32)
        data = data.replace('numWavelengths="1000"', f'numWavelengths="{numWavelengths}"')

        begin_str = '<FullInstrument'
        end_str = '</FullInstrument>'
        offset = len(end_str)
        
        idx_begin = data.index(begin_str)
        idx_end = data.index(end_str) + offset
        
        instrumentInfo = data[idx_begin: idx_end]
        remainingInfo = data[idx_end:]
        
        data = data.replace(instrumentInfo, '')
        
        instrumentInfo = '\n' + instrumentInfo + '\n'
        fov = self.boxLength * 10**3 # in pc

        spatialResol = self.properties['baseRes']
        inclinations = self.properties['inclinations']
        azimuths = self.properties['azimuths']
        
        numViews = self.properties['numViews']

        numPixels = int(fov / spatialResol)

        recordComponents = 'true' if self.__get_value('recordComponents', bool) else 'false'
        
        insert_begin_idx = idx_begin
        for i, (inclination, azimuth) in enumerate(zip(inclinations, azimuths)):
            info = instrumentInfo.replace('view', f'view_{i:02d}')
            info = info.replace('inclination="0 deg"', f'inclination="{inclination} deg"')
            info = info.replace('azimuth="0 deg"', f'azimuth="{azimuth} deg"')
            info = info.replace('fieldOfViewX="1e5 pc"', f'fieldOfViewX="{fov} pc"')
            info = info.replace('fieldOfViewY="1e5 pc"', f'fieldOfViewY="{fov} pc"')
            info = info.replace('numPixelsX="1000"', f'numPixelsX="{numPixels}"')
            info = info.replace('numPixelsY="1000"', f'numPixelsY="{numPixels}"')
            info = info.replace('recordComponents="false"', f'recordComponents="{recordComponents}"')
            
            data = data[:insert_begin_idx] + info
            insert_begin_idx = insert_begin_idx + len(info)

        data = data + remainingInfo
        with open(self.workingDir + '/skirt.ski', 'w') as file:
            file.write(data)
        
        with open(self.workingDir + '/properties.json', 'w') as file:
            json.dump(self.properties, file, 
                      default=custom_serialier, indent=4)

        print('------estimate memory usage------')
        print(f'numViews: {numViews}')
        print(f'numSpatialPixels: {numPixels}')
        print(f'numWavelengthPixels: {numWavelengths}')
        factor = 7 if recordComponents else 1
        numPixels = np.float64(numPixels) # avoid overflow
        dataCubeSize = np.int64(numPixels ** 2 * numWavelengths * numViews)
        dataSize_in_GB = np.around(dataCubeSize * 8 * factor * 1e-9, 3)
        print(f'Estimated memory usage: {dataSize_in_GB} GB')
        
    def __get_attributes(self, directory):
        
        attributes = vars(self)
        attr_keys = list(attributes.keys())
        
        config = attributes['config']
        attr_keys.remove('config')
        
        for key in attr_keys:
            if key in config.keys():
                config[key] = attributes[key]

        with open(os.path.join(directory, 'config.json'), 'w') as file:
            json.dump(config, file, default=custom_serialier, indent=4)
    
    def prepare(self, data: Union[dict, NoneType]=None):
        
        if data is not None:
            exist_keys = self.config.keys()
            for key, value in data.items():
                if key in exist_keys:
                    setattr(self, key, value)
            
        os.makedirs(self.workingDir, exist_ok=True)
        self.__get_particles()
        self.__create_ski()
        self.__get_attributes(self.workingDir)
        
    def inputs(self, data: dict):
        '''
        Inputs for custom hydrodynamical simulations
        
        data: dict
        data['snapRedshift']: float
        data['cosmology']: astropy.cosmology.Cosmology
        data['stellarMass']: float
        data['subhaloID']: int
        data['boxLength']: float
        
        '''
        
        keys = data.keys()
        
        self.snapz = data['snapRedshift'] # include
        self.a = 1 / (1 + self.snapz)
        self.cosmology = data['cosmology'] # include
        self.mass = data['stellarMass'] # include 
        self.id = data['subhaloID'] # include
        self.boxLength = data['boxLength'] # include
        
        keys.remove('snapRedshift')
        keys.remove('cosmology')
        keys.remove('stellarMass')
        keys.remove('subhaloID')
        keys.remove('boxLength')
        
        # reset attributes
        for key in keys:
            setattr(self, key, data[key])
        
        mass_in_10 = np.log10(self.mass)
        
        print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [Msun]')
        
        os.makedirs(self.workingDir, exist_ok=True)
        self.__create_ski()
        self.__get_attributes(self.workingDir)
