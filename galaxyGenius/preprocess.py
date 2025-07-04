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
import sys
import time

class PartSmoothing():
    def __init__(self, position: np.ndarray, smoothLength: np.ndarray, 
                 mass: np.ndarray, subhaloPos: np.ndarray, a: float, h: float):
        
        '''
        Store particles with smoothing lengths, convert to physical units without hubble parameter
        
        Args:
            position: positions of particles
            smoothLength: smoothing lengths of particles
            mass: masses of particles
            subhaloPos: positions of subhalos
            a: scale factor
            h: hubble parameter
        '''
        
        self.pos = position * a / h - subhaloPos
        self.smoothLength = smoothLength * a / h
        self.mass = mass * 10**10 / h
        
class PartNoSmoothing():
    def __init__(self, position: np.ndarray, mass: np.ndarray, subhaloPos: np.ndarray, a: float, h: float):
        
        '''
        Store particles without smoothing lengths, convert to physical units without hubble parameter
        
        Args:
            position: positions of particles
            mass: masses of particles
            subhaloPos: positions of subhalos
            a: scale factor
            h: hubble parameter
        '''
        
        self.pos = position * a / h - subhaloPos
        self.mass = mass * 10**10 / h
        
class PreProcess:
    
    def __init__(self, config: dict):
        
        self.config = config
        self.cosmology = Planck15
        self.fage = self.__fage()
        self.__init()
        
    def __fage(self) -> interp1d:
        z = np.linspace(0, 4, 1000)
        t = self.cosmology.age(z).to(u.Myr).value
        fage = interp1d(z, t, kind='cubic', fill_value='extrapolate')
        return fage
        
    def __init(self):
        
        self.snapRedshift = self.config['snapRedshift']
        self.viewRedshift = self.config['viewRedshift']
        
        # to avoid potential error in SKIRT execution
        if np.isclose(self.viewRedshift, 0, atol=0.005):
            self.viewRedshift = 0.005
        
        self.h = self.cosmology.h
        self.a = 1 / (1 + self.snapRedshift)
        self.workingDir = self.config['workingDir']
        os.makedirs(self.workingDir, exist_ok=True)
        
        self.dataDir = self.config['dataDir']        
        self.simulation = self.config['simulation']
        self.snapnum = self.config['snapNum']
        
        if 'TNG' in self.config['simulation']:
            name = 'tng'
        else:
            name = None
            print('Simulation name is unrecognized. Please manually input data by calling inputs().')
            
        if self.config['requests']:
            self.base_url = f"https://www.{name}-project.org/api/"
            self.headers = {'api-key': self.config['apiKey']}
            
    def __make_request_with_retry(self, url: str, headers: dict = None, 
                                params: dict = None, max_retries: int = 5) -> requests.Response:
        """
        Make an HTTP request with retry logic.
        
        Args:
            url: The URL to request
            headers: Optional headers for the request
            params: Optional parameters for the request
            max_retries: Maximum number of retry attempts
        
        Returns:
            Response object if successful
            
        Raises:
            SystemExit if all retries fail
        """
        
        start_time = time.time()
        content = None
        
        for attempt in range(max_retries):
            try:
                response = requests.get(url, headers=headers, params=params)
                response.raise_for_status()  # Raises an HTTPError for bad responses
                content = response
            except requests.exceptions.RequestException as e:
                if attempt == max_retries - 1:  # Last attempt
                    print(f"Error: Failed to make request after {max_retries} attempts.")
                    print(f"URL: {url}")
                    print(f"Error message: {str(e)}")
                else:
                    print(f"Request failed (attempt {attempt + 1}/{max_retries}). Retrying...")
                    time.sleep(2 ** attempt)  # Exponential backoff        
                    
        end_time = time.time()
        print(f"Requesting time taken: {end_time - start_time:.2f} seconds")
        
        if content is None:
            sys.exit(1)
        
        return content
                    
    def __read_subhalos(self) -> dict:
        snap_subhalos = ill.groupcat.loadSubhalos(self.config['TNGPath'], self.config['snapNum'])
        return snap_subhalos
    
    def __retrieve_subhalos_with_requests(self) -> list:
        
        cache_dir = os.path.join('.', 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        
        minStellarMass = np.float32(self.config['minStellarMass']) / 10**10 * self.h
        maxStellarMass = np.float32(self.config['maxStellarMass']) / 10**10 * self.h
        snapnum = self.config['snapNum']
        
        if maxStellarMass == np.inf:
            url = f'{self.base_url}{self.simulation}/snapshots/{snapnum}' \
                + f'/subhalos/?mass_stars__gt={minStellarMass}&subhaloflag=1&limit=1000000'
            
            filename = f'subhalos_snap_{snapnum}_mass_gt_{minStellarMass}_1e10_Msun_over_h_subhaloflag1.json'
        else:
            url = f'{self.base_url}{self.simulation}/snapshots/{snapnum}' \
                + f'/subhalos/?mass_stars__gt={minStellarMass}&mass_stars' \
                + f'__lt={maxStellarMass}&subhaloflag=1&limit=1000000'
                
            filename = f'subhalos_snap_{snapnum}_mass_gt_{minStellarMass}' \
                + f'_lt_{maxStellarMass}_1e10_Msun_over_h_subhaloflag1.json'
        
        cache_file = os.path.join(cache_dir, filename)
        
        response = self.__make_request_with_retry(url, headers=self.headers)
        data = response.json()
        results = data.get('results', [])
        
        with open(cache_file, 'w') as f:
            json.dump(results, f, indent=4)
        
        return results
    
    def get_subhalos(self) -> dict:
        
        '''
        Get subhalos
        
        Returns:
            subhalos: dictionary containing number of subhalos, subhaloIDs, subhaloSFRs
        '''
        
        minStellarMass = np.float32(self.config['minStellarMass'])
        maxStellarMass = np.float32(self.config['maxStellarMass'])
        
        minStellarMass_in_10 = np.around(np.log10(minStellarMass), 2)
        maxStellarMass_in_10 = np.around(np.log10(maxStellarMass), 2)
        
        subhalos = {}
        
        if not self.config['requests']:
            
            snap_subhalos = self.__read_subhalos()
            stellarMass = snap_subhalos['SubhaloMassType'][:, 4] * 10**10 / self.h # Msun
            
            subhalo_indices = np.where((stellarMass > minStellarMass) \
                                        & (stellarMass < maxStellarMass) \
                                        & (snap_subhalos['SubhaloFlag'] == 1))[0]
            
            self.subhaloIDs = subhalo_indices # indices are same as subhaloIDs
            self.stellarMasses = stellarMass[self.subhaloIDs] # in Msun
            halfMassRad = snap_subhalos['SubhaloHalfmassRadType'][:, 4] * self.a / self.h # in kpc
            self.halfMassRadii = halfMassRad[self.subhaloIDs]
            subhaloNum = self.subhaloIDs.shape[0]
            subhaloPos = snap_subhalos['SubhaloPos'] * self.a / self.h # in kpc
            self.subhaloPos = subhaloPos[self.subhaloIDs]
            self.subhaloSFR = snap_subhalos['SubhaloSFR'][self.subhaloIDs] # in Msun/yr
            
            subhalos['subhaloNum'] = subhaloNum
            subhalos['subhaloIDs'] = self.subhaloIDs
            subhalos['subhaloSFR'] = self.subhaloSFR
            
            subhalos['units'] = ['1', '1', 'Msun/yr']
            
            if maxStellarMass == np.inf:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass higher than 10^{minStellarMass_in_10:.2f} [M_sun]'
            else:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass from 10^{minStellarMass_in_10:.2f} to 10^{maxStellarMass_in_10:.2f} [M_sun]'
            
        else:
            
            results = self.__retrieve_subhalos_with_requests()
            self.results = results # to access in following functions
            
            subhaloNum = len(results)
            self.subhaloIDs = [result['id'] for result in results]
            self.subhaloSFR = [result['sfr'] for result in results]
            
            subhalos['subhaloNum'] = subhaloNum
            subhalos['subhaloIDs'] = self.subhaloIDs
            subhalos['subhaloSFR'] = self.subhaloSFR
            
            subhalos['units'] = ['1', '1', 'Msun/yr']
            
            if maxStellarMass == np.inf:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass higher than 10^{minStellarMass_in_10:.2f} [M_sun]'
            else:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass from 10^{minStellarMass_in_10:.2f} to 10^{maxStellarMass_in_10:.2f} [M_sun]'
        
        print(print_info)
        
        return subhalos
    
    def subhalo(self, subhaloID: int):
        
        '''
        Specify one subhalo to be processed
        
        Args:
            subhaloID: subhaloID of the subhalo to be processed
        '''
        
        if not self.config['requests']:
            idx = idx = list(self.subhaloIDs).index(subhaloID)
            self.id = self.subhaloIDs[idx]
            self.mass = self.stellarMasses[idx]
            self.radius = self.halfMassRadii[idx]
            self.pos = self.subhaloPos[idx]
            self.sfr = self.subhaloSFR[idx]

            mass_in_10 = np.log10(self.mass)
            
            print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [M_sun].')
            print(f'SFR of Subhalo {self.id} is {self.sfr:.2f} [M_sun/yr].')
            
        else:
            subhaloIDs = [result['id'] for result in self.results]
            idx = list(subhaloIDs).index(subhaloID)
            self.id = subhaloIDs[idx]
            subhalo = self.results[idx]
            subhalo_url = subhalo['url']
            print('Subhalo URL: ', subhalo_url)
            subhalo_response = self.__make_request_with_retry(subhalo_url, headers=self.headers)
            self.data = subhalo_response.json()
            self.mass = self.data['mass_stars'] * 10**10 / self.h # Msun
            self.radius = self.data['halfmassrad_stars'] * self.a / self.h # kpc
            self.pos_x = self.data['pos_x'] * self.a / self.h # kpc
            self.pos_y = self.data['pos_y'] * self.a / self.h # kpc
            self.pos_z = self.data['pos_z'] * self.a / self.h # kpc
            self.pos = np.array([self.pos_x, self.pos_y, self.pos_z]) # kpc
            self.sfr = self.subhaloSFR[idx]
            
            mass_in_10 = np.log10(self.mass)
            
            print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [M_sun].')
            print(f'SFR of Subhalo {self.id} is {self.sfr:.2f} [M_sun/yr].')
            
            
    def __calculate_angular_momentum_and_angles(self) -> tuple:
        
        print('------Calculating face-on and edge-on viewing angles------')
        
        stars = np.loadtxt(self.workingDir + '/stars.txt')
        
        headers = []
        with open(self.workingDir + '/stars.txt', 'r') as f:
            for line in f:
                if line.startswith('# Column'):
                    headers.append(line)
                
        idx_coords = []
        idx_vels = []
        idx_masses = []
        
        for i, header in enumerate(headers):
            if 'coordinate' in header.lower():
                idx_coords.append(i)
            elif 'velocity' in header.lower():
                idx_vels.append(i)
            elif 'Mass' in header:
                idx_masses.append(i)
                
        positions = stars[:, idx_coords] # positions
        velocities = stars[:, idx_vels] # velocities
        masses = stars[:, idx_masses] # current mass
        
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
    
    def __read_temp_file(self, response_content, partType='star', suffix='.hdf5', fields=None):
        """
        Handle temporary file operations with retry logic.
        
        Args:
            response_content: Content to write to temporary file
            suffix: File suffix
            fields: Fields to read
        
        Returns:
            Dictionary containing the particle data
        """
        try:
            with tempfile.NamedTemporaryFile(suffix=suffix, dir=self.workingDir) as tmp:
                tmp.write(response_content)
                tmp.flush()
                
                with h5py.File(tmp.name, 'r') as f:
                    
                    if partType == 'star':
                        particle_data = {key: f[f'PartType4/{key}'][:] for key in fields}
                    elif partType == 'gas':
                        particle_data = {key: f[f'PartType0/{key}'][:] for key in fields}
                    
            particle_data['count'] = len(particle_data['Coordinates'])
                    
        except Exception as e:
            print(f'Error reading temporary file: {e}')
            particle_data = {}
            particle_data['count'] = 0
            
        return particle_data
    
    def __get_particles(self):
        
        print('Retrieving Stellar and Gas particles.')
        
        boxLengthScale = np.float32(self.config['boxLengthScale'])
        maxBoxLength = np.float32(self.config['maxBoxLength'])
        
        boxLength = self.radius * boxLengthScale
        boxLength = np.min([boxLength, maxBoxLength])
        self.boxLength = boxLength # in kpc
        
        region = boxLength / 2.
        
        fields = ['GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime',
                'Coordinates', 'Velocities', 'Masses']
        
        if not self.config['requests']:
            starPart = ill.snapshot.loadSubhalo(self.config['TNGPath'], self.config['snapNum'],
                                                self.id, 'star', fields)
        else:
            cutout_url = f"{self.base_url}{self.simulation}/snapshots/{self.config['snapNum']}/subhalos/{self.id}/cutout.hdf5"
            
            params = {'stars': ','.join(fields)}
            response = self.__make_request_with_retry(cutout_url, headers=self.headers, params=params)
            
            starPart = self.__read_temp_file(response.content, partType='star', fields=fields)
            # starPart['count'] = len(starPart['Coordinates'])
            
        if 'TNG50' in self.simulation:
            stellar_smoothLength = 288
        elif 'TNG100' in self.simulation:
            stellar_smoothLength = 740
        elif 'TNG300' in self.simulation:
            stellar_smoothLength = 1480
            
        if self.snapRedshift <= 1:
            stellar_smoothLength = stellar_smoothLength * 10**-3 / self.a * self.h # in ckpc/h, comoving scale
        else:
            stellar_smoothLength = stellar_smoothLength * 10**-3 / 0.5 * self.h # in ckpc/h, comoving scale
        
        starPart['smoothLength'] = np.array([stellar_smoothLength] * starPart['count'])
        
        mask = np.where(starPart['GFM_StellarFormationTime'] > 0)[0] # select stellar particles instead of wind
        for key in starPart.keys():
            if key == 'count':
                pass
            else:
                starPart[key] = starPart[key][mask] 
            
        snapshot_age = self.fage(self.snapRedshift)
        starPart['age'] = snapshot_age - self.fage(1/starPart['GFM_StellarFormationTime'] - 1) # in Myr
        g = PartSmoothing(starPart['Coordinates'], starPart['smoothLength'], 
                   starPart['GFM_InitialMass'], self.pos, a=self.a, h=self.h)
        
        ageThreshold = np.float32(self.config['ageThreshold'])
        
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
        part['smoothLength'] = g.smoothLength[mask]
        part['sfr'] = g.mass[mask] / (ageThreshold * 10**6) # constant SFR in Msun/yr
        part['Z'] = starPart['GFM_Metallicity'][mask]
        logCompactnessMean = np.float32(self.config['logCompactnessMean'])
        logCompactnessStd = np.float32(self.config['logCompactnessStd'])
        part['compactness'] = np.random.normal(loc=logCompactnessMean, 
                                               scale=logCompactnessStd, size=size) # from Kapoor et al. 2021
        logPressure = np.float32(self.config['logPressure'])
        pressure = (10**logPressure * const.k_B) * (u.J * u.cm**-3).to(u.J * u.m**-3) # J * m**3 == Pa
        pressure = pressure.value
        part['pressure'] = np.array([pressure] * size)
        age = starPart['age'][mask] # in Myr
        constantCoveringFactor = self.config['constantCoveringFactor']
        if constantCoveringFactor:
            coveringFactor = np.float32(self.config['coveringFactor'])
            part['covering'] = np.array([coveringFactor] * size)
        else:
            PDRClearingTimescale = np.float32(self.config['PDRClearingTimescale'])
            part['covering'] = np.array(np.exp(-age / PDRClearingTimescale)) # from Baes, M., et al. 2024
            
        part['velocities'] = starPart['Velocities'][mask] * np.sqrt(self.a) # in km/s
        part['masses'] = starPart['Masses'][mask] * 10**10 / self.h # Msun
        
        header = 'Starforming regions\n' \
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

        info = np.column_stack((part['x'], part['y'], part['z'], part['smoothLength'],
                            part['sfr'], part['Z'], part['compactness'],
                            part['pressure'], part['covering'],
                            part['velocities'], part['masses'])) 

        np.savetxt(self.workingDir + '/starforming_regions.txt', info, header=header)
        
        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                        & (np.abs(g.pos[:, 1]) < region)\
                        & (np.abs(g.pos[:, 2]) < region)\
                        & (starPart['age'] > ageThreshold))[0]
        size = mask.shape[0]
        print('Stars: ', size)
        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask]
        part['z'] = g.pos[:, 2][mask]
        part['smoothLength'] = g.smoothLength[mask]
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
                
        info = np.column_stack((part['x'], part['y'], part['z'], part['smoothLength'],
                                part['initialMass'], part['Z'], part['age'],
                                part['velocities'], part['masses']))
        np.savetxt(self.workingDir + '/stars.txt', info, header=header)
        
        part = {}
        
        if self.config['includeDust']:
            
            fields = ['GFM_Metallicity', 'Coordinates', 'Masses', 
                      'InternalEnergy', 'StarFormationRate', 'ElectronAbundance', 
                      'Velocities', 'Density']
            
            if not self.config['requests']:
                gasPart = ill.snapshot.loadSubhalo(self.config['TNGPath'], self.config['snapNum'],
                                                   self.id, 'gas', fields)
            else:
                cutout_url = f"{self.base_url}{self.simulation}/snapshots/{self.config['snapNum']}/subhalos/{self.id}/cutout.hdf5"
                
                params = {'gas': ','.join(fields)}
                response = self.__make_request_with_retry(cutout_url, headers=self.headers, params=params)
                
                gasPart = self.__read_temp_file(response.content, partType='gas', fields=fields)
                
            if gasPart['count'] == 0:
                # print('No gas particles found, simulationMode falls to NoMedium')
                
                # self.config['simulationMode'] = 'NoMedium'
                # self.config['includeDust'] = False
                
                part['x'] = np.array([])
                part['y'] = np.array([])
                part['z'] = np.array([])
                part['mass'] = np.array([])
                part['Z'] = np.array([])
                part['Temperature'] = np.array([])
                part['velocities'] = np.array([])
                
            else:
                # gasPart['count'] = len(gasPart['Coordinates'])
                
                gasPart['Temperature'] = u2temp(gasPart['InternalEnergy'], 
                                                gasPart['ElectronAbundance'])
                gasPart['velocities'] = gasPart['Velocities'] * np.sqrt(self.a) # in km/s
                
                gas = PartNoSmoothing(gasPart['Coordinates'], gasPart['Masses'], 
                            self.pos, a=self.a, h=self.h)
                
                spatialmask = np.where((np.abs(gas.pos[:, 0]) < region)\
                            & (np.abs(gas.pos[:, 1]) < region)\
                            & (np.abs(gas.pos[:, 2]) < region))[0]
                
                DISMModel = self.config['DISMModel']
                if DISMModel == 'Camps_2016':
                    temperatureThreshold = np.float32(self.config['temperatureThreshold'])
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
                part['mass'] = gas.mass[mask]
                part['Z'] = gasPart['GFM_Metallicity'][mask]
                part['Temperature'] = gasPart['Temperature'][mask]
                part['velocities'] = gasPart['Velocities'][mask]
                
        else:
            part['x'] = np.array([])
            part['y'] = np.array([])
            part['z'] = np.array([])
            part['mass'] = np.array([])
            part['Z'] = np.array([])
            part['Temperature'] = np.array([])
            part['velocities'] = np.array([])
            
        header = 'Dusts\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n' \
                + 'Column 2: y-coordinate (kpc)\n' \
                + 'Column 3: z-coordinate (kpc)\n' \
                + 'Column 4: mass (Msun)\n' \
                + 'Column 5: metallicity (1)\n' \
                + 'Column 6: temperature (K)\n' \
                + 'Column 7: x-velocity (km/s)\n' \
                + 'Column 8: y-velocity (km/s)\n' \
                + 'Column 9: z-velocity (km/s)'
                
        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['mass'], part['Z'], part['Temperature'],
                                part['velocities']))

        np.savetxt(self.workingDir + '/dusts.txt', info, header=header)
        
    def __get_properties_survey(self, distance: float, properties: dict) -> dict:
        
        postProcessing = self.config['postProcessing']
        if postProcessing:
            surveys = self.config['surveys']
            spatialResols = []
            for survey in surveys:
                filters = self.config[f'filters_{survey}']
                numfilters = len(filters)
                if self.config[f'resolFromPix_{survey}']:
                    pixelScales = self.config[f'pixelScales_{survey}']
                    pixelScales = extend(pixelScales, numfilters)
                    resol = [distance * 10**6 * np.deg2rad(ang / 3600) for ang in pixelScales]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                else:
                    resol = np.float32(self.config[f'resolution_{survey}'])
                    resol = extend(resol, numfilters)
                    pixelScales = [np.rad2deg(res / (distance * 10**6)) * 3600 for res in resol]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                    
                properties[f'resolution_{survey}'] = resol
                properties[f'angleRes_{survey}'] = pixelScales
                properties[f'ps_in_sr_{survey}'] = ps_in_sr
                
                spatialResols += resol
                
                numExposure = np.int32(self.config[f'numExposure_{survey}'])
                numExposure = extend(numExposure, numfilters)
                properties[f'numExposure_{survey}'] = numExposure
                
                exposureTime = np.float32(self.config[f'exposureTime_{survey}'])
                exposureTime = extend(exposureTime, numfilters)
                properties[f'exposureTime_{survey}'] = exposureTime
                
                aperture = np.float32(self.config[f'aperture_{survey}'])
                properties[f'aperture_{survey}'] = aperture
                
                properties[f'numfilters_{survey}'] = numfilters
                properties[f'filters_{survey}'] = filters
                
            spatialResol = np.min(spatialResols)
            
            for survey in surveys:
                ratios = [spatialResol / res for res in properties[f'resolution_{survey}']]
                properties[f'ratios_{survey}'] = ratios
                
        else:
            spatialResol = np.float32(self.config['spatialResol'])
            
        properties['baseRes'] = spatialResol
        
        return properties
    
    def __get_properties(self) -> dict:
        
        properties = {}
        properties['subhaloID'] = self.id
        properties['stellarMass'] = np.log10(self.mass)
        properties['redshift'] = self.viewRedshift
        
        if self.config['inLocal']:
            distance = np.float32(self.config['viewDistance'])
            properties['cosmology'] = 'LocalUniverseCosmology'
        else:
            distance = self.cosmology.luminosity_distance(self.viewRedshift).value
            properties['cosmology'] = 'FlatUniverseCosmology'
            
        properties['lumiDis'] = distance
        
        properties = self.__get_properties_survey(distance, properties)
        
        faceAndEdge = self.config['faceAndEdge']
        
        inclinations = self.config['inclinations']
        azimuths = self.config['azimuths']
        numViews = np.int32(self.config['numViews'])
        
        estimateMorph = self.config['estimateMorph']
        if estimateMorph:
            galmorph = self.__estimate_morphology()
            properties['morph'] = galmorph
        
        if faceAndEdge:
            inclinations, azimuths = self.__calculate_angular_momentum_and_angles()
            numViews = 2
            self.config['numViews'] = numViews
            
        if not faceAndEdge and self.config['randomViews']:
            
            inclinations = np.random.uniform(0, 180, numViews)
            azimuths = np.random.uniform(-360, 360, numViews)
            
        self.config['inclinations'] = inclinations
        self.config['azimuths'] = azimuths
        
        for i, inc, azi in zip(range(numViews), inclinations, azimuths):
            print(f'View {i}: Inclination = {inc:.2f} deg, Azimuth = {azi:.2f} deg')
        
        properties['numViews'] = numViews
        properties['inclinations'] = inclinations
        properties['azimuths'] = azimuths
        
        boxLength = self.boxLength * 10**3 # in pc
        properties['boxLength'] = boxLength
        properties['boxlength_in_arcsec'] = np.rad2deg(boxLength / (distance * 10**6)) * 3600 # in arcsec
        
        return properties
        
    def __create_ski(self):
        
        print('Creating .ski file.')
        
        mode = self.config['simulationMode']
        
        ski_file = os.path.join(self.dataDir,  f'ski_templates/{mode}_template.ski')
        
        with open(ski_file, 'r') as file:
            data = file.read()
            
        if mode in ['DustEmission', 'ExtinctionOnly']:
            
            begin_str = '<VoronoiMeshMedium'
            end_str = '</VoronoiMeshMedium>'
            offset = len(end_str)
            
            idx_begin = data.index(begin_str)
            idx_end = data.index(end_str) + offset
            
            voronoiMeshMediumInfo = data[idx_begin: idx_end]
            remainingInfo = data[idx_end:]
            
            if self.config['hydrodynamicSolver'] == 'smoothParticle':
                replace_str = '<ParticleMedium filename="dusts.txt" '
                replace_str += 'massType="Mass" massFraction="0.3" '
                replace_str += 'importMetallicity="true" importTemperature="true" '
                replace_str += 'maxTemperature="0 K" importVelocity="false" importMagneticField="false" importVariableMixParams="false" useColumns="">\n'
                replace_str += '<smoothingKernel type="SmoothingKernel">\n'
                replace_str += '<CubicSplineSmoothingKernel/>\n'
                replace_str += '</smoothingKernel>\n'
                replace_str += '<materialMix type="MaterialMix">\n'
                replace_str += '<ZubkoDustMix numSilicateSizes="15" numGraphiteSizes="15" numPAHSizes="15"/>\n'
                replace_str += '</materialMix>\n'
                replace_str += '</ParticleMedium>\n'

                data = data.replace(voronoiMeshMediumInfo, replace_str)
                
        begin_str = '<cosmology'
        end_str = '</cosmology>'
        offset = len(end_str)
        
        idx_begin = data.index(begin_str)
        idx_end = data.index(end_str) + offset
        
        cosmologyInfo = data[idx_begin: idx_end]
        remainingInfo = data[idx_end:]
        
        if self.config['inLocal']:
            replace_str = '<cosmology type="Cosmology">\n'
            replace_str += '<LocalUniverseCosmology/>\n'
            replace_str += '</cosmology>\n'
            
            data = data.replace(cosmologyInfo, replace_str)
        else:
            data = data.replace('redshift="0.008"', f'redshift="{self.viewRedshift}"')
            
        properties = self.__get_properties()
        
        SEDFamily = {'BC03': 'BruzualCharlotSEDFamily',
                     'FSPS': 'FSPSSEDFamily'}
        SEDFamily = SEDFamily[self.config['SEDFamily']]
        
        if SEDFamily == 'FSPSSEDFamily':
            data = data.replace('resolution="High"', '')
        
        if mode in ['DustEmission', 'ExtinctionOnly']:
            dustEmissionType = self.config['dustEmissionType']
            data = data.replace('dustEmissionType="Equilibrium"',
                                f'dustEmissionType="{dustEmissionType}"')
        
        initialMassFunction = self.config['initialMassFunction']
        
        data = data.replace('BruzualCharlotSEDFamily', SEDFamily)
        data = data.replace('Chabrier', initialMassFunction)
        
        numPackets = np.float32(self.config['numPackets'])
        data = data.replace('numPackets="1e7"', f'numPackets="{numPackets}"')
        
        if self.config['inLocal']:
            minWavelength = np.float32(self.config['minWavelength'])
            maxWavelength = np.float32(self.config['maxWavelength'])
        else:
            minWavelength = np.float32(self.config['minWavelength']) * (1 + self.viewRedshift)
            maxWavelength = np.float32(self.config['maxWavelength']) * (1 + self.viewRedshift)
            
        data = data.replace('minWavelength="0.01 micron"', f'minWavelength="{minWavelength} micron"')
        data = data.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength} micron"')
        
        # print(minWavelength, maxWavelength)
        
        properties['minWavelength'] = minWavelength # micron
        properties['maxWavelength'] = maxWavelength # micron
        
        massFraction = np.float32(self.config['massFraction'])
        data = data.replace('massFraction="0.3"', f'massFraction="{massFraction}"')
        
        dustConfig = '<ZubkoDustMix numSilicateSizes="15" numGraphiteSizes="15" numPAHSizes="15"/>'
        dustModel = self.config['dustModel']
        
        numSilicateSizes = np.int32(self.config['numSilicateSizes'])
        numGraphiteSizes = np.int32(self.config['numGraphiteSizes'])
        numPAHSizes = np.int32(self.config['numPAHSizes'])
        numHydrocarbonSizes = np.int32(self.config['numHydrocarbonSizes'])
        
        if dustModel == 'ThemisDustMix':
            dustConfigNew = f'<{dustModel} numHydrocarbonSizes="{numHydrocarbonSizes}" numSilicateSizes="{numSilicateSizes}"/>'
        else:
            dustConfigNew = f'<{dustModel} numSilicateSizes="{numSilicateSizes}" numGraphiteSizes="{numGraphiteSizes}" numPAHSizes="{numPAHSizes}"/>'
            
        data = data.replace(dustConfig, dustConfigNew)
        
        minLevel = np.int32(self.config['minLevel'])
        maxLevel = np.int32(self.config['maxLevel'])
        data = data.replace('minLevel="8"', f'minLevel="{minLevel}"')
        data = data.replace('maxLevel="12"', f'maxLevel="{maxLevel}"')
        
        spatialRange = self.boxLength / 2. * 10**3 # in pc
        data = data.replace('minX="-5e4 pc"', f'minX="{-spatialRange} pc"')
        data = data.replace('maxX="5e4 pc"', f'maxX="{spatialRange} pc"')
        data = data.replace('minY="-5e4 pc"', f'minY="{-spatialRange} pc"')
        data = data.replace('maxY="5e4 pc"', f'maxY="{spatialRange} pc"')
        data = data.replace('minZ="-5e4 pc"', f'minZ="{-spatialRange} pc"')
        data = data.replace('maxZ="5e4 pc"', f'maxZ="{spatialRange} pc"')
        
        numWavelengths = np.int32(self.config['numWavelengths'])
        data = data.replace('numWavelengths="1000"', f'numWavelengths="{numWavelengths}"')
        
        properties['numWavelengths'] = numWavelengths
        
        begin_str = '<FullInstrument'
        end_str = '</FullInstrument>'
        offset = len(end_str)
        
        idx_begin = data.index(begin_str)
        idx_end = data.index(end_str) + offset
        
        instrumentInfo = data[idx_begin: idx_end]
        remainingInfo = data[idx_end:]
        
        data = data.replace(instrumentInfo, '')
        
        if self.config['inLocal']:
            distance = np.float32(self.config['viewDistance'])
        else:
            distance = self.cosmology.luminosity_distance(self.viewRedshift).value
            
        fieldOfView = np.float32(self.config['fieldOfView']) # in arcsec
        
        if fieldOfView == 0: 
            fovSize = self.boxLength * 10**3 # in pc
            fieldOfView = np.rad2deg(fovSize / (distance * 10**6)) * 3600 # in arcsec
        else:
            fovSize = distance * 10**6 * np.deg2rad(fieldOfView / 3600) # in pc
            
        properties['fovSize'] = fovSize
        properties['fieldOfView'] = fieldOfView
        
        if self.config['outputSEDOnly']:
            replace_str = '<SEDInstrument instrumentName="view" distance="0 Mpc" inclination="0 deg" azimuth="0 deg" '
            replace_str += 'roll="0 deg" radius="0 pc" recordComponents="false" numScatteringLevels="0" recordPolarization="false" recordStatistics="false">\n'
            replace_str += '<wavelengthGrid type="WavelengthGrid">\n'
            replace_str += '<LinWavelengthGrid minWavelength="0.1 micron" maxWavelength="1.2 micron" numWavelengths="200"/>\n'
            replace_str += '</wavelengthGrid>\n'
            replace_str += '</SEDInstrument>\n'
            
            replace_str = replace_str.replace('minWavelength="0.1 micron"', f'minWavelength="{minWavelength} micron"')
            replace_str = replace_str.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength} micron"')
            replace_str = replace_str.replace('numWavelengths="200"', f'numWavelengths="{numWavelengths}"')
            
            replace_str = replace_str.replace('radius="0 pc"', f'radius="{fovSize / 2} pc"')
            
            instrumentInfo = replace_str
            
        instrumentInfo = '\n' + instrumentInfo + '\n'
        
        spatialResol = properties['baseRes']
        inclinations = properties['inclinations']
        azimuths = properties['azimuths']
        
        numViews = properties['numViews']
        
        numPixels = int(fovSize / spatialResol)
        
        recordComponents = 'true' if self.config['recordComponents'] else 'false'
        
        insert_begin_idx = idx_begin
        for i, (inclination, azimuth) in enumerate(zip(inclinations, azimuths)):
            info = instrumentInfo.replace('view', f'view_{i:02d}')
            info = info.replace('inclination="0 deg"', f'inclination="{inclination} deg"')
            info = info.replace('azimuth="0 deg"', f'azimuth="{azimuth} deg"')
            info = info.replace('fieldOfViewX="1e5 pc"', f'fieldOfViewX="{fovSize} pc"')
            info = info.replace('fieldOfViewY="1e5 pc"', f'fieldOfViewY="{fovSize} pc"')
            info = info.replace('numPixelsX="1000"', f'numPixelsX="{numPixels}"')
            info = info.replace('numPixelsY="1000"', f'numPixelsY="{numPixels}"')
            info = info.replace('recordComponents="false"', f'recordComponents="{recordComponents}"')
            
            data = data[:insert_begin_idx] + info
            insert_begin_idx = insert_begin_idx + len(info)
            
        if self.config['inLocal']:
            data = data.replace('distance="0 Mpc"', f'distance="{distance} Mpc"')

        data = data + remainingInfo
        
        wavelengthGrid = self.config['wavelengthGrid']
        grid_type = {'Linear': 'LinWavelengthGrid',
                'Log': 'LogWavelengthGrid'}
        grid_type = grid_type[wavelengthGrid]
        data = data.replace('LinWavelengthGrid', grid_type)
        
        with open(self.workingDir + '/skirt.ski', 'w') as file:
            file.write(data)
        
        with open(self.workingDir + '/properties.json', 'w') as file:
            json.dump(properties, file, 
                      default=custom_serialier, indent=4)

        print('------estimate memory usage------')
        print(f'numViews: {numViews}')
        print(f'numSpatialPixels: {numPixels}')
        print(f'numWavelengthPixels: {numWavelengths}')
        
        factor = 7 if recordComponents == 'true' else 1
        numPixels = np.int64(numPixels) # avoid overflow
        dataCubeSize = np.int64(numPixels ** 2 * numWavelengths * numViews)
        dataSize_in_GB = np.around(dataCubeSize * 8 * factor * 1e-9, 3)
        print(f'Estimated memory usage: {dataSize_in_GB} GB')
        
    def __save_configs(self):
        
        with open(self.workingDir + '/config.json', 'w') as file:
            json.dump(self.config, file, default=custom_serialier, indent=4)
            
    def __estimate_morphology(self):
    
        '''
        Estimate the morphology of the galaxy by the S/T ratio.
        see https://arxiv.org/pdf/1904.12860, https://arxiv.org/pdf/2407.19152
        and https://www.tng-project.org/data/docs/specifications/#sec5c for more details.
        
        Returns:
            gal: 'spiral' or 'nonspiral'
        '''
        
        print('------Estimating morphology of galaxy------')
        
        stars = np.loadtxt(self.workingDir + '/stars.txt')
        
        headers = []
        with open(self.workingDir + '/stars.txt', 'r') as f:
            for line in f:
                if line.startswith('# Column'):
                    headers.append(line)
                    
        idx_coords = []
        idx_vels = []
        idx_masses = []
                    
        for i, header in enumerate(headers):
            if 'coordinate' in header.lower():
                idx_coords.append(i)
            elif 'velocity' in header.lower():
                idx_vels.append(i)
            elif 'Mass' in header:
                idx_masses.append(i)
                
        positions = stars[:, idx_coords]
        velocities = stars[:, idx_vels]
        masses = stars[:, idx_masses]
        
        offset = self.pos
        radius = self.radius
        
        boxsize = 3 * radius
        
        mask = np.where((np.abs(positions[:, 0]) < boxsize) \
                    & (np.abs(positions[:, 1]) < boxsize) \
                    & (np.abs(positions[:, 2]) < boxsize))[0]
        
        positions = positions[mask] + offset # original positions
        velocities = velocities[mask]
        masses = masses[mask].reshape(-1, 1)

        mass_center = np.sum(masses * positions, axis=0) / np.sum(masses)
        vel_center = np.sum(velocities * masses, axis=0) / np.sum(masses)
        
        dpos = positions - mass_center
        dvel = velocities - vel_center
        
        J3d_star = np.cross(dpos * masses, dvel)
        J_total = J3d_star.sum(axis=0)
        n_rotate = J_total / np.linalg.norm(J_total, keepdims=True, axis=-1)
        
        J_star = np.linalg.norm(J3d_star, axis=1)
        Jz_star = np.abs(J3d_star@n_rotate)
        
        disk_mask = J_star > 0
        disk_mask[disk_mask] = (Jz_star[disk_mask] / J_star[disk_mask]) > 0.7
        disk_mass = np.sum(masses[:, 0] * disk_mask)
        total_mass = np.sum(masses)
        sph_mass = total_mass - disk_mass
        
        S_T_ratio = sph_mass / total_mass
        
        print(f'S/T ratio: {S_T_ratio:.2f}')
        if S_T_ratio < 0.5:
            print('This galaxy is probably a spiral one.')
            gal = 'spiral'
        else:
            print('This galaxy is probably not a spiral one.')
            gal = 'nonspiral'

        return gal
        
            
    def prepare(self, arguments: Union[dict, NoneType]=None):
        
        '''
        Prepare the data and ski file for simulation
        
        Args:
            arguments: dictionary containing modifications for configurations to be used for the simulation
        '''
        
        if arguments is not None:
            exist_keys = self.config.keys()
            for key, value in arguments.items():
                if key in exist_keys:
                    self.config[key] = value
                    
        self.__init()
        self.__get_particles()
        self.__create_ski()
        self.__save_configs()
                    
    def inputs(self, data: dict):
        
        '''
        Input the data to create the ski file, used for simulations that are not TNG
        
        Args:
            data: dictionary containing data to be used for the simulation
        '''
        
        self.snapRedshift = data['snapRedshift']
        self.a = 1 / (1 + self.snapRedshift)
        self.cosmology = data['cosmology']
        self.mass = data['stellarMass']
        self.id = data['subhaloID']
        self.boxLength = data['boxLength']
        
        exist_keys = self.config.keys()
        for key, value in data.items():
            if key in exist_keys:
                self.config[key] = value
        
        mass_in_10 = np.log10(self.mass)
        
        print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [Msun]')
        
        self.__init()
        self.__create_ski()
        self.__save_configs()