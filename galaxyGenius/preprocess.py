import numpy as np
import h5py
import os
import illustris_python as ill
from astropy.cosmology import Planck15, Cosmology
import astropy.units as u
import astropy.constants as const
import requests
import json
import tempfile
from scipy.interpolate import interp1d
from typing import Union
import sys
import time
import gc
import numba
from typing import Union, Optional

from .utils import Units, u2temp, extend, custom_serialier, setup_logging, galaxygenius_data_dir

class PreProcess:
    
    def __init__(self, config: dict):
        
        self.config = config
        self.workingDir = self.config['workingDir']
        os.makedirs(self.workingDir, exist_ok=True)
        self.logger = setup_logging(os.path.join(self.workingDir, 'galaxyGenius.log'))
        self.logger.info(f'Initializing PreProcess class.')
        
        self.inputMethod = 'snapshot' # subhaloFile, snapshot, provided
        
        self.__init()
        self.__precompile_numba()
        
    def __precompile_numba(self):
        dummy_pos = np.ones((10, 3), dtype=np.float64)
        dummy_vel = np.ones((10, 3), dtype=np.float64)
        dummy_mass = np.ones(10, dtype=np.float64)
        self.angular_momentum(dummy_pos, dummy_vel, dummy_mass)
        
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
            self.logger.info('View redshift cannot be 0, setting to 0.005.')
            self.viewRedshift = 0.005
        
        if Units._instance is None:
            # if not defined, use Planck15
            self.cosmology = Planck15
            self.units = Units(cosmology=self.cosmology, 
                                snapRedshift=self.snapRedshift)
        else:
            # if defined, use cosmology from the units
            self.units = Units()
            self.cosmology = self.units.get_cosmology()
        
        self.logger.info(f'Use cosmology {self.cosmology.name} at redshift {self.snapRedshift}')
        
        self.h = self.cosmology.h
        self.a = 1 / (1 + self.snapRedshift)
        
        self.fage = self.__fage()
        
        self.dataDir = galaxygenius_data_dir() 
        self.simulation = self.config['simulation']
        self.snapnum = self.config['snapNum']
    
        self.__init_params()
        
        if 'TNG' in self.config['simulation']:
            name = 'tng'
        else:
            name = None
            self.logger.warning('Simulation name is unrecognized. Please manually input data by calling inputParticles() or inputs().')
            
        if self.config['requests']:
            self.base_url = f"https://www.{name}-project.org/api/"
            self.headers = {'api-key': self.config['apiKey']}
            
    def __make_request_with_retry(self, url: str, headers: dict = None, 
                                params: dict = None, max_retries: int = 5) -> requests.Response:
        
        start_time = time.time()
        content = None
        
        for attempt in range(max_retries):
            try:
                response = requests.get(url, headers=headers, params=params)
                response.raise_for_status()  # Raises an HTTPError for bad responses
                content = response
                break
            except requests.exceptions.RequestException as e:
                if attempt == max_retries - 1:  # Last attempt
                    self.logger.info(f"Error: Failed to make request after {max_retries} attempts.")
                    self.logger.info(f"URL: {url}")
                    self.logger.error(f"Error message: {str(e)}")
                else:
                    self.logger.info(f"Request failed (attempt {attempt + 1}/{max_retries}). Retrying...")
                    time.sleep(2 ** attempt)  # Exponential backoff        
                    
        end_time = time.time()
        self.logger.info(f"Requesting time taken: {end_time - start_time:.2f} seconds")
        
        if content is None:
            sys.exit(1)
        
        return content
                    
    def __read_subhalos(self) -> dict:
        snap_subhalos = ill.groupcat.loadSubhalos(self.config['TNGPath'], self.config['snapNum'])
        return snap_subhalos
    
    def __retrieve_subhalos_with_requests(self) -> list:
        
        cache_dir = os.path.join('.', 'cache')
        os.makedirs(cache_dir, exist_ok=True)
    
        minStellarMass = self.config['minStellarMass'].value / 10**10 * self.h
        maxStellarMass = self.config['maxStellarMass'].value / 10**10 * self.h
        
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
        
        """
        Retrieve subhalos from TNG snapshot or from requests.

        Returns:
            dict: A dictionary containing subhaloIDs, SFRs.
        """
        
        
        minStellarMass = self.config['minStellarMass'].to(u.Msun).value
        maxStellarMass = self.config['maxStellarMass'].to(u.Msun).value
        
        minStellarMass_in_10 = np.around(np.log10(minStellarMass), 2)
        maxStellarMass_in_10 = np.around(np.log10(maxStellarMass), 2)
        
        subhalos = {}
        
        if not self.config['requests']:
            
            snap_subhalos = self.__read_subhalos()
            stellarMass = (snap_subhalos['SubhaloMassType'][:, 4] * self.units.mass).to(u.Msun)
            
            subhalo_indices = np.where((stellarMass.value > minStellarMass) \
                                        & (stellarMass.value < maxStellarMass) \
                                        & (snap_subhalos['SubhaloFlag'] == 1))[0]
            
            self.subhaloIDs = subhalo_indices # indices are same as subhaloIDs
            self.stellarMasses = stellarMass[self.subhaloIDs] # with u.Msun
            halfMassRad = snap_subhalos['SubhaloHalfmassRadType'][:, 4] * self.units.distance
            self.halfMassRadii = halfMassRad[self.subhaloIDs] # with u.kpc
            subhaloNum = self.subhaloIDs.shape[0]
            subhaloPos = snap_subhalos['SubhaloPos'] * self.units.distance
            self.subhaloPos = subhaloPos[self.subhaloIDs]
            self.subhaloSFR = snap_subhalos['SubhaloSFR'][self.subhaloIDs] * self.units.sfr
            
            subhalos['subhaloNum'] = subhaloNum
            subhalos['subhaloIDs'] = self.subhaloIDs
            subhalos['subhaloSFR'] = self.subhaloSFR.to(u.Msun/u.yr)
            
            # subhalos['units'] = ['1', '1', 'Msun/yr']
            
            if maxStellarMass == np.inf:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass higher than 10^{minStellarMass_in_10:.2f} [M_sun]'
            else:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass from 10^{minStellarMass_in_10:.2f} to 10^{maxStellarMass_in_10:.2f} [M_sun]'
            
        else:
            
            self.logger.info('Using Web-based API to retrieve data.')
            
            results = self.__retrieve_subhalos_with_requests()
            self.results = results # to access in following functions
            
            subhaloNum = len(results)
            self.subhaloIDs = [result['id'] for result in results]
            self.subhaloSFR = [result['sfr'] for result in results] * self.units.sfr
            
            subhalos['subhaloNum'] = subhaloNum
            subhalos['subhaloIDs'] = self.subhaloIDs
            subhalos['subhaloSFR'] = self.subhaloSFR.to(u.Msun/u.yr)
            
            # subhalos['units'] = ['1', '1', 'Msun/yr']
            
            if maxStellarMass == np.inf:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass higher than 10^{minStellarMass_in_10:.2f} [M_sun]'
            else:
                print_info = f'{subhaloNum} subhalos in snapshot {self.snapnum} ' \
                    f'in stellar mass from 10^{minStellarMass_in_10:.2f} to 10^{maxStellarMass_in_10:.2f} [M_sun]'
        
        self.logger.info(print_info)
        
        return subhalos
    
    def subhalo(self, subhaloID: int):
        
        '''
        Specify one subhalo to be processed
        
        Args:
            subhaloID: subhaloID of the subhalo to be processed
        '''
        
        
        if not self.config['requests']:
            idx = list(self.subhaloIDs).index(subhaloID)
            self.id = self.subhaloIDs[idx]
            self.mass = self.stellarMasses[idx]
            self.radius = self.halfMassRadii[idx]
            self.pos = self.subhaloPos[idx]
            self.sfr = self.subhaloSFR[idx]

            mass_in_10 = np.log10(self.mass.to(u.Msun).value)
            
            self.logger.info(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [M_sun].')
            self.logger.info(f'SFR of Subhalo {self.id} is {self.sfr.to(u.Msun/u.yr).value:.2f} [M_sun/yr].')
            
        else:
            subhaloIDs = [result['id'] for result in self.results]
            idx = list(subhaloIDs).index(subhaloID)
            self.id = subhaloIDs[idx]
            subhalo = self.results[idx]
            subhalo_url = subhalo['url']
            self.logger.info(f'Subhalo URL: {subhalo_url}')
            subhalo_response = self.__make_request_with_retry(subhalo_url, headers=self.headers)
            self.data = subhalo_response.json()
            
            self.mass = self.data['mass_stars'] * self.units.mass
            self.radius = self.data['halfmassrad_stars'] * self.units.distance
            self.pos_x = self.data['pos_x'] * self.units.distance
            self.pos_y = self.data['pos_y'] * self.units.distance
            self.pos_z = self.data['pos_z'] * self.units.distance
            self.pos = u.Quantity([self.pos_x, self.pos_y, self.pos_z])
            self.sfr = self.data['sfr'] * self.units.sfr
            
            mass_in_10 = np.log10(self.mass.to(u.Msun).value)
            
            self.logger.info(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10:.2f} [M_sun].')
            self.logger.info(f'SFR of Subhalo {self.id} is {self.sfr.to(u.Msun/u.yr).value:.2f} [M_sun/yr].')

        self.partRegion = (self.config['boxLengthScale'] * self.radius).to(u.kpc) # boxlength, aperture
        self.partRegion = np.min([self.partRegion.value, self.config['maxBoxLength'].value])
        self.partRegion = self.partRegion * u.kpc
        
    @staticmethod
    @numba.njit(parallel=True, fastmath=True, cache=True)
    def angular_momentum(positions, velocities, masses):
        n = len(positions)
        partial = np.zeros((n, 3), dtype=positions.dtype)
        
        for i in numba.prange(n):
            
            r = positions[i]
            v = velocities[i]
            m = masses[i]

            lx = r[1] * m * v[2] - r[2] * m * v[1]
            ly = r[2] * m * v[0] - r[0] * m * v[2]
            lz = r[0] * m * v[1] - r[1] * m * v[0]
            
            partial[i, 0] = lx.item()
            partial[i, 1] = ly.item()
            partial[i, 2] = lz.item()

        return partial.sum(axis=0)
    
    def __calculate_angular_momentum_and_angles(self) -> tuple:
        
        self.logger.info('------Calculating face-on and edge-on viewing angles------')
        
        if self.inputMethod == 'input':
            particle_file = os.path.join(self.config['workingDir'], 
                                         'stars.txt')
            
            headers = []
            with open(particle_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        if 'column' in line.lower():
                            headers.append(line.strip().lower())
                    else:
                        break
                    
            idx_coordinates = [i for i, header in enumerate(headers) if 'coordinate' in header]
            idx_velocities = [i for i, header in enumerate(headers) if 'velocity' in header]
            idx_masses = [i for i, header in enumerate(headers) if 'mass' in header and 'initial' not in header]
            
            idx_columns = idx_coordinates + idx_velocities + idx_masses
            
            particles = np.loadtxt(particle_file, usecols=idx_columns)
            positions = particles[:, :3]
            velocities = particles[:, 3:6]
            masses = particles[:, 6]
        else:
            positions = self.starPart['Coordinates'].to(u.kpc).value
            velocities = self.starPart['Velocities'].to(u.km/u.s).value
            masses = self.starPart['Masses'].to(u.Msun).value
        
        # Initialize total angular momentum vector
        total_angular_momentum = np.zeros(3)
        
        mask = np.where((np.abs(positions[:, 0]) < 30) \
                    & (np.abs(positions[:, 1]) < 30) \
                    & (np.abs(positions[:, 2]) < 30))[0]
        
        positions = positions[mask]
        velocities = velocities[mask]
        masses = masses[mask]

        total_angular_momentum = self.angular_momentum(positions, velocities, masses)

        # Calculate angular momentum for each particle
        # for i in range(len(positions)):
        #     r = positions[i]
        #     v = velocities[i]
        #     m = masses[i]
            
        #     # Calculate angular momentum for the current particle
        #     angular_momentum = np.cross(r, m * v)
            
        #     # Add to total angular momentum
        #     total_angular_momentum += angular_momentum

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
        
        self.logger.info(f'Face-on angle (inc, azi): {incs[0]:.2f} deg, {azis[0]:.2f} deg')
        self.logger.info(f'Edge-on angle (inc, azi): {incs[1]:.2f} deg, {azis[1]:.2f} deg')
        
        return incs, azis
    
    def __read_temp_file(self, response_content, partType='star', params=None):

        try:
            with tempfile.NamedTemporaryFile(suffix='.hdf5', dir=self.workingDir, delete=True) as tmp:
                tmp.write(response_content)
                tmp.flush()
                
                with h5py.File(tmp.name, 'r') as f:
                    
                    if partType == 'star':
                        particle_data = {key: f[f'PartType4/{key}'][:] for key in params}
                    elif partType == 'gas':
                        particle_data = {key: f[f'PartType0/{key}'][:] for key in params}
                    
            particle_data['count'] = len(particle_data['Coordinates'])
                    
        except Exception as e:
            raise Exception(f'Error reading temporary file: {e}')
            
        return particle_data
    
    def __init_params(self):
        
        self.starParams = ['GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime',
                    'Coordinates', 'Velocities', 'Masses']
        self.starParamsUnits = [self.units.mass, self.units.dimless, self.units.dimless, 
                                self.units.distance, self.units.velocity, self.units.mass]

        
        if self.config['includeDust']:
            self.gasParams = ['GFM_Metallicity', 'Coordinates', 'Masses', 
                      'InternalEnergy', 'StarFormationRate', 'ElectronAbundance', 
                      'Velocities', 'Density']
            self.gasParamsUnits = [self.units.dimless, self.units.distance, self.units.mass, 
                                   self.units.energy, self.units.sfr, self.units.dimless, 
                                   self.units.velocity, self.units.density]
    
    def __deduce_unit(self, key):
        if 'coordinate' in key.lower():
            return self.units.distance
        elif 'x' in key.lower():
            return self.units.distance
        elif 'y' in key.lower():
            return self.units.distance
        elif 'z' in key.lower():
            return self.units.distance
        elif 'length' in key.lower():
            return self.units.distance
        elif 'radius' in key.lower():
            return self.units.distance
        elif 'density' in key.lower():
            return self.units.density
        elif 'mass' in key.lower():
            return self.units.mass
        elif 'sfr' in key.lower():
            return self.units.sfr
        elif 'velocity' in key.lower():
            return self.units.velocity
        elif 'potential' in key.lower():
            return self.units.potential
        elif 'temperature' in key.lower():
            return self.units.temperature
        elif 'energy' in key.lower():
            return self.units.energy
        elif 'metallicity' in key.lower():
            return self.units.dimless
        elif 'time' in key.lower():
            return self.units.dimless
        elif 'age' in key.lower():
            return self.units.time

    def __align_param_units(self, params, paramsUnits):
        
        seen = set()
        aligned_params = []
        aligned_paramsUnits = []
        for param, unit in zip(params, paramsUnits):
            if param not in seen:
                aligned_params.append(param)
                aligned_paramsUnits.append(u.Unit(unit)) # convert to astropy.unit.Unit
                seen.add(param)
                
        return aligned_params, aligned_paramsUnits
    
    def get_default_params(self):
        
        """
        Returns default parameters and their units for stars and gas (if includeDust is True).

        Parameters
        ----------
        None

        Returns
        -------
        returns : tuple
            Returns a tuple, where the first element is a dictionary of star parameters and their units,
            and the second element is a dictionary of gas parameters and their units (if includeDust is True).
            If includeDust is False, then the tuple contains only the star parameters and their units.
            The format of the dictionaries is {parameter: unit}.
        """
        
        
        starParamsUnits = [str(paramUnit) for paramUnit in self.starParamsUnits]
        starParamsUnit_dict = dict(zip(self.starParams, starParamsUnits))
        if self.config['includeDust']:
            gasParamsUnits = [str(paramUnit) for paramUnit in self.gasParamsUnits]             
            gasParamsUnit_dict = dict(zip(self.gasParams, gasParamsUnits))
            returns = (starParamsUnit_dict, gasParamsUnit_dict)
        else:
            returns = starParamsUnit_dict
        
        return returns
    
    def includeParams(self, starParams: Optional[list | dict] = None,
                          gasParams: Optional[list | dict] = None, **kwargs):
        """
        Adds parameters to the list of parameters to include in the data.

        Parameters
        ----------
        starParams : list or dict, optional
            List of parameters to include for stars. If a list is given, it should contain the
            parameter names as strings. If a dict is given, it should be in the format
            {parameter: unit}. In this case, the unit is not used to deduce the unit of the
            parameter but is directly used.
        gasParams : list or dict, optional
            List of parameters to include for gas particles. If a list is given, it should contain
            the parameter names as strings. If a dict is given, it should be in the format
            {parameter: unit}. In this case, the unit is not used to deduce the unit of the
            parameter but is directly used.
        **kwargs
            Additional keyword arguments are not used.

        Returns
        -------
        None
            This function does not return anything.
        """

        if starParams is not None:
            
            if isinstance(starParams, list):
                
                self.starParams = starParams
                self.starParamsUnits = [self.__deduce_unit(param) for param in starParams]
                
            if isinstance(starParams, dict):
                
                self.starParams = list(starParams.keys())
                self.starParamsUnits = list(starParams.values())
            
            self.starParams, self.starParamsUnits =\
                self.__align_param_units(self.starParams, self.starParamsUnits)
                
            starParamsUnits = [str(paramUnit) for paramUnit in self.starParamsUnits]
            starParamsUnit_dict = dict(zip(self.starParams, starParamsUnits))
            returns = starParamsUnit_dict
            
            self.logger.info(','.join(self.starParams) + ' for star particles included.')
        
        if self.config['includeDust']:
            if gasParams is not None:
                
                if isinstance(gasParams, list):
                    
                    self.gasParams = gasParams
                    self.gasParamsUnits = [self.__deduce_unit(param) for param in gasParams]
                    
                elif isinstance(gasParams, dict):
                    self.gasParams = list(gasParams.keys())
                    self.gasParamsUnits = list(gasParams.values())
            
                self.gasParams, self.gasParamsUnits =\
                    self.__align_param_units(self.gasParams, self.gasParamsUnits)
                    
                gasParamsUnits = [str(paramUnit) for paramUnit in self.gasParamsUnits]             
                gasParamsUnit_dict = dict(zip(self.gasParams, gasParamsUnits))
                returns = (starParamsUnit_dict, gasParamsUnit_dict)

            self.logger.info(','.join(self.gasParams) + ' for gas particles included.')
            
        return returns
    
    def __add_unit(self, parts, units):
        
        for i, param in enumerate(parts.keys()):
            
            parts[param] = np.array(parts[param]) * units[i]
            
        return parts
        
    def inputSubhaloParticleFile(self, subhaloParticleFile: str,
                                 subhaloInfo: dict, 
                                 centerPosition: Union[list, u.Quantity, None] = None):
        """
        Input subhalo particle file in h5 format.
        
        Args:
            subhaloParticleFile (str): path to the subhalo particle file in h5 format.
            subhaloInfo (dict): dictionary containing information about the subhalo, required: SubhaloID, stellarMass, halfStellarMassRadius.
            centerPosition (Union[list, u.Quantity, None], optional): center position of the galaxy. If None, the particles are assumed to be in center coordinates.
        """
        
        self.inputMethod = 'subhaloFile'
        
        self.logger.info('------Inputting subhalo particle file------')
        required_params = ['SubhaloID', 'stellarMass', 'halfStellarMassRadius']
        for par in required_params:
            if par not in subhaloInfo.keys():
                raise KeyError(f'Key {par} not found in subhalo info.')
        
        for key, value in subhaloInfo.items():
            if key == 'SubhaloID':
                self.id = value
            elif key == 'stellarMass':
                try:
                    self.mass = u.Quantity(value).to(u.Msun)
                except:
                    raise ValueError(f'Unable to convert {value} to u.Quantity.')
            elif key == 'halfStellarMassRadius':
                try:
                    self.radius = u.Quantity(value).to(u.kpc)
                except:
                    raise ValueError(f'Unable to convert {value} to u.Quantity.')
        
        self.partRegion = (self.config['boxLengthScale'] * self.radius).to(u.kpc)
        self.partRegion = np.min([self.partRegion.value, self.config['maxBoxLength'].to(u.kpc).value])
        self.partRegion = self.partRegion * u.kpc
    
        if not os.path.exists(subhaloParticleFile):
            raise FileNotFoundError(f'Subhalo particle file {subhaloParticleFile} does not exist.')
            
        with h5py.File(subhaloParticleFile, 'r') as file:
            
            for par in self.starParams:
                if f'PartType4/{par}' not in file:
                    raise KeyError(f'Key PartType4/{par} not found in subhalo particle file.')
            
            self.starPart = {key: file[f'PartType4/{key}'][:] for key in self.starParams}
            self.starPart = self.__add_unit(self.starPart, self.starParamsUnits)
            
            if self.config['includeDust']:
                
                for par in self.gasParams:
                    if f'PartType0/{par}' not in file:
                        raise KeyError(f'Key PartType0/{par} not found in subhalo particle file.')
                
                self.gasPart = {key: file[f'PartType0/{key}'][:] for key in self.gasParams}
                self.gasPart = self.__add_unit(self.gasPart, self.gasParamsUnits)
            
            if hasattr(self, 'dmParams'):

                for par in self.dmParams:
                    if f'PartType1/{par}' not in file:
                        raise KeyError(f'Key PartType1/{par} not found in subhalo particle file.')
                
                self.dmPart = {key: file[f'PartType1/{key}'][:] for key in self.dmParams}
                self.dmPart = self.__add_unit(self.dmPart, self.dmParamsUnits)
                
            if hasattr(self, 'bhParams'):
                
                for par in self.bhParams:
                    if f'PartType5/{par}' not in file:
                        raise KeyError(f'Key PartType5/{par} not found in subhalo particle file.')
                
                self.bhPart = {key: file[f'PartType5/{key}'][:] for key in self.bhParams}
                self.bhPart = self.__add_unit(self.bhPart, self.bhParamsUnits)
                
            self.logger.info(f'Particles from {subhaloParticleFile} loaded.')
            
        if centerPosition is not None:
            
            if isinstance(centerPosition, (u.Quantity, list)):
                if len(centerPosition) == 0 or len(centerPosition) != 3:
                    raise ValueError("Center position is empty or has wrong dimension.")
            
            if isinstance(centerPosition, list):
                
                if not all(isinstance(item, type(centerPosition[0])) for item in centerPosition):
                    raise TypeError("All elements in center position list must have the same type.")
                
                if isinstance(centerPosition[0], str):
                    try:
                        centerPosition = [u.Quantity(pos) for pos in centerPosition]
                        centerPosition = u.Quantity(centerPosition).to(u.kpc)
                    except:
                        raise ValueError("All elements in center position list must be convertable to u.kpc.")
                elif isinstance(centerPosition[0], u.Quantity):
                    try:
                        centerPosition = u.Quantity(centerPosition).to(u.kpc)
                    except:
                        raise ValueError("All elements in center position list must be convertable to u.kpc.")
                else:
                    raise TypeError('Type of center position list unrecognized.')
            
            elif isinstance(centerPosition, u.Quantity):
                try:
                    centerPosition = centerPosition.to(u.kpc)
                except:
                    raise ValueError("Center position must be convertable to u.kpc.")
                
        self.starPart = self.__in_box_mask(self.starPart, centerPosition)
        
        if self.config['includeDust']:
            self.gasPart = self.__in_box_mask(self.gasPart, centerPosition)
        
        if hasattr(self, 'dmParams'):
            self.dmPart = self.__in_box_mask(self.dmPart, centerPosition)
            # self.dmPart['Coordiantes'] = (self.dmPart['Coordinates'] - centerPosition).to(u.kpc)
            
        if hasattr(self, 'bhParams'):
            self.bhPart = self.__in_box_mask(self.bhPart, centerPosition)
            # self.bhPart['Coordinates'] = (self.bhPart['Coordinates'] - centerPosition).to(u.kpc)
            
    def __retrieveParts(self, partType, params):
        
        if not self.config['requests']:
            Part = ill.snapshot.loadSubhalo(
                self.config['TNGPath'], self.config['snapNum'],
                self.id, partType, params
            )
        else:
            cutout_url = f"{self.base_url}{self.simulation}/snapshots/{self.config['snapNum']}/subhalos/{self.id}/cutout.hdf5"
            self.logger.info(f'Retrieving {partType} particles.')
            read_params = params.copy()
            if partType == 'star':
                key = 'stars'
            elif partType == 'gas':
                key = 'gas'
            params = {f'{key}': ','.join(params)}
            
            response = self.__make_request_with_retry(cutout_url, headers=self.headers, params=params)
            
            Part = self.__read_temp_file(response.content, partType=partType, params=read_params)

        if 'count' in Part.keys():
            Part.pop('count')
        
        return Part
    
    def __get_particles(self):
        
        self.logger.info('Retrieving particles.')
        # if not self.input_subhalo_filename or not hasattr(self, 'input_subhalo_filename'):
        
        if self.inputMethod == 'snapshot':
            
            self.starPart = self.__retrieveParts('star', self.starParams)
            self.starPart = self.__add_unit(self.starPart, self.starParamsUnits)
            self.starPart = self.__in_box_mask(self.starPart, self.pos)
            
            if self.config['includeDust']:
                self.gasPart = self.__retrieveParts('gas', self.gasParams)
                self.gasPart = self.__add_unit(self.gasPart, self.gasParamsUnits)
                self.gasPart = self.__in_box_mask(self.gasPart, self.pos)
                
            if hasattr(self, 'dmParams'):
                self.dmPart = self.__retrieveParts('dm', self.dmParams)
                self.dmPart = self.__add_unit(self.dmPart, self.dmParamsUnits)
                self.dmPart = self.__in_box_mask(self.dmPart, self.pos)
                
            if hasattr(self, 'bhParams'):
                self.bhPart = self.__retrieveParts('bh', self.bhParams)
                self.bhPart = self.__add_unit(self.bhPart, self.bhParamsUnits)
                self.bhPart = self.__in_box_mask(self.bhPart, self.pos)
                
    def __get_smoothingLength(self):
        
        # from https://www.tng-project.org/about/
        if 'TNG50' in self.simulation:
            tng_smoothLength = 288 # in pc
        elif 'TNG100' in self.simulation:
            tng_smoothLength = 740 # in pc
        elif 'TNG300' in self.simulation:
            tng_smoothLength = 1480 # in pc
            
        # from https://www.tng-project.org/data/forum/topic/265/spatial-resolution-softening-length-mum-length/
        if self.snapRedshift <= 1:
            smoothLength = (tng_smoothLength * 10**-3 / self.a * self.h * self.units.distance).to(u.kpc)
            # smoothLength = (tng_smoothLength / self.a * self.h).to(u.kpc)
        else:
            smoothLength = (tng_smoothLength * 10**-3 / 0.5 * self.h * self.units.distance).to(u.kpc)
            # smoothLength = (tng_smoothLength / 0.5 * self.h).to(u.kpc)
        
        return smoothLength
    
    def __in_box_mask(self, particles, centerPosition=None):
        
        if centerPosition is not None:
            particles['Coordinates'] = (particles['Coordinates'].to(u.kpc) - centerPosition.to(u.kpc)).to(u.kpc)
        
        mask = np.where((np.abs(particles['Coordinates'][:, 0]) < self.partRegion / 2)\
                        & (np.abs(particles['Coordinates'][:, 1]) < self.partRegion / 2)\
                        & (np.abs(particles['Coordinates'][:, 2]) < self.partRegion / 2))[0]
        
        if 'count' in particles.keys():
            particles.pop('count')
        
        for key in particles.keys():
            particles[key] = particles[key][mask]
        
        return particles
    
    def __processing_particles(self):
        
        def starFunction(particles):
            
            smoothLength = self.__get_smoothingLength() # with kpc
            
            starFormationMask = np.where(particles['GFM_StellarFormationTime'] > 0)[0]
            snapshotAge = self.fage(self.snapRedshift) * u.Myr
            particles['age'] = snapshotAge - self.fage(1/particles['GFM_StellarFormationTime'] - 1) * u.Myr
            ageMask = np.where(particles['age'] > self.config['ageThreshold'])[0]
            
            idx = np.intersect1d(starFormationMask, ageMask)
            
            self.logger.info(f'Star particles: {len(idx)}')
            
            if len(idx) == 0:
                self.logger.info('No star particles found.')
                return {}
            
            properties = {}
            properties['x-coordinate'] = particles['Coordinates'][:, 0][idx]
            properties['y-coordinate'] = particles['Coordinates'][:, 1][idx]
            properties['z-coordinate'] = particles['Coordinates'][:, 2][idx]
            properties['smoothing length'] = np.full(idx.shape[0], smoothLength.value) * smoothLength.unit
            properties['initial mass'] = particles['GFM_InitialMass'][idx]
            properties['metallicity'] = particles['GFM_Metallicity'][idx]
            properties['age'] = particles['age'][idx]
            properties['x-velocity'] = particles['Velocities'][:, 0][idx]
            properties['y-velocity'] = particles['Velocities'][:, 1][idx]
            properties['z-velocity'] = particles['Velocities'][:, 2][idx]
            properties['mass'] = particles['Masses'][idx]
            
            return properties
        
        paramNames = ['x-coordinate', 'y-coordinate', 'z-coordinate', 'smoothing length',
                      'initial mass', 'metallicity', 'age', 'x-velocity', 'y-velocity', 'z-velocity', 'mass']
        
        self.createFile(paramNames, 'star', 
                        os.path.join(self.workingDir, 'stars.txt'), 
                        starFunction)
        
        def starFormingFunction(particles):
            
            smoothLength = self.__get_smoothingLength()
            
            starFormationMask = np.where(particles['GFM_StellarFormationTime'] > 0)[0]
            snapshotAge = self.fage(self.snapRedshift) * u.Myr
            particles['age'] = snapshotAge - self.fage(1/particles['GFM_StellarFormationTime'] - 1) * u.Myr
            ageMask = np.where(particles['age'] < self.config['ageThreshold'])[0]
            
            idx = np.intersect1d(starFormationMask, ageMask)
            
            self.logger.info(f'Star forming regions: {len(idx)}')
            
            if len(idx) == 0:
                self.logger.info('No star forming regions found.')
                return {}
            
            properties = {}
            properties['x-coordinate'] = particles['Coordinates'][:, 0][idx]
            properties['y-coordinate'] = particles['Coordinates'][:, 1][idx]
            properties['z-coordinate'] = particles['Coordinates'][:, 2][idx]
            properties['smoothing length'] = np.full(idx.shape[0], smoothLength.value) * smoothLength.unit
            properties['star formation rate'] = (particles['GFM_InitialMass'][idx] / self.config['ageThreshold']).to(u.Msun / u.yr)
            properties['metallicity'] = particles['GFM_Metallicity'][idx]
            # from Kapoor et al. 2021
            properties['compactness'] = np.random.normal(loc=self.config['logCompactnessMean'], 
                                                        scale=self.config['logCompactnessStd'], size=idx.shape[0]) * u.dimensionless_unscaled
            
            pressure = (10**self.config['logPressure'] * const.k_B * u.K * u.cm**-3).to(u.J / u.m**3) # J / m**3 == Pa
            properties['pressure'] = np.full(idx.shape[0], pressure.value) * pressure.unit
            if self.config['constantCoveringFactor']:
                properties['covering factor'] = np.full(idx.shape[0], self.config['coveringFactor']) * u.dimensionless_unscaled
            else:
                # from Baes, M., et al. 2024
                properties['covering factor'] = np.exp(-particles['age'][idx] / self.config['PDRClearingTimescale']) * u.dimensionless_unscaled
            
            properties['x-velocity'] = particles['Velocities'][:, 0][idx] # in km/s
            properties['y-velocity'] = particles['Velocities'][:, 1][idx] # in km/s
            properties['z-velocity'] = particles['Velocities'][:, 2][idx] # in km/s
            properties['mass'] = particles['Masses'][idx] # in Msun
            
            return properties
        
        paramNames = ['x-coordinate', 'y-coordinate', 'z-coordinate', 'smoothing length',
                      'star formation rate', 'metallicity', 'compactness', 'pressure', 'covering factor',
                      'x-velocity', 'y-velocity', 'z-velocity', 'mass']
        # paramUnits = ['kpc', 'kpc', 'kpc', 'kpc', 'Msun/yr',
        #               '1', '1', 'J/m3', '1', 'km/s', 'Msun']
        
        self.createFile(paramNames, 'starforming', 
                        os.path.join(self.workingDir, 'starforming_regions.txt'), 
                        starFormingFunction)
        
        if self.config['includeDust']:
            
            def dustFunction(particles):
                
                particles['Temperature'] = u2temp(particles['InternalEnergy'],
                                                particles['ElectronAbundance']) * u.K
                if self.config['DISMModel'] == 'Camps_2016':
                    temperatureThreshold = self.config['temperatureThreshold']
                    idx = np.where((particles['StarFormationRate'] > 0) \
                        | (particles['Temperature'] < temperatureThreshold))[0]
                elif self.config['DISMModel'] == 'Torrey_2012':
                    left = np.log10(particles['Temperature'].to(u.K).value)
                    right_log = np.log10(particles['Density'].to(10**10 * self.h**2 * u.Msun * u.kpc**-3).value)
                    right = 6 + 0.25 * right_log
                    idx = np.where(left < right)[0]
                
                self.logger.info(f'Dust: {len(idx)}')
                
                if len(idx) == 0:
                    self.logger.info('No dust found.')
                    return {}
                
                properties = {}
                properties['x-coordinate'] = particles['Coordinates'][:, 0][idx]
                properties['y-coordinate'] = particles['Coordinates'][:, 1][idx]
                properties['z-coordinate'] = particles['Coordinates'][:, 2][idx]
                properties['mass'] = particles['Masses'][idx]
                properties['metallicity'] = particles['GFM_Metallicity'][idx]
                properties['temperature'] = particles['Temperature'][idx]
                properties['x-velocity'] = particles['Velocities'][:, 0][idx]
                properties['y-velocity'] = particles['Velocities'][:, 1][idx]
                properties['z-velocity'] = particles['Velocities'][:, 2][idx]
                    
                return properties
        
            paramNames = ['x-coordinate', 'y-coordinate', 'z-coordinate', 'mass',
                          'metallicity', 'temperature', 'x-velocity', 'y-velocity', 'z-velocity']
            # paramUnits = ['kpc', 'kpc', 'kpc', 'Msun', '1', 'K', 'km/s']
            
            self.createFile(paramNames, 'dust', 
                            os.path.join(self.workingDir, 'dusts.txt'), 
                            dustFunction)
        
    def createFile(self, paramNames, partType, saveFilename, function):
        
        """
        Create a file containing the parameters of particles of a given type.

        Parameters
        ----------
        paramNames : list
            List of parameter names.
        partType : str
            Type of particles.
        saveFilename : str
            Path to save the file.
        function : function
            Function used to process particles.

        """
        
    
        if partType in ['star', 'stars', 'stellar', 'starforming', 'starformingRegions']:
            particles = self.starPart
        elif partType in ['gas', 'gases', 'dust']:
            particles = self.gasPart
        elif partType in ['dm', 'darkMatter']:
            if hasattr(self, 'dmPart'):
                particles = self.dmPart
            else:
                raise ValueError('DM particles not included.')
        elif partType in ['bh', 'bhs', 'blackhole', 'blackholes']:
            if hasattr(self, 'bhPart'):
                particles = self.bhPart
            else:
                raise ValueError('BH particles not included.')
        
        properties = function(particles)
        
        paramUnits = []
        for name in paramNames:
            unit = properties[name].unit
            paramUnits.append(unit)
            properties[name] = properties[name].to(unit).value
            size = properties[name].shape[0]

        header = f'{os.path.basename(saveFilename).split(".")[0]}\n'
        for i, (key, unit) in enumerate(zip(paramNames, paramUnits)):
            if unit == u.dimensionless_unscaled:
                unit = '1'
            elif 'solMass' in str(unit):
                unit = str(unit).replace('solMass', 'Msun')
            unit = str(unit).replace(' ', '') # remove empty space
            header = header + f'\ncolumn {i + 1}: {key} ({unit})'
        
        if len(properties) == 0:
            arr_size = (0, len(paramNames))
        else:
            arr_size = (size, len(paramNames))
            
        info_array = np.zeros(arr_size)
        for i, key in enumerate(properties):
            info_array[:, i] = properties[key]
        
        np.savetxt(f'{saveFilename}', info_array, header=header)
    
    
    def __get_properties_survey(self, distance: Union[float, u.Quantity], properties: dict) -> dict:
        
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
                    # note that u.rad is not dimensionless quantity
                    resol = u.Quantity([(distance * ang.to(u.rad).value).to(u.pc) for ang in pixelScales])
                    ps_in_sr = u.Quantity([(ang ** 2).to(u.sr) for ang in pixelScales])
                else:
                    resol = self.config[f'resolution_{survey}']
                    resol = extend(resol, numfilters)
                    pixelScales = u.Quantity([(res / distance * u.rad).to(u.deg) for res in resol])
                    ps_in_sr = u.Quantity([(ang**2).to(u.sr) for ang in pixelScales])
                    
                properties[f'resolution_{survey}'] = resol
                properties[f'angleRes_{survey}'] = pixelScales
                properties[f'ps_in_sr_{survey}'] = ps_in_sr
                
                spatialResols += resol.to(u.pc).value.tolist()
                
                numExposure = self.config[f'numExposure_{survey}']
                numExposure = extend(numExposure, numfilters)
                properties[f'numExposure_{survey}'] = numExposure.value
                
                exposureTime = self.config[f'exposureTime_{survey}']
                exposureTime = extend(exposureTime, numfilters)
                properties[f'exposureTime_{survey}'] = exposureTime
                
                aperture = self.config[f'aperture_{survey}']
                properties[f'aperture_{survey}'] = aperture
                
                properties[f'numfilters_{survey}'] = numfilters
                properties[f'filters_{survey}'] = filters
            
            spatialResol = np.min(spatialResols) * u.pc
            
            for survey in surveys:
                ratios = [(spatialResol / res).value for res in properties[f'resolution_{survey}']]
                properties[f'ratios_{survey}'] = ratios
                
        else:
            spatialResol = self.config['spatialResol']
            
        properties['baseRes'] = spatialResol
        
        return properties
    
    def __get_properties(self) -> dict:
        
        properties = {}
        properties['snapRedshift'] = self.snapRedshift
        properties['snapNum'] = self.snapnum
        properties['subhaloID'] = self.id
        properties['stellarMass'] = self.mass
        properties['radius'] = self.radius
        
        if self.config['inLocal']:
            distance = self.config['viewDistance']
            properties['cosmology'] = 'LocalUniverseCosmology'
            properties['viewRedshift'] = None
        else:
            distance = self.cosmology.luminosity_distance(self.viewRedshift)
            properties['cosmology'] = 'FlatUniverseCosmology'
            properties['viewRedshift'] = self.viewRedshift
            
        properties['lumiDis'] = distance
        
        properties = self.__get_properties_survey(distance, properties)
        
        faceAndEdge = self.config['faceAndEdge']
        
        inclinations = self.config['inclinations']
        azimuths = self.config['azimuths']
        numViews = self.config['numViews']
        
        estimateMorph = self.config['estimateMorph']
        if estimateMorph:
            galmorph = self.__estimate_morphology()
            properties['morph'] = galmorph
        
        if faceAndEdge:
            
            inclinations, azimuths = self.__calculate_angular_momentum_and_angles()
            numViews = 2
            self.config['numViews'] = numViews
            
            self.logger.info('Using calculated face-on and edge-on angles')
            
        if not faceAndEdge and self.config['randomViews']:
            
            inclinations = np.random.uniform(0, 180, numViews)
            azimuths = np.random.uniform(-360, 360, numViews)
            numViews = len(inclinations)
            
            
            self.logger.info(f'Using {numViews} random views')
            
        elif not faceAndEdge and not self.config['randomViews']:
            
            inclinations = self.config['inclinations']
            azimuths = self.config['azimuths']
            numViews = len(inclinations)
            
            self.logger.info(f'Using {numViews} specified views')
            
        self.config['inclinations'] = inclinations
        self.config['azimuths'] = azimuths
        self.config['numViews'] = numViews
        
        for i, inc, azi in zip(range(numViews), inclinations, azimuths):
            self.logger.info(f'View {i}: Inclination = {inc:.2f} deg, Azimuth = {azi:.2f} deg')
        
        properties['numViews'] = numViews
        properties['inclinations'] = inclinations
        properties['azimuths'] = azimuths
        
        properties['boxLength'] = self.partRegion
        properties['boxlength_in_arcsec'] = (self.partRegion / distance * u.rad).to(u.arcsec) 
        # properties['boxlength_in_arcsec'] = np.rad2deg(boxLength / (distance * 10**6)) * 3600 # in arcsec
        
        return properties
        
    def __create_ski(self):
        
        self.logger.info('Creating .ski file.')
        
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
        
        numPackets = self.config['numPackets']
        data = data.replace('numPackets="1e7"', f'numPackets="{numPackets}"')
        
        if self.config['inLocal']:
            minWavelength = self.config['minWavelength'].to(u.um)
            maxWavelength = self.config['maxWavelength'].to(u.um)
        else:
            minWavelength = self.config['minWavelength'].to(u.um) * (1 + self.viewRedshift)
            maxWavelength = self.config['maxWavelength'].to(u.um) * (1 + self.viewRedshift)
            
        data = data.replace('minWavelength="0.01 micron"', f'minWavelength="{minWavelength.value} micron"')
        data = data.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength.value} micron"')
    
        properties['minWavelength'] = minWavelength # micron
        properties['maxWavelength'] = maxWavelength # micron
        
        massFraction = self.config['massFraction']
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
        
        spatialRange = self.partRegion.to(u.pc).value / 2
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
        
        distance = properties['lumiDis']
        
        fieldOfView = self.config['fieldOfView'] # in arcsec
        
        if fieldOfView == 0:
            fovSize = self.partRegion.to(u.pc) # in pc
            fieldOfView = (fovSize / distance * u.rad).to(u.arcsec)
            # fieldOfView = np.rad2deg(fovSize / (distance * 10**6)) * 3600 # in arcsec
        else:
            fovSize = (distance * fieldOfView.to(u.rad).value).to(u.pc)
            # fovSize = distance * 10**6 * np.deg2rad(fieldOfView / 3600) # in pc
            
        properties['fovSize'] = fovSize
        properties['fieldOfView'] = fieldOfView
        
        if self.config['outputSEDOnly']:
            replace_str = '<SEDInstrument instrumentName="view" distance="0 Mpc" inclination="0 deg" azimuth="0 deg" '
            replace_str += 'roll="0 deg" radius="0 pc" recordComponents="false" numScatteringLevels="0" recordPolarization="false" recordStatistics="false">\n'
            replace_str += '<wavelengthGrid type="WavelengthGrid">\n'
            replace_str += '<LinWavelengthGrid minWavelength="0.1 micron" maxWavelength="1.2 micron" numWavelengths="200"/>\n'
            replace_str += '</wavelengthGrid>\n'
            replace_str += '</SEDInstrument>\n'
            
            replace_str = replace_str.replace('minWavelength="0.1 micron"', f'minWavelength="{minWavelength.value} micron"')
            replace_str = replace_str.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength.value} micron"')
            replace_str = replace_str.replace('numWavelengths="200"', f'numWavelengths="{numWavelengths}"')
            
            replace_str = replace_str.replace('radius="0 pc"', f'radius="{fovSize.value / 2} pc"')
            
            instrumentInfo = replace_str
            
        instrumentInfo = '\n' + instrumentInfo + '\n'
        
        spatialResol = properties['baseRes']
        inclinations = properties['inclinations']
        azimuths = properties['azimuths']
        
        numViews = properties['numViews']
        
        numPixels = int(fovSize / spatialResol)
        
        properties['numPixels'] = numPixels
        
        recordComponents = 'true' if self.config['recordComponents'] else 'false'
        
        insert_begin_idx = idx_begin
        for i, (inclination, azimuth) in enumerate(zip(inclinations, azimuths)):
            info = instrumentInfo.replace('view', f'view_{i:02d}')
            info = info.replace('inclination="0 deg"', f'inclination="{inclination} deg"')
            info = info.replace('azimuth="0 deg"', f'azimuth="{azimuth} deg"')
            info = info.replace('fieldOfViewX="1e5 pc"', f'fieldOfViewX="{fovSize.value} pc"')
            info = info.replace('fieldOfViewY="1e5 pc"', f'fieldOfViewY="{fovSize.value} pc"')
            info = info.replace('numPixelsX="1000"', f'numPixelsX="{numPixels}"')
            info = info.replace('numPixelsY="1000"', f'numPixelsY="{numPixels}"')
            info = info.replace('recordComponents="false"', f'recordComponents="{recordComponents}"')
            
            data = data[:insert_begin_idx] + info
            insert_begin_idx = insert_begin_idx + len(info)
            
        if self.config['inLocal']:
            data = data.replace('distance="0 Mpc"', f'distance="{distance.value} Mpc"')

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

        self.logger.info('------estimate memory usage------')
        self.logger.info(f'numViews: {numViews}')
        self.logger.info(f'numSpatialPixels: {numPixels}')
        self.logger.info(f'numWavelengthPixels: {numWavelengths}')
        
        factor = 7 if recordComponents == 'true' else 1
        numPixels = np.int64(numPixels) # avoid overflow
        dataCubeSize = np.int64(numPixels ** 2 * numWavelengths * numViews)
        dataSize_in_GB = np.around(dataCubeSize * 8 * factor * 1e-9, 3)
        self.logger.info(f'Estimated memory usage: {dataSize_in_GB} GB')
        
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
        
        self.logger.info('------Estimating morphology of galaxy------')
        
        if self.inputMethod == 'input':
            particle_file = os.path.join(self.config['workingDir'], 
                                         'stars.txt')
            
            headers = []
            with open(particle_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        headers.append(line.strip().lower())
                    else:
                        break
                    
            idx_coordinates = [i for i, header in enumerate(headers) if 'coordinate' in header]
            idx_velocities = [i for i, header in enumerate(headers) if 'velocity' in header]
            idx_masses = [i for i, header in enumerate(headers) if 'mass' in header and 'initial' not in header]
            
            idx_columns = idx_coordinates + idx_velocities + idx_masses
            
            particles = np.loadtxt(particle_file, usecols=idx_columns)
            positions = particles[:, :3]
            velocities = particles[:, 3:6]
            masses = particles[:, 6]
            
        else:
            positions = self.starPart['Coordinates'].to(u.kpc).value
            velocities = self.starPart['Velocities'].to(u.km/u.s).value
            masses = self.starPart['Masses'].to(u.Msun).value
        
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
        
        self.logger.info(f'S/T ratio: {S_T_ratio:.2f}')
        if S_T_ratio < 0.5:
            self.logger.info('This galaxy is probably a spiral one.')
            gal = 'spiral'
        else:
            self.logger.info('This galaxy is probably not a spiral one.')
            gal = 'nonspiral'

        return gal
    
    def inputParticles(self, partType: str, particles: dict, 
                       subhaloInfo: dict, cosmology: Union[Cosmology, None] = None,
                       centerPosition: Union[list, u.Quantity, None] = None):
        
        self.inputMethod = 'partInput'
        self.logger.info(f'------Inputing {partType} particles------')
        required_params = ['SubhaloID', 'stellarMass', 'halfStellarMassRadius']
        for par in required_params:
            if par not in subhaloInfo.keys():
                raise KeyError(f'Key {par} not found in subhalo info.')
        
        for key, value in subhaloInfo.items():
            if key == 'SubhaloID':
                self.id = value
            elif key == 'stellarMass':
                try:
                    self.mass = u.Quantity(value).to(u.Msun)
                except:
                    raise ValueError(f'Unable to convert {value} to u.Quantity.')
            elif key == 'halfStellarMassRadius':
                try:
                    self.radius = u.Quantity(value).to(u.kpc)
                except:
                    raise ValueError(f'Unable to convert {value} to u.Quantity.')
        
        self.partRegion = (self.config['boxLengthScale'] * self.radius).to(u.kpc)
        self.partRegion = np.min([self.partRegion.value, self.config['maxBoxLength'].to(u.kpc).value])
        self.partRegion = self.partRegion * u.kpc
        
        for key, value in particles.items():
            try:
                unit = value.unit
                particles[key] = np.array(value) * unit # convert to numpy array with unit
            except:
                raise ValueError(f'{key} data is not a u.Quantity object.')
            
        if centerPosition is not None:
            
            if isinstance(centerPosition, (u.Quantity, list)):
                if len(centerPosition) == 0 or len(centerPosition) != 3:
                    raise ValueError("Center position is empty or has wrong dimension.")

            if isinstance(centerPosition, list):
                
                if not all(isinstance(item, type(centerPosition[0])) for item in centerPosition):
                    raise TypeError("All elements in center position list must have the same type.")
                
                if isinstance(centerPosition[0], str):
                    try:
                        centerPosition = [u.Quantity(pos) for pos in centerPosition]
                        centerPosition = u.Quantity(centerPosition).to(u.kpc)
                    except:
                        raise ValueError("All elements in center position list must be convertable to u.kpc.")
                elif isinstance(centerPosition[0], u.Quantity):
                    try:
                        centerPosition = u.Quantity(centerPosition).to(u.kpc)
                    except:
                        raise ValueError("All elements in center position list must be convertable to u.kpc.")
                else:
                    raise TypeError('Type of center position list unrecognized.')
                
            elif isinstance(centerPosition, u.Quantity):
                try:
                    centerPosition = centerPosition.to(u.kpc)
                except:
                    raise ValueError("Center position must be convertable to u.kpc.")
        
        if cosmology is not None:
            try:
                self.cosmology = cosmology
                h = self.cosmology.h
            except:
                raise ValueError('cosmology is not an astropy cosmology object.')

        if partType in ['star', 'stars']:
            self.starPart = particles
            self.starPart = self.__in_box_mask(self.starPart, centerPosition)
            
        if self.config['includeDust']:
            if partType in ['gas', 'gases', 'dust', 'dusts']:
                self.gasPart = particles
                self.gasPart = self.__in_box_mask(self.gasPart, centerPosition)

    def modify_configs(self, arguments: dict):
        
        exist_keys = self.config.keys()
        for key, value in arguments.items():
            if key in exist_keys:
                original_value = self.config[key]
                if isinstance(original_value, u.Quantity):
                    try:
                        value = u.Quantity(value)
                        unit = original_value.unit
                        value = value.to(unit)
                    except:
                        raise ValueError(f'{value} cannot be converted to u.Quantity with the same unit as {original_value}')
                self.config[key] = value
                self.logger.info(f'Changed {key} from {original_value} to {value}')
    
    def inputs(self, data: dict):
        
        self.inputMethod = 'input'
        
        # in case changing configs
        # exist_keys = self.config.keys()
        # for key, value in data.items():
        #     if key in exist_keys:
        #         self.config[key] = value
        
        self.modify_configs(data)
        
        required_params = ['SubhaloID', 'stellarMass', 'halfStellarMassRadius']
        for par in required_params:
            if par not in data.keys():
                raise ValueError(f'{par} is required but not found in the input data.')
            
            if par == 'SubhaloID':
                self.id = data['SubhaloID']
            
            if par == 'stellarMass':
                try:
                    self.mass = u.Quantity(data['stellarMass']).to(u.Msun)
                except:
                    raise ValueError('stellarMass is not with astropy Unit.')
                    
            if par == 'halfStellarMassRadius':
                try:
                    self.radius = u.Quantity(data['halfStellarMassRadius']).to(u.kpc)
                except:
                    raise ValueError('halfStellarMassRadius is not with astropy Unit.')
                
        self.partRegion = (self.config['boxLengthScale'] * self.radius).to(u.kpc)
        self.partRegion = np.min([self.partRegion.value, self.config['maxBoxLength'].value])
        self.partRegion = self.partRegion * u.kpc
        
        if not os.path.join(self.workingDir, 'stars.txt'):
            raise FileNotFoundError(f'stars.txt not found in {self.workingDir}.')
        
        if not os.path.join(self.workingDir, 'starforming_regions.txt'):
            raise FileNotFoundError(f'starforming_regions.txt not found in {self.workingDir}.')
        
        if self.config['includeDust']:
            if not os.path.join(self.workingDir, 'dusts.txt'):
                raise FileNotFoundError(f'dusts.txt not found in {self.workingDir}.')
    
    def prepare(self, arguments: Union[dict, None]=None):
        
        # if arguments is not None:
        #     exist_keys = self.config.keys()
        #     for key, value in arguments.items():
        #         if key in exist_keys:
        #             self.config[key] = value
        if arguments is not None:
            self.modify_configs(arguments) 
        
        self.__init()
        if self.inputMethod == 'snapshot':
            self.__get_particles()
            self.__processing_particles()
        elif self.inputMethod == 'partInput':
            data = {'SubhaloID': self.id,
                    'stellarMass': self.mass,
                    'halfStellarMassRadius': self.radius}
            self.inputs(data)
        elif self.inputMethod == 'subhaloFile':
            self.__processing_particles()
        elif self.inputMethod == 'input':
            pass
        
        self.__create_ski()
        self.__save_configs()