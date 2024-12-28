import numpy as np
import h5py
import os
import illustris_python as ill
from astropy.cosmology import Planck15
import astropy.units as u
import subprocess
from .utils import *
from scipy.spatial import cKDTree
import sys
import pickle

class Galaxy:

    def __init__(self, position, smoothLength, mass, subhaloPos, a, h):
        # convert comoving distance (ckpc/h) to physical distance (kpc)
        self.pos = position * a / h - subhaloPos # kpc 
        self.smoothLength = smoothLength * a / h # kpc
        self.mass = mass * 10**10 / h # Msun
        

class PreProcess:
    
    def __init__(self, config):
        self.config = config
        if self.config['flags'] > 0:
            print('Conflictions in configuration. Please edit them !')
            sys.exit()
        
        self.snapnum = np.int32(self.config['snapNum'])
        self.cosmology = Planck15
        self.snapz = self.__get_snapz()
        self.h = self.cosmology.h
        self.a = 1 / (1 + self.snapz)
        self.workingDir = self.config['workingDir']
        self.dataDir = self.config['dataDir']
        
        # logger = get_logger('logger', log_file="log.txt")

    def __get_snapz(self):
        snap = h5py.File(os.path.join(self.config['filePath'],
                                      f'snapdir_{self.snapnum:03d}/snap_{self.snapnum:03d}.0.hdf5'), 'r')
        snapz = dict(snap['Header'].attrs.items())['Redshift']
        return snapz
        
    def __read_subhalos(self):
        snap_subhalos = ill.groupcat.loadSubhalos(self.config['filePath'], self.snapnum)
        return snap_subhalos

    def get_subhalos(self):

        minStellarMass = np.float32(self.config['minStellarMass'])
        if self.config['maxStellarMass'] == 'inf':
            maxStellarMass = np.inf
        else:
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

        minStellarMass_in_10 = np.around(np.log10(minStellarMass), 2)
        maxStellarMass_in_10 = np.around(np.log10(maxStellarMass), 2)
        
        if maxStellarMass == np.inf:
            print(f'{self.subhaloNum} subhalos in snapshot {self.snapnum} in stellar mass higher than 10^{minStellarMass_in_10} [M_sun]')
        else:
            print(f'{self.subhaloNum} subhalos in snapshot {self.snapnum} in stellar mass from 10^{minStellarMass_in_10} to 10^{maxStellarMass_in_10} [M_sun]')
 
    def get_subhaloIDs(self):
        return self.subhaloIDs

    def get_stellarMasses(self):
        return self.stellarMasses * u.Msun

    def subhalo(self, subhaloID):

        self.idx = list(self.subhaloIDs).index(subhaloID)
        self.id = self.subhaloIDs[self.idx]
        self.mass = self.stellarMasses[self.idx]
        self.radius = self.halfMassRadii[self.idx]
        self.pos = self.subhaloPos[self.idx]

        mass_in_10 = np.log10(self.mass)
        print(f'Stellar Mass of Subhalo {self.id} is 10^{mass_in_10} [M_sun].')
        
    def __get_TNG(self):
        filePath = self.config['filePath']
        if filePath.endswith('/'):
            filePath = '/'.join(filePath.split('/')[:-1])
        tng = filePath.split('/')[-1]
        return tng

    def __get_particles(self):
        
        print('Retrieving Stellar and Gas particles.')

        self.boxLength = self.radius * np.float32(self.config['boxLengthScale'])
        self.boxLength = np.min([self.boxLength, np.float32(self.config['maxBoxLength'])])

        region = self.boxLength / 2.

        fields = ['GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime',
                'Coordinates']
        
        starPart = ill.snapshot.loadSubhalo(self.config['filePath'], self.snapnum,
                                            self.id, 'star', fields)
        
        if not self.config['smoothLengthFromNeighbors']:
            # see https://www.tng-project.org/data/docs/background/ for smoothing length
            # see https://www.tng-project.org/data/forum/topic/265/spatial-resolution-softening-length-minimum-length/
            stellar_smoothLength_dict = {'TNG-50': 288, 'TNG-100': 740, 'TNG-300': 1480} # in pc, physical
            tng = self.__get_TNG()
            stellar_smoothLength = stellar_smoothLength_dict[tng]

            if self.snapz <= 1:
                # z = 0 ~ 1, smoothing length are changing
                stellar_smoothLength = stellar_smoothLength * 10**-3 / self.a * self.h # in ckpc/h, comoving scale
            else:
                # z > 1, smoothing length are constant as z = 1 value
                stellar_smoothLength = stellar_smoothLength * 10**-3 / 0.5 * self.h # in ckpc/h, comoving scale
            
            starPart['smoothLength'] = np.array([stellar_smoothLength] * starPart['count'])
            
        else:
            
            # obtained from the nearest neighbor
            ckdtree = cKDTree(starPart['Coordinates'])
            k = np.int32(self.config['NthNeighbor']) + 1 
            distances, _ = ckdtree.query(starPart['Coordinates'], k=k) # in ckpc/h
            distances_to_32nd_neighbor = distances[:, k]
            starPart['smoothLength'] = distances_to_32nd_neighbor
            
            max_smooth_length = np.float32(self.config['maxSmoothingLength']) / self.a * self.h # kpc to ckpc/h
            starPart['smoothLength'][starPart['smoothLength'] > max_smooth_length] = max_smooth_length
        
        mask = np.where(starPart['GFM_StellarFormationTime'] > 0)[0]
        for key in starPart.keys():
            if key == 'count':
                pass
            else:
                starPart[key] = starPart[key][mask] # select stellar particles instead of wind
        
        snapshot_age = fage(self.snapz)
        starPart['age'] = snapshot_age - fage(1/starPart['GFM_StellarFormationTime'] - 1) # in Myr
        g = Galaxy(starPart['Coordinates'], starPart['smoothLength'], starPart['GFM_InitialMass'],
                   self.pos, a=self.a, h=self.h)


        ageThreshold = np.float32(self.config['ageThreshold'])

        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                    & (np.abs(g.pos[:, 1]) < region)\
                    & (np.abs(g.pos[:, 2]) < region)\
                    & (starPart['age'] < ageThreshold))[0]
        size = mask.shape[0]
        print('MAPPING III particles:', size)

        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask]
        part['z'] = g.pos[:, 2][mask]
        part['smoothLength'] = g.smoothLength[mask] # smoothing length in kpc
        part['sfr'] = g.mass[mask] / (ageThreshold * 10**6) # constant SFR in Msun/yr
        part['Z'] = starPart['GFM_Metallicity'][mask]
        part['compactness'] = np.random.normal(loc=np.float32(self.config['logCompactnessMean']), 
                                               scale=np.float32(self.config['logCompactnessStd']), size=size) # from Kapoor et al. 2021
        pressure = (10**np.float32(self.config['logPressure']) * const.k_B) * (u.J * u.cm**-3).to(u.J * u.m**-3) # J * m**3 == Pa
        pressure = pressure.value
        part['pressure'] = np.array([pressure] * size)
        age = starPart['age'][mask] # in Myr
        if self.config['constantCoveringFactor']:
            part['covering'] = np.array([np.float32(self.config['coveringFactor'])] * size)
        else:
            part['covering'] = np.array(np.exp(-age / np.float32(self.config['PDRClearingTimescale']))) # from Baes, M., et al. 2024

        header = 'star forming particles\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n' \
                + 'Column 2: y-coordinate (kpc)\n' \
                + 'Column 3: z-coordinate (kpc)\n' \
                + 'Column 4: smoothing length (kpc)\n' \
                + 'Column 5: star formation rate (Msun/yr)\n' \
                + 'Column 6: metallicity (1)\n' \
                + 'Column 7: compactness (1)\n' \
                + 'Column 8: pressure (J/m3)\n' \
                + 'Column 9: coveringFactor (1)\n'

        info = np.column_stack((part['x'], part['y'], part['z'],
                            part['smoothLength'], part['sfr'],
                            part['Z'], part['compactness'],
                            part['pressure'], part['covering'])) # 9 params

        np.savetxt(self.workingDir + '/starforming_stars.txt', info, header=header)

        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                        & (np.abs(g.pos[:, 1]) < region)\
                        & (np.abs(g.pos[:, 2]) < region)\
                        & (starPart['age'] > ageThreshold))[0]
        size = mask.shape[0]
        sf = self.config['SEDFamily']
        print(f'{sf} particles:', size)
        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask] # in kpc
        part['z'] = g.pos[:, 2][mask] # in kpc
        part['smoothLength'] = g.smoothLength[mask] # smoothing length in kpc
        part['initialMass'] = g.mass[mask]
        part['Z'] = starPart['GFM_Metallicity'][mask]
        part['age'] = starPart['age'][mask]

        header = 'quenched particles\n' \
                    + '\n' \
                    + 'Column 1: x-coordinate (kpc)\n'\
                    + 'Column 2: y-coordinate (kpc)\n'\
                    + 'Column 3: z-coordinate (kpc)\n'\
                    + 'Column 4: smoothing length (kpc)\n'\
                    + 'Column 5: initial mass (Msun)\n'\
                    + 'Column 6: metallicity (1)\n'\
                    + 'Column 7: age (Myr)\n'
        
        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['smoothLength'], part['initialMass'],
                                part['Z'], part['age']))
        np.savetxt(self.workingDir + '/quenched_stars.txt', info, header=header)

        part = {}
        if self.config['includeDust']:
            try:
                fields = ['GFM_Metallicity', 'Coordinates', 'Masses',
                        'InternalEnergy', 'StarFormationRate', 'ElectronAbundance', 
                        'Density']
                gasPart = ill.snapshot.loadSubhalo(self.config['filePath'], 
                        self.snapnum, self.id, 'gas', fields)
                
                # smoothing length of gas cells is 2.5 times cell radius
                # Density in (10^10 Msun/h) / (ckpc/h)^3, Masses in 10^10 Msun/h
                # see https://arxiv.org/pdf/1707.03395
                volumes = gasPart['Masses'] / gasPart['Density']
                cellRadius = (3 * volumes / 4 / np.pi)**(1/3)
                gasPart['smoothLength'] = 2.5 * cellRadius # in ckpc/h
                
                gasPart['Temperature'] = u2temp(gasPart['InternalEnergy'], 
                                                gasPart['ElectronAbundance'])
            
                gas = Galaxy(gasPart['Coordinates'], gasPart['smoothLength'], gasPart['Masses'], 
                            self.pos, a=self.a, h=self.h)
                
                spatialmask = np.where((np.abs(gas.pos[:, 0]) < region)\
                            & (np.abs(gas.pos[:, 1]) < region)\
                            & (np.abs(gas.pos[:, 2]) < region))[0]
                
                if self.config['DISMModel'] == 'Camps_2016':
        
                    othermask = np.where((gasPart['StarFormationRate'] > 0) \
                                    | (gasPart['Temperature'] < np.float32(self.config['temperatureThreshold'])))[0]
                    mask = np.intersect1d(spatialmask, othermask)
                    size = mask.shape[0]
                
                elif self.config['DISMModel'] == 'Torrey_2012':
                    # Density from (10^10 Msun/h) / (ckpc/h)^3 to 10^10 h^2 Msun / kpc^3
                    density = gasPart['Density'] * self.a**-3
                    othermask = np.where(np.log10(gasPart['Temperature']) < (6 + 0.25 * np.log10(density)))[0]
                    mask = np.intersect1d(spatialmask, othermask)
                    size = mask.shape[0]
                                         
                print('Dust particles:', size)
                part['x'] = gas.pos[:, 0][mask] # in kpc
                part['y'] = gas.pos[:, 1][mask]
                part['z'] = gas.pos[:, 2][mask]
                part['smoothLength'] = gas.smoothLength[mask]
                part['mass'] = gas.mass[mask]
                part['Z'] = gasPart['GFM_Metallicity'][mask]
                part['Temperature'] = gasPart['Temperature'][mask]
            except:
                part['x'] = np.array([])
                part['y'] = np.array([])
                part['z'] = np.array([])
                part['smoothLength'] = np.array([])
                part['mass'] = np.array([])
                part['Z'] = np.array([])
                part['Temperature'] = np.array([])
        else:
            part['x'] = np.array([])
            part['y'] = np.array([])
            part['z'] = np.array([])
            part['smoothLength'] = np.array([])
            part['mass'] = np.array([])
            part['Z'] = np.array([])
            part['Temperature'] = np.array([])

        header = 'dust particles\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n' \
                + 'Column 2: y-coordinate (kpc)\n' \
                + 'Column 3: z-coordinate (kpc)\n' \
                + 'Column 4: smoothing length (kpc)\n' \
                + 'Column 5: mass (Msun)\n' \
                + 'Column 6: metallicity (1)\n' \
                + 'Column 7: temperature (K)\n'
                
        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['smoothLength'], part['mass'],
                                part['Z'], part['Temperature']))

        np.savetxt(self.workingDir + '/dusts.txt', info, header=header)

    def __get_properties(self, z):
        properties = {}
        properties['subhaloID'] = self.id
        properties['subhaloMass'] = np.log10(self.mass)
        properties['redshift'] = z
        distance = self.cosmology.luminosity_distance(z).value
        properties['lumiDis'] = distance
        if self.config['postProcessing']:
            surveys = split(self.config['surveys'])
            spatialResols = []
            for survey in surveys:
                filters = split(self.config[f'filters_{survey}'])
                numfilters = len(filters)
                if self.config[f'resolFromPix_{survey}']:
                    pixelScales = split(self.config[f'pixelScales_{survey}'], float)
                    pixelScales = extend(pixelScales, numfilters)
                    resol = [distance * 10**6 * np.deg2rad(ang / 3600) for ang in pixelScales]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                    
                else:
                    resol = split(self.config[f'resolution_{survey}'], float)
                    resol = extend(resol, numfilters)
                    pixelScales = [np.rad2deg(res / (distance * 10**6)) * 3600 for res in resol]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                
                properties[f'resolution_{survey}'] = resol
                properties[f'angleRes_{survey}'] = pixelScales
                properties[f'ps_in_sr_{survey}'] = ps_in_sr

                spatialResols += resol

                numExposure = split(self.config[f'numExposure_{survey}'], int)
                numExposure = extend(numExposure, numfilters)
                properties[f'numExposure_{survey}'] = numExposure

                exposureTime = split(self.config[f'exposureTime_{survey}'], float)
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

        numViews = np.int32(self.config['numViews'])
        properties['numViews'] = numViews
        randomViews = self.config['randomViews']
        if randomViews:
            inclinations = np.random.uniform(0, 180, numViews).tolist()
            azimuths = np.random.uniform(-360, 360, numViews).tolist()
        else:
            inclinations = split(self.config['inclinations'], float)
            azimuths = split(self.config['azimuths'], float)
        
        properties['inclinations'] = inclinations
        properties['azimuths'] = azimuths
        fovscale = np.float32(self.config['FoVboxLengthRatio'])
        fov = fovscale * self.boxLength * 10**3
        properties['FoV'] = fov
        
        return properties


    def __create_ski(self):

        print('Creating .ski file.')

        mode = self.config['simulationMode']

        ski_file = os.path.join(self.dataDir,  f'ski_templates/{mode}_template.ski')

        with open(ski_file, 'r') as file:
            data = file.read()

        z = np.float32(self.config['fixedRedshift'])

        # self.properties['fixed']
        data = data.replace('redshift="0.008"', f'redshift="{z}"')

        self.properties = self.__get_properties(z=z)
        self.properties['fixedRedshift'] = z

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
        
        numPackets = np.float32(self.config['numPackets'])
        data = data.replace('numPackets="1e7"', f'numPackets="{numPackets}"')
        
        minWavelength = np.float32(self.config['minWavelength']) * (1 + z)
        maxWavelength = np.float32(self.config['maxWavelength']) * (1 + z)
        
        data = data.replace('minWavelength="0.01 micron"', f'minWavelength="{minWavelength} micron"')
        data = data.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength} micron"')
        
        self.properties['minWavelength'] = minWavelength # micron
        self.properties['maxWavelength'] = maxWavelength # micron    

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
        
        wavelengthGrid = self.config['wavelengthGrid']
        grid_type = {'Linear': 'LinWavelengthGrid',
                'Log': 'LogWavelengthGrid'}
        grid_type = grid_type[wavelengthGrid]
        data = data.replace('LinWavelengthGrid', grid_type)
        
        numWavelengths = np.int32(self.config['numWavelengths'])
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
        fovscale = np.float32(self.config['FoVboxLengthRatio'])
        fov = fovscale * self.boxLength * 10**3 # in pc

        spatialResol = self.properties['baseRes']
        inclinations = self.properties['inclinations']
        azimuths = self.properties['azimuths']
        numViews = self.properties['numViews']

        numPixels = int(fov / spatialResol)

        recordComponents = 'true' if self.config['recordComponents'] else 'false'
        
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

        print('------estimate memory usage------')
        print(f'numViews: {numViews}')
        print(f'numSpatialPixels: {numPixels}')
        print(f'numWavelengthPixels: {numWavelengths}')
        factor = 7 if self.config['recordComponents'] else 1
        numPixels = np.float64(numPixels) # avoid overflow
        dataCubeSize = np.int64(numPixels ** 2 * numWavelengths * numViews)
        dataSize_in_GB = np.around(dataCubeSize * 8 * factor * 1e-9, 3)
        print(f'Estimated memory usage: {dataSize_in_GB} GB')

    def prepare(self):
        os.makedirs(self.workingDir, exist_ok=True)
        self.__get_particles()
        self.__create_ski()

    def runSKIRT(self):
        self.__run_skirt()
        prop = self.get_properties()
        dict_path = os.path.join(self.workingDir, 'properties.pkl')
        with open(dict_path, 'wb') as f:
            pickle.dump(prop, f)
        

    def __run_skirt(self):

        print('Running SKIRT')

        base = os.getcwd()
        os.chdir(self.workingDir)
        numThreads = int(self.config['numThreads'])
        if numThreads > 24:
            numThreads = 24
        command = f'skirt -t {numThreads} skirt.ski'
        process = subprocess.Popen(command, shell=True)
        flag = process.wait()

        if flag != 0:
            print('SKIRT exited with error.')
            sys.exit()

        os.chdir(base)
        return flag
    
    def get_properties(self):
        return self.properties