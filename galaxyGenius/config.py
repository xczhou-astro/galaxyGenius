import os
from types import NoneType
from termcolor import colored
import numpy as np
from pathlib import Path
from tomlkit import parse, dumps
from typing import Union
from copy import deepcopy
from .utils import *
import subprocess

class Configuration:
    def __init__(self, surveys: Union[str, list[str], NoneType] = None):
        
        if surveys is not None:
            if isinstance(surveys, str):
                self.surveys = split(surveys)
            elif isinstance(surveys, list):
                self.surveys = surveys
        else:
            self.surveys = []
        
        if os.environ.get('GALAXYGENIUS_DATA_DIR') is not None:
            self.dataDir = os.environ.get('GALAXYGENIUS_DATA_DIR').split(':')[0]
        else:
            print('GALAXYGENIUS_DATA_DIR not set in environment variables. ' + 'Data directory falling to default path: ../Data')
            self.dataDir = '../Data'
            if not os.path.exists(self.dataDir):
                self.__issue('Data directory not found. Please set GALAXYGENIUS_DATA_DIR environment variable.')
                sys.exit()
            
        
        self.main_config_template = self.__read_config(
            os.path.join(self.dataDir, 'config/config_main.toml')
        )
        
        self.survey_config_template = self.__read_config(
            os.path.join(self.dataDir, 'config/config_survey.toml')
        )
    
    def __read_config(self, file: str) -> dict:
        
        with open(file, 'r') as f:
            config = parse(f.read())
        
        return config

    def __modify_main_config(self):
        if len(self.surveys) == 0:
            self.main_config_template['postProcessing'] = False
            self.main_config_template['surveys'] = []
            
            self.config['postProcessing'] = False
            self.config['surveys'] = []
        else:
            self.main_config_template['postProcessing'] = True
            self.main_config_template['surveys'] = self.surveys
            
            self.config['postProcessing'] = True
            self.config['surveys'] = self.surveys
            
            
    def init(self):
        
        """
        Initialize configuration. Only use by init.py
        """
        
        self.__update_config()
        self.check_config()
        self.save_config()
        
    def add_survey(self, surveys: Union[str, list[str], NoneType] = None):
        
        '''
        Add surveys. Call ``get_config()`` to get the updated configs.

        Parameters
        ----------
        surveys : str or list[str] or None
            Survey names to add. If None, no surveys will be added.
        '''
        
        if surveys is not None:
            if isinstance(surveys, str):
                surveys_add = split(surveys)
            else:
                surveys_add = surveys
        else:
            print('Please provide survey names.')
            surveys_add = []
            
        added_surveys = []
        for survey in surveys_add:
            if survey in self.surveys:
                print(f'{survey} already exists.')
            else:
                print(f'{survey} added.')
                added_surveys.append(survey)
        
        self.surveys += added_surveys
    
    def remove_survey(self, surveys: Union[str, list[str], NoneType] = None):
        
        '''
        Remove surveys. Call ``get_config()`` to get the updated configs.

        Parameters
        ----------
        surveys : str or list[str] or None
            Survey names to remove. If None, no surveys will be removed.
        '''
        
        if surveys is not None:
            if isinstance(surveys, str):
                surveys_remove = split(surveys)
            else:
                surveys_remove = surveys
        else:
            print('Please provide survey names.')
            surveys_remove = []
        
        for survey in surveys_remove:
            if survey in self.surveys:
                print(f'{survey} removed.')
                if os.path.exists(f'config_{survey}.toml'):
                    os.remove(f'config_{survey}.toml')
                self.surveys.remove(survey)
            else:
                print(f'{survey} does not exist.')
    
    def __load_config(self):
        self.config = self.__read_config('config.toml')
        
        self.survey_configs = {}
        survey_config_path = Path('.')
        survey_config_files = list(survey_config_path.glob('config_*.toml'))
        detected_surveys = [
            file.name.split('_')[1].split('.')[0] for file in survey_config_files
        ]
        
        
        for survey in detected_surveys:
            if survey not in self.surveys:
                self.surveys += [survey]
            
            config = self.__read_config(f'config_{survey}.toml')
            
            self.survey_configs[survey] = config
            
        
        
            
        self.__modify_main_config()
    
    def __create_survey_config(self, survey: str) -> dict:
        
        survey_config = self.survey_config_template.copy()
        adapted_survey_config = self.__adapt_survey(survey, survey_config)
        return adapted_survey_config
    
    def __update_config(self):
        
        if len(self.surveys) == 0:
            if os.path.exists('config.toml'):
                self.__load_config()
            else:
                self.config = self.main_config_template
                self.__modify_main_config()
        
        else:
            if os.path.exists('config.toml'):
                self.__load_config()
                
                self.survey_configs = {}
                for survey in self.surveys:
                    if not os.path.exists(f'config_{survey}.toml'):
                        survey_config = self.__create_survey_config(survey)
                        self.survey_configs[survey] = survey_config
                    else:
                        survey_config = self.__read_config(f'config_{survey}.toml')
                        self.survey_configs[survey] = survey_config
                        
            else:
                self.config = self.main_config_template
                self.__modify_main_config()
                self.survey_configs = {}
                for survey in self.surveys:
                    survey_config = self.__create_survey_config(survey)
                    self.survey_configs[survey] = survey_config
                        
    def __config_to_dict(self):
        
        config = self.config
        
        if hasattr(self, 'survey_configs'):
            survey_configs = deepcopy(self.survey_configs)
            
            survey_config = {}
            for survey in self.surveys:
                for key in list(survey_configs[survey].keys()):
                    key_suffix = key + f'_{survey}'
                    survey_config[key_suffix] = survey_configs[survey][key]
                
            config = config | survey_config
        
        return config
              
    
    def save_config(self, conf: Union[dict, NoneType] = None):
        
        """
        Manually save the current configuration to the specified TOML files.

        Parameters:
        conf (dict, optional): A dictionary containing configuration 
        values.
        """
        
        if conf is not None:
        
            for key in list(self.config.keys()):
                if conf[key] != self.config[key]:
                    self.config[key] = conf[key]
                    
            self.config['dataDir'] = self.dataDir
            
                    
            for survey in self.surveys:
                for key in list(self.survey_configs[survey]):
                    if conf[key + f'_{survey}'] != self.survey_configs[survey][key]:
                        self.survey_configs[survey][key] = conf[key + f'_{survey}']
        
        with open('config.toml', 'w') as f:
            f.write(dumps(self.config))
            
        for survey in self.surveys:
            with open(f'config_{survey}.toml', 'w') as f:
                f.write(dumps(self.survey_configs[survey]))
        
                
    def get_config(self):
        
        """
        Update and return the current configuration.

        Returns:
            dict: The updated configuration dictionary.
        """

        self.__update_config()
        self.check_config()
        config = self.__config_to_dict()
        config['dataDir'] = self.dataDir
        
        return dict(config)
    
    def __adapt_survey(self, survey: str, survey_config: dict) -> dict:
        if survey == 'HST':
            adapted_survey_config = self.__for_HST(survey_config)
        elif survey == 'JWST':
            adapted_survey_config = self.__for_JWST(survey_config)
        elif survey == 'Euclid':
            adapted_survey_config = self.__for_Euclid(survey_config)
        elif survey == 'Roman':
            adapted_survey_config = self.__for_Roman(survey_config)
        elif survey == 'HSC':
            adapted_survey_config = self.__for_HSC(survey_config)
        elif survey == 'CSST':
            adapted_survey_config = survey_config
        else:
            adapted_survey_config = survey_config
            print(f'Please provide your own Throughputs and PSFs of filters for {survey}')
        
        return adapted_survey_config
    
    def __for_HST(self, survey_config: dict) -> dict:
        
        survey_config['filters'] = ["UV_F275W", "UV_F390W", "UV_F814W", "IR_F105W", "IR_F140W"]
        survey_config['resolFromPix'] = True
        survey_config['pixelScales'] = [0.04, 0.04, 0.04, 0.13, 0.13]
        survey_config['numExposure'] = 1
        survey_config['exposureTime'] = 600
        survey_config['aperture'] = 2.4
        survey_config['RGBImg'] = False
        survey_config['RGBFilters'] = ''
        survey_config['displayFilter'] = "UV_F814W"
        survey_config['skyBkg'] = [0.0004, 0.0181, 0.0562, 1.9397, 2.2485]
        survey_config['darkCurrent'] = [0.0008, 0.0008, 0.0008, 0.048, 0.048]
        survey_config['readOut'] = [3.1, 3.1, 3.1, 21, 21]
        
        return survey_config
    
    def __for_HSC(self, survey_config: dict) -> dict:
        
        survey_config['filters'] = ["g", "r", "i", "z", "y"]
        survey_config['resolFromPix'] = True
        survey_config['pixelScales'] = 0.168
        survey_config['numExposure'] = 1
        survey_config['exposureTime'] = 600
        survey_config['aperture'] = 8.5
        survey_config['RGBImg'] = True
        survey_config['RGBFilters'] = ["g", "r", "z"]
        survey_config['displayFilter'] = "r"
        survey_config['skyBkg'] = [49.983, 79.8137, 91.8011, 107.205, 198.5708]
        survey_config['darkCurrent'] = 0
        survey_config['readOut'] = 4.5
        
        return survey_config

    
    def __for_JWST(self, survey_config: dict) -> dict:
        
        survey_config['filters'] = ["F356W", "F444W", "F070W", "F150W", "F200W", "F182M"]
        survey_config['resolFromPix'] = True
        survey_config['pixelScales'] = [0.063, 0.063, 0.031, 0.031, 0.031, 0.031]
        survey_config['numExposure'] = 1
        survey_config['exposureTime'] = 600
        survey_config['aperture'] = 6.5
        survey_config['RGBImg'] = False
        survey_config['RGBFilters'] = ''
        survey_config['displayFilter'] = 'F200W'
        survey_config['skyBkg'] = [1.1707, 3.2627, 0.3307, 0.4721, 0.4326, 0.2502]
        survey_config['darkCurrent'] = [0.0342, 0.0342, 0.0019, 0.0019, 0.0019, 0.0019]
        survey_config['readOut'] = [13.25, 13.25, 15.77, 15.77, 15.77, 15.77]
        
        self.config['minWavelength'] = 0.5
        self.config['maxWavelength'] = 5
        self.config['simulationMode'] = 'DustEmission'
        self.config['includeDust'] = True
        
        return survey_config
    
    def __for_Euclid(self, survey_config: dict) -> dict:

        survey_config['filters'] = ["VIS", "Y", "J", "H"]
        survey_config['resolFromPix'] = True
        survey_config['pixelScales'] = [0.1, 0.3, 0.3, 0.3]
        survey_config['numExposure'] = 1
        survey_config['exposureTime'] = 600
        survey_config['aperture'] = 1.2
        survey_config['PSFFromFile'] = False
        survey_config['PSFFWHM'] = [0.204, 0.493, 0.515, 0.553]
        survey_config['RGBImg'] = False
        survey_config['RGBFilters'] = ''
        survey_config['displayFilter'] = 'VIS'
        survey_config['skyBkg'] = [0.7162, 3.6863, 4.4556, 3.8228]
        survey_config['darkCurrent'] = [0.005, 0.02, 0.02, 0.02]
        survey_config['readOut'] = [2.2, 6, 6, 6]
        
        self.config['minWavelength'] = 0.5
        self.config['maxWavelength'] = 5
        self.config['simulationMode'] = 'DustEmission'
        self.config['includeDust'] = True
        
        return survey_config
        
    def __for_Roman(self, survey_config: dict) -> dict:
        
        survey_config['filters'] = ["F062", "F087", "F106", "F129", "F158", "F184", "F213", "F146"]
        survey_config['resolFromPix'] = True
        survey_config['pixelScales'] = 0.11
        survey_config['numExposure'] = 1
        survey_config['exposureTime'] = 600
        survey_config['aperture'] = 2.36
        survey_config['RGBImg'] = False
        survey_config['RGBFilters'] = ''
        survey_config['displayFilter'] = 'F106'
        survey_config['skyBkg'] = [7.1684, 6.6275, 7.0874, 7.1373, 7.0349, 4.1791, 3.515, 21.3173]
        survey_config['darkCurrent'] = 0.1
        survey_config['readOut'] = 20
        
        self.config['minWavelength'] = 0.5
        self.config['maxWavelength'] = 5
        self.config['simulationMode'] = 'DustEmission'
        self.config['includeDust'] = True
        
        return survey_config
        
    def __issue(self, message: str):
        print(colored(message, 'red'))

    def __exist(self, key: str, threshold: Union[float, NoneType] = None) -> bool:
        if key in self.all_config:
            if threshold is not None:
                values = self.all_config[f'{key}']
                if isinstance(values, (int, float)):
                    values = [values]
                values = np.array([np.float32(val) for val in values])
                check = (values > threshold).all()
                if check == False:
                    self.__issue(f'{key} value must > {threshold}.')
                    self.flag_count += 1
                    return False
                else:
                    return True
            else:
                return True
        
        else:
            self.__issue(f'{key} is not provided.')
            self.flag_count += 1
            return False
        
    def __exist_return(self, key: str, threshold: float = None) -> bool:
        value = self.__exist(key, threshold=threshold)
        if value:
            return self.all_config[key]
        else:
            return False
        
    def __match(self, base_key: str, key: str, homo: bool = False) -> bool:
        
        if self.__exist(key):
            
            base = self.all_config[base_key]
            values = self.all_config[key]
            
            if isinstance(values, (int, float)):
                values = [values]
            
            if len(values) == len(base):
                return True
            elif (len(values) == 1) & (homo):
                return True
            else:
                self.__issue(f'number of {key} inconsistent with number of {base_key}.')
                self.flag_count += 1
                return False
    
    def __logic(self, base_key: str, key: str, relation: str = 'match') -> bool:
        if relation == 'match':
            if self.all_config[base_key] != self.all_config[key]:
                self.__issue(f'{key} must be {self.all_config[base_key]} if {base_key} is {self.all_config[base_key]}.')
                self.flag_count += 1
                return False
        elif relation == 'diff':
            if self.all_config[base_key] == self.all_config[key]:
                self.__issue(f'{key} must be {not self.all_config[base_key]} if {base_key} is {self.all_config[base_key]}.')
                self.flag_count += 1
                return False
        else:
            return True

    def check_config(self) -> int:
        
        """
        Check the configuration for any inconsistencies and issues.

        Returns:
            int: The number of issues found in the configuration.
        """
        
        self.all_config = self.__config_to_dict()
        
        self.flag_count = 0
        
        # skirtPATH = self.__exist_return('skirtPath')
        # if skirtPATH == 'PATH':
        #     sys_paths = os.environ['PATH'].split(':')
        #     for path in sys_paths:
        #         if 'SKIRT' in path:
        #             skirtPATH = path
        #             break
                
        #     if skirtPATH == 'PATH':
        #         self.__issue('SKIRT not found in PATH. Please edit ~/.bashrc or provide the absolute path of SKIRT.')
        #         self.flag_count += 1

        # executable = os.path.join(skirtPATH, 'skirt')
        # if not os.path.exists(executable):
        #     self.__issue('SKIRT executable not found. SKIRT is not installed correctly.')
        #     self.flag_count += 1
            
        simName = self.__exist_return('simulation')
        
        requests = self.__exist_return('requests')
        
        solver = self.__exist_return('hydrodynamicSolver')
        if solver not in ['VoronoiMesh', 'smoothParticle']:
            self.__issue('hydrodynamicSolver unrecognized.')
            self.flag_count += 1
        
        if 'TNG' in simName:
            if requests:
                self.__exist('apiKey')
            else:
                self.__exist('TNGPath')
                self.__exist('postprocessingPath')
        
                tngpath = self.__exist_return('TNGPath')
                postpath = self.__exist_return('postprocessingPath')
            
                # if os.path.exists(tngpath) and os.path.exists(postpath):
                #     if tngpath.split('/')[:-1] != postpath.split('/')[:-1]:
                #         self.__issue('TNGPath and postprocessingPath should be in the same directory.')
                #         self.flag_count += 1
                # else:
                #     if not os.path.exists(tngpath):
                #         self.__issue('TNGPath not found.')
                #         self.flag_count += 1
                #     if not os.path.exists(postpath):
                #         self.__issue('postprocessingPath not found.')
                #         self.flag_count += 1
        
        self.__exist('workingDir')
        
        if self.__exist('inLocal'):
            self.__exist('viewDistance', 0)
            
        self.__exist('snapRedshift')
        
        simulationMode = self.__exist_return('simulationMode')
        if simulationMode:
            if not simulationMode in ['NoMedium', 'DustEmission']:
                self.__issue('simulationMode unrecognized.')
                self.flag_count += 1

                if simulationMode == 'DustEmission':
                    includeDust = self.__exist_return('includeDust')
                    if not includeDust:
                        self.__issue('includeDust should be True if simulationMode is DustEmission.')
                        self.flag_count += 1
                    
                    dustEmissionType = self.__exist_return('dustEmissionType')
                    if dustEmissionType:
                        if not dustEmissionType in ['Equilibrium', 'Stochastic']:
                            self.__issue('dustEmissionType unrecognized.')
                            self.flag_count += 1
            
        
        dustModel = self.__exist_return('dustModel')
        if dustModel:
            if (dustModel != 'ZubkoDustMix') & (dustModel != 'DraineLiDustMix') & (dustModel != 'ThemisDustMix'):
                self.__issue('dustModel unrecognized.')
                self.flag_count += 1
                
        self.__exist('minWavelength')
        self.__exist('maxBoxLength')
        self.__exist('boxLengthScale')
                
        wavelengthGrid = self.__exist_return('wavelengthGrid')
        if wavelengthGrid:
            if (wavelengthGrid != 'Linear') & (wavelengthGrid != 'Log'):
                self.__issue('wavelengthGrid unrecognized.')
                self.flag_count += 1
                
        self.__exist('numWavelengths')
        self.__exist('minLevel')
        self.__exist('maxLevel')
        self.__exist('numPackets')
        
        if self.__exist_return('faceAndEdge'):
            pass
        else:
            if self.__exist('numViews'):
                if not self.__exist_return('randomViews'):
                    if self.__exist('inclinations') & self.__exist('azimuths'):
                    
                        numViews = self.all_config['numViews']
                        inclinations = self.all_config['inclinations']
                        azimuths = self.all_config['azimuths']
                        
                        if (len(inclinations) != numViews) | (len(azimuths) != numViews):
                            self.__issue('number of inclinations or azimuth inconsistent with numViews.')
                            self.flag_count += 1
                        
                        inc_in_cond = [True if 0 <= inc <= 180 else False for inc in inclinations]
                        azi_in_cond = [True if -360 <= azi <= 360 else False for azi in azimuths]

                        if not all(inc_in_cond):
                            self.__issue('inclinations must be in 0 to 180 degree.')
                            self.flag_count += 1

                        if not all(azi_in_cond):
                            self.__issue('azimuths must be in -360 to 360 degree.')
                            self.flag_count += 1
                    
        SEDFamily = self.__exist_return('SEDFamily')
        if SEDFamily:
            if (SEDFamily != 'BC03') & (SEDFamily != 'FSPS'):
                self.__issue('SEDFamily unrecognized.')
                self.flag_count += 1
                
            elif SEDFamily == 'BC03':
                initialMassFunction = self.__exist_return('initialMassFunction')
                if initialMassFunction:
                    if (initialMassFunction != 'Chabrier') & (initialMassFunction != 'Salpeter'):
                        self.__issue('initialMassFunction unrecognized for BC03 SEDFamily.')
                        self.flag_count += 1
                    
            elif SEDFamily == 'FSPS':
                initialMassFunction = self.__exist_return('initialMassFunction')
                if initialMassFunction:
                    if  (initialMassFunction != 'Chabrier') & (initialMassFunction != 'Salpeter') & (initialMassFunction != 'Kroupa'):
                        self.__issue('initialMassFunction unrecognized for FSPS SEDFamily.')
                        self.flag_count += 1
                        
        imageUnit = self.__exist_return('imageUnit')
        if (imageUnit != 'electron') & (imageUnit != 'flux'):
            self.__issue('imageUnit unrecognized.')
            self.flag_count += 1
                        
        self.__exist('minStellarMass')
        self.__exist('maxStellarMass')
                        
        if self.__exist_return('displaySED'):
            self.__exist('displaySEDxlogscale')
            
        self.__exist('snapNum')
        
        self.__exist('spatialResol')
        
        if self.__exist_return('postProcessing'):
            surveys = self.__exist_return('surveys')
            pivots = []
            if surveys:
                for survey in surveys:
                    filters = self.__exist_return(f'filters_{survey}')
                    if filters:
                        for filter in filters:
                            if not os.path.exists(f'{self.dataDir}/filters/{survey}/{filter}.fil'):
                                self.__issue(f'Throughput file {survey}/{filter} not found! Please specify correct filters!')
                                self.flag_count += 1
                            else:
                                pivots.append(calc_pivot(self.dataDir, survey, filter))
                                
                        if self.__exist(f'numExposure_{survey}'):
                            self.__match(f'filters_{survey}', f'numExposure_{survey}', True)
                        
                        if self.__exist(f'exposureTime_{survey}'):
                            self.__match(f'filters_{survey}', f'exposureTime_{survey}', True)
                        
                        self.__exist(f'aperture_{survey}')
                        
                        if self.__exist_return(f'resolFromPix_{survey}'):
                            pixelScales = self.__exist_return(f'pixelScales_{survey}', 0)
                            if pixelScales:
                                self.__match(f'filters_{survey}', f'pixelScales_{survey}', True)
                            
                            if self.__exist_return(f'includePSF_{survey}'):
                                PSFFromFile = self.__exist_return(f'PSFFromFile_{survey}')
                                if not PSFFromFile:
                                    PSFFWHM = self.__exist_return(f'PSFFWHM_{survey}')
                                    if PSFFWHM:
                                        self.__match(f'filters_{survey}', f'PSFFWHM_{survey}', True)
                                else:
                                    if os.path.exists(f'{self.dataDir}/PSFs/{survey}'):
                                        psffiles = os.listdir(f'{self.dataDir}/PSFs/{survey}')
                                        psffilters = [name.split('.')[0] for name in psffiles]
                                        for filter in filters:
                                            if not filter in psffilters:
                                                self.__issue(f'PSF file {survey}/{filter} not found.')
                                                self.__issue(f'Please add PSF file in Data/PSFs.')
                                                self.flag_count += 1
                                    else:
                                        self.__issue(f'PSF folder of {survey} not found.')
                                        self.flag_count += 1
                            
                            if self.__exist_return(f'includeBkg_{survey}'):
                                
                                if self.__exist_return(f'noiseType_{survey}') not in ['instrument', 'limitingMagnitude']:
                                    self.__issue('noiseType unrecognized.')
                                    self.flag_count += 1
                                    
                                else:
                                    if self.__exist_return(f'noiseType_{survey}') == 'instrument':
                                
                                        self.__match(f'filters_{survey}', f'skyBkg_{survey}', True)
                                        self.__match(f'filters_{survey}', f'darkCurrent_{survey}', True)
                                        self.__match(f'filters_{survey}', f'readOut_{survey}', True)
                                        
                                    elif self.__exist_return(f'noiseType_{survey}') == 'limitingMagnitude':
                                        self.__match(f'filters_{survey}', f'limitMag_{survey}', True)
                                        self.__match(f'filters_{survey}', f'limitSNR_{survey}', True)
                                        self.__match(f'filters_{survey}', f'limitAperture_{survey}', True)
                                        self.__match(f'filters_{survey}', f'zeroPoint_{survey}', True)

                                        
                            
                            if self.__exist_return(f'imgDisplay_{survey}'):
                                if self.__exist_return(f'RGBImg_{survey}'):
                                    RGBFilters = self.__exist_return(f'RGBFilters_{survey}')
                                    if not all([True if fil in filters else False for fil in RGBFilters]):
                                        self.__issue(f'RBGFilters not in {survey} filter set.')
                                        self.flag_count += 1
                                else:
                                    displayFilter = self.__exist_return(f'displayFilter_{survey}')
                                    if displayFilter:
                                        if not displayFilter in filters:
                                            self.__issue(f'displayFilter not in {survey} filter set.')
                                            self.flag_count += 1
                        else:
                            if self.__exist_return(f'includePSF_{survey}'):
                                self.__issue(f'includePSF_{survey} can only be False if resolFromPix_{survey} is False.')
                                self.flag_count += 1
                            
                            if self.__exist_return(f'includeBkg_{survey}'):
                                self.__issue(f'includeBkg_{survey} can only be False if resolFromPix_{survey} is False.')
                                self.flag_count += 1
                            
                            resolution = self.__exist_return(f'resolution_{survey}', 0)
                            if resolution:
                                self.__match(f'filters_{survey}', f'resolution_{survey}', True)
                                
                            if self.__exist_return(f'imgDisplay_{survey}'):
                                if self.__exist_return(f'RGBImg_{survey}'):
                                    RGBFilters = self.__exist_return(f'RGBFilters_{survey}')
                                    if not all([True if fil in filters else False for fil in RGBFilters]):
                                        self.__issue(f'RBGFilters not in {survey} filter set.')
                                        self.flag_count += 1
                                else:
                                    displayFilter = self.__exist_return(f'displayFilter_{survey}')
                                    if displayFilter:
                                        if not displayFilter in filters:
                                            self.__issue(f'displayFilter not in {survey} filter set.')
                                            self.flag_count += 1
            
                             
            if (len(pivots) > 0) & (self.__exist('maxWavelength')) & (self.__exist('minWavelength')):
                max_pivot = np.max(pivots)
                min_pivot = np.min(pivots)

                redshift = self.__exist_return('viewRedshift')
                if np.isclose(redshift, 0, atol=0.005):
                    redshift = 0.005
                redshift = np.float32(redshift)
                
                
                maxWavelength = np.float32(self.all_config['maxWavelength']) * 10**4 * (1 + redshift)
                minWaveLength = np.float32(self.all_config['minWavelength']) * 10**4 * (1 + redshift)
                
                if (max_pivot > maxWavelength) | (min_pivot < minWaveLength):
                    self.__issue('Please make sure the redshiftted wavelength range covers all filters.')
                    self.flag_count += 1
                
                if (max_pivot > 2 * 10**4) & (self.all_config['simulationMode'] != 'DustEmission'):
                    self.__issue('Filters reaching infrared, simulationMode should be DustEmission.')
            
        
        self.__exist('outputSEDOnly')
        self.__exist('fieldOfView')
        
        if self.__exist('numThreads'):
            if self.all_config['numThreads'] > 24:
                self.__issue('No speed improvement with numThreads larger than, 24, falling back to 24.')

        self.__exist('recordComponents')
        self.__exist('ageThreshold')
        self.__exist('logCompactnessMean')
        self.__exist('logCompactnessStd')
        self.__exist('massFraction')
        covering = self.__exist_return('constantCoveringFactor')
        if covering:
            self.__exist('coveringFactor')
        else:
            self.__exist('PDRClearingTimescale')

        self.__exist('logPressure')
        self.__exist('temperatureThreshold')
        DISMModel = self.__exist_return('DISMModel')
        if DISMModel:
            if (DISMModel != 'Camps_2016') & (DISMModel != 'Torrey_2012'):
                self.__issue('DISMModel unrecognized.')
                self.flag_count += 1
        
        self.__exist('numSilicateSizes')
        self.__exist('numGraphiteSizes')
        self.__exist('numPAHSizes')
        self.__exist('numHydrocarbonSizes')

        if self.flag_count == 0:
            print(colored('No conflicts in configurations. ', 'green') + '\U0001F44D')
        else:
            print(colored('Conflicts exist in configurations. ', 'red') + '\U0001F914')
            sys.exit()
        
        return self.flag_count