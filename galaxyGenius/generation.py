import os
import subprocess
import sys
import json
from shutil import copyfile, move, rmtree
from typing import Union

from .utils import setup_logging, read_properties
import logging

class DataGeneration:
    
    def __init__(self, config: dict):
        
        self.workingDir = config['workingDir']
        self.logger = setup_logging(os.path.join(self.workingDir, 'galaxyGenius.log'))
        self.logger.info(f'Initializing DataGeneration class.')
        self.properties = read_properties(self.workingDir)
        self.config = self.__read_configs()
        
    def __read_configs(self):
        
        with open(self.workingDir + '/config.json', 'r') as file:
            config = json.load(file)
            
        return config
            
    def __save_basics(self, directory: str):
        copyfile(os.path.join(self.workingDir, 'skirt_parameters.xml'),
                os.path.join(directory, 'skirt_parameters.xml'))
        copyfile(os.path.join(self.workingDir, 'stars.txt'),
                os.path.join(directory, 'stars.txt'))
        copyfile(os.path.join(self.workingDir, 'starforming_regions.txt'),
                os.path.join(directory, 'starforming_regions.txt'))
        copyfile(os.path.join(self.workingDir, 'dusts.txt'),
                os.path.join(directory, 'dusts.txt'))
        copyfile(os.path.join(self.workingDir, 'properties.json'),
                os.path.join(directory ,'properties.json'))
        copyfile(os.path.join(self.workingDir, 'skirt_log.txt'),
                os.path.join(directory, 'skirt_log.txt'))
        copyfile(os.path.join(self.workingDir, 'config.json'),
                os.path.join(directory, 'config.json'))
        copyfile(os.path.join(self.workingDir, 'galaxyGenius.log'),
                os.path.join(directory, 'galaxyGenius.log'))
        
    def __saveDataCube(self):
        numViews = int(self.properties['numViews'])
        
        dataCubeDir = f'dataCubes/Subhalo_{self.properties["subhaloID"]}'
        self.dataCubeDir = dataCubeDir
        
        os.makedirs(dataCubeDir, exist_ok=True)
        self.__save_basics(dataCubeDir)
        for i in range(numViews):
            
            if self.config['outputSEDOnly']:
                move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_sed.dat'), 
                    os.path.join(dataCubeDir, f'skirt_view_{i:02d}_sed.dat'))
            else:
                move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_total.fits'),
                    os.path.join(dataCubeDir, f'skirt_view_{i:02d}_total.fits'))
                move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_sed.dat'), 
                    os.path.join(dataCubeDir, f'skirt_view_{i:02d}_sed.dat'))
            
    # def __get_properties(self) -> dict:
    #     with open(self.workingDir + '/properties.json', 'r') as file:
    #         properties = json.load(file)
            
    #     for key in properties.keys():
    #         if 'value' in properties[key] and 'unit' in properties[key]:
    #             properties[key] = u.Quantity(properties[key]['value']) * properties[key]['unit']
                
    #     return properties
    
    def __check_files(self):
        
        if not os.path.exists(os.path.join(self.workingDir, 'stars.txt')):
            raise FileNotFoundError('stars.txt not found.')
            
        if not os.path.exists(os.path.join(self.workingDir, 'starforming_regions.txt')):
            raise FileNotFoundError('starforming_regions.txt not found.')
            
        if not os.path.exists(os.path.join(self.workingDir, 'dusts.txt')):
            raise FileNotFoundError('dusts.txt not found.')
        else:
            with open(os.path.join(self.workingDir, 'dusts.txt'), 'r') as file:
                lines = ''
                for _ in range(20):
                    lines += file.readline()
            
            if self.config['hydrodynamicSolver'] == 'smoothParticle':
                
                if not 'smoothing length' in lines:
                    raise ValueError('Smoothing length must be provided for particle-based gas representation.')
        
    def __run_skirt(self, skirtPath: Union[str, None]=None):
        self.logger.info('Running SKIRT')
        self.logger.info(f'Subhalo ID: {self.properties["subhaloID"]}')
        
        base = os.getcwd()
        
        if skirtPath is None:
            executable = 'skirt'
        else:
            executable = skirtPath
        
        os.chdir(self.workingDir)
        numThreads = int(self.config['numThreads'])
        if numThreads > 24:
            numThreads = 24
        command = f'{executable} -t {numThreads} skirt.ski'
        
        result = subprocess.run(command, shell=True, check=True)
        
        run_flag = 0
        if result.returncode != 0:
            self.logger.error('SKIRT exited with error.')
            run_flag = 1
            
        os.chdir(base)
        
        return run_flag
    
    def __exit(self):
        sys.exit()
    
    def runSKIRT(self, skirtPath: Union[str, None]=None):
        
        """Run SKIRT radiative transfer simulation.

        This method executes the SKIRT radiative transfer simulation using the configured parameters.
        It checks required files, runs SKIRT, saves the output data cube, and cleans up temporary files.

        Args:
            skirtPath (str, optional): Path to SKIRT executable. If None, assumes 'skirt' is in PATH.

        Returns:
            int: 0 if successful, exits with error otherwise.
        """
        
        self.__check_files()
        run_flag = self.__run_skirt(skirtPath)
        if run_flag == 1:
            return self.__exit()
        else:
            self.logger.info('Cleaning up working directory')
            
            root_logger = logging.getLogger()
            for handler in root_logger.handlers[:]:
                handler.close()
                root_logger.removeHandler(handler)
            
            self.__saveDataCube()
            rmtree(self.workingDir, ignore_errors=True)
            
            new_log_path = os.path.join(self.dataCubeDir, 'galaxyGenius.log')
            self.logger = setup_logging(new_log_path, force_reconfigure=True)
        
            return 0