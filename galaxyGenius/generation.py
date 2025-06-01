import os
import subprocess
import sys
import json
from shutil import copyfile, move, rmtree
from types import NoneType
from typing import Union

class DataGeneration:
    
    def __init__(self, config: dict):
        
        self.workingDir = config['workingDir']
        self.properties = self.__get_properties()
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
        
    def __saveDataCube(self):
        numViews = int(self.properties['numViews'])
        
        dataCubeDir = f'dataCubes/Subhalo_{self.properties["subhaloID"]}'
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
            
    def __get_properties(self) -> dict:
        with open(self.workingDir + '/properties.json', 'r') as file:
            properties = json.load(file)
        return properties
    
    def __check_files(self):
        
        if not os.path.exists(os.path.join(self.workingDir, 'stars.txt')):
            print('stars.txt not found.')
            sys.exit()
            
        if not os.path.exists(os.path.join(self.workingDir, 'starforming_regions.txt')):
            print('starforming_regions.txt not found.')
            sys.exit()
            
        if not os.path.exists(os.path.join(self.workingDir, 'dusts.txt')):
            print('dusts.txt not found.')
            sys.exit()
        else:
            with open(os.path.join(self.workingDir, 'dusts.txt'), 'r') as file:
                lines = ''
                for _ in range(20):
                    lines += file.readline()
            
            if self.config['hydrodynamicSolver'] == 'smoothParticle':
                
                if not 'smoothing length' in lines:
                    print('Smoothing length must be provided for particle-based gas representation.')
                    sys.exit()
        
    def __run_skirt(self):
        print('Running SKIRT')
        print('Subhalo ID: ', self.properties['subhaloID'])
        
        base = os.getcwd()
        # if self.config['skirtPath'] == 'PATH':
        #     executable = 'skirt'
        # else:
        #     executable = self.config['skirtPath']
        
        os.chdir(self.workingDir)
        numThreads = int(self.config['numThreads'])
        if numThreads > 24:
            numThreads = 24
        command = f'skirt -t {numThreads} skirt.ski'
        
        result = subprocess.run(command, shell=True, check=True)
        
        run_flag = 0
        if result.returncode != 0:
            print('SKIRT exited with error.')
            run_flag = 1
            
        os.chdir(base)
        
        return run_flag
    
    def __exit(self):
        sys.exit()
    
    def runSKIRT(self):
        
        """
        This method runs the SKIRT radiative transfer simulation using the 
        configuration specified in the working directory. After execution, 
        it saves the resulting data cubes and cleans up the working directory.
        """
        
        self.__check_files()
        run_flag = self.__run_skirt()
        if run_flag == 1:
            return self.__exit()
        else:
            self.__saveDataCube()
            
            print('Cleaning up working directory')
            rmtree(self.workingDir)
        
            return 0