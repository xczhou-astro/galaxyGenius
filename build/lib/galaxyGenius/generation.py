import os
import subprocess
import sys
import json
from shutil import copyfile, move, rmtree

class DataGeneration:
    
    def __init__(self, config: dict):
        
        self.config = config
        self.workingDir = config['workingDir']
        self.properties = self.__get_properties()
        
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
            move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_total.fits'),
                os.path.join(dataCubeDir, f'skirt_view_{i:02d}_total.fits'))
            move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_sed.dat'), 
                os.path.join(dataCubeDir, f'skirt_view_{i:02d}_sed.dat'))
            
    def __get_properties(self) -> dict:
        with open(self.workingDir + '/properties.json', 'r') as file:
            properties = json.load(file)
        return properties
        
    def __run_skirt(self):
        print('Running SKIRT')

        base = os.getcwd()
        os.chdir(self.workingDir)
        numThreads = int(self.config['numThreads'])
        if numThreads > 24:
            numThreads = 24
        command = f'skirt -t {numThreads} skirt.ski'
        try:
            result = subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError:
            print('SKIRT exited with error.')
            os.chdir(base)
            sys.exit()

        os.chdir(base)
    
    def runSKIRT(self):
        self.__run_skirt()
        self.__saveDataCube()
        
        print('Cleaning up working directory')
        rmtree(self.workingDir)
    
    
        