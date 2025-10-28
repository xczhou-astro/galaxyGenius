import os

# Set the GALAXYGENIUS_DATA_DIR environment variable
os.environ['GALAXYGENIUS_DATA_DIR'] = 'your/galaxygenius/Data'

import json
from astropy.cosmology import Planck15

from galaxyGenius.config import Configuration
from galaxyGenius.generation import DataGeneration
from galaxyGenius.preprocess import PreProcess
from galaxyGenius.postprocess import PostProcess
from galaxyGenius.utils import Units

config = Configuration()
config.add_survey('CSST,HSC') # add two surveys, CSST and HSC

conf = config.get_config()
config.save_config(conf)

# The units singleton for specific hydrodymical simulation
# If not defined, the default unit convention is for Planck15 at snapshot redshift
units = Units(cosmology=Planck15, 
              snapRedshift=conf['snapRedshift'])

# Initialize PreProcess
# The calculation in PreProcess will use the unit conventions in units singleton
preprocess = PreProcess(conf)

subhaloParticleFile = 'data/TNG_100_snap_94_subhalo_31.h5'

with open('data/TNG_100_snap_94_subhalo_31.json', 'r') as f:
    data = json.load(f)

stellarMass = data['mass_stars'] * units.mass # units.mass is 10^10 * Msun / h
halfStellarMassRadius = data['halfmassrad_stars'] * units.distance # units.distance is ckpc / h

subhaloInfo = {
    'SubhaloID': 31,
    'stellarMass': stellarMass,
    'halfStellarMassRadius': halfStellarMassRadius
}

pos_x = data['pos_x']
pos_y = data['pos_y']
pos_z = data['pos_z']
centerPosition = [pos_x, pos_y, pos_z] * units.position # units.position is ckpc / h

preprocess.inputSubhaloParticleFile(subhaloParticleFile, subhaloInfo, centerPosition)

# specify views to be generated
arguments = {'faceAndEdge': False,
             'randomViews': False,
             'inclinations': [0], 
             'azimuths': [0]}
preprocess.prepare(arguments)

dataGeneration = DataGeneration(conf)
dataGeneration.runSKIRT()

postprocess = PostProcess(subhaloID=31)
postprocess.runPostprocess()