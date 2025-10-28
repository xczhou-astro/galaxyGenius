
from galaxyGenius.config import Configuration
from galaxyGenius.generation import DataGeneration
from galaxyGenius.preprocess import PreProcess
from galaxyGenius.postprocess import PostProcess

# initialize Configuration class
config = Configuration()

# add surveys
config.add_survey('HSC')
config.add_survey('CSST')

# get configuration dict
conf = config.get_config()

# initialize PreProcess class
preprocess = PreProcess(conf)

# get basic informations for subhalos 
subhalos = preprocess.get_subhalos()

# get subhaloIDs
subhaloIDs = subhalos['subhaloIDs']

for i, ID in enumerate(subhaloIDs):
    
    # specify subhalo by subhaloID
    preprocess.subhalo(ID)
    # prepare necessary files, including particle files and ski file
    preprocess.prepare()
    
    # initialize DataGeneration class
    dataGeneration = DataGeneration(config=conf)
    # run SKIRT
    dataGeneration.runSKIRT()
    
    # initialize PostProcess class for mock observations
    postprocess = PostProcess(subhaloID=ID)
    # perform mock observations
    postprocess.runPostprocess()
