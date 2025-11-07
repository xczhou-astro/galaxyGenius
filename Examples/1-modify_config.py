
from galaxyGenius.config import Configuration
from galaxyGenius.generation import DataGeneration
from galaxyGenius.preprocess import PreProcess
from galaxyGenius.postprocess import PostProcess

config = Configuration()

config.add_survey('HSC')
config.add_survey('CSST')

conf = config.get_config()

# Modify config
conf['requests'] = True # Use Web-based API of TNG simulation
conf['apiKey'] = 'your_api_key' # specify your own key
conf['simulationMode'] = 'DustEmission'
conf['includeDust'] = True
conf['numPackets'] = 1e6
conf['minStellarMass'] = 1e11

# Modify survey config
conf['includeBkg_CSST'] = False
conf['includePSF_CSST'] = False

preprocess = PreProcess(conf)

subhalos = preprocess.get_subhalos()

subhaloIDs = subhalos['subhaloIDs']

for i, ID in enumerate(subhaloIDs):
    
    # Dynamically modify the config for each subhalo
    arguments = {}
    arguments['faceAndEdge'] = False
    arguments['numViews'] = 2
    arguments['randomViews'] = True
    
    preprocess.subhalo(ID)
    preprocess.prepare(arguments)
    # The modifications for config will be recorded and saved
    
    dataGeneration = DataGeneration(config=conf)
    dataGeneration.runSKIRT()
    
    postprocess = PostProcess(subhaloID=ID)
    postprocess.runPostprocess()