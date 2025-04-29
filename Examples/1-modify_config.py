
from galaxyGenius.config import Configuration
from galaxyGenius.generation import DataGeneration
from galaxyGenius.preprocess import PreProcess
from galaxyGenius.postprocess import PostProcess

config = Configuration()

config.add_survey('HSC')
config.add_survey('CSST')

conf = config.get_config()

# Modify config
conf['requests'] = True
conf['apiKey'] = 'your_api_key'
conf['simulationMode'] = 'DustEmission'
conf['includeDust'] = True
conf['numPackets'] = 1e6

# Modift survey config
conf['includeBkg_CSST'] = False
conf['includePSF_CSST'] = False

preprocess = PreProcess(conf)

subhalos = preprocess.get_subhalos()

subhaloIDs = subhalos['subhaloIDs']

for i, ID in enumerate(subhaloIDs):
    
    # Dynamically modify the config for each subhalo
    data = {}
    data['faceAndEdge'] = False
    data['numViews'] = 2
    data['randomViews'] = True
    
    preprocess.subhalo(ID)
    preprocess.prepare(data)
    
    dataGeneration = DataGeneration(config=conf)
    dataGeneration.runSKIRT()
    
    postprocess = PostProcess(subhaloID=ID, config=conf)
    postprocess.runPostprocess()

# The modifications by config and data will be recorded in workingDir/config.json