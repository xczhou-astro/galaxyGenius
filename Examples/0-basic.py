
from galaxyGenius.config import Configuration
from galaxyGenius.generation import DataGeneration
from galaxyGenius.preprocess import PreProcess
from galaxyGenius.postprocess import PostProcess

config = Configuration()

config.add_survey('HSC')
config.add_survey('CSST')

conf = config.get_config()

preprocess = PreProcess(conf)

subhalos = preprocess.get_subhalos()

subhaloIDs = subhalos['subhaloIDs']

for i, ID in enumerate(subhaloIDs):
    
    preprocess.subhalo(ID)
    preprocess.prepare()
    
    dataGeneration = DataGeneration(config=conf)
    dataGeneration.runSKIRT()
    
    postprocess = PostProcess(subhaloID=ID, config=conf)
    postprocess.runPostprocess()

    