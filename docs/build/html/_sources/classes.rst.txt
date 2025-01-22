Classes Reference
===============

Configuration
-----------

``config.Configuration(surveys=None)``

``surveys`` (list)
    considered surveys. If None, configs will be read from current directory, or ``config.toml`` files will be created from template if there is no config files.

Instance Methods
~~~~~~~~~~~~~~

``conf = config.get_config()``
    update and return configs.

``config.add_survey(surveys)``
    ``surveys`` (list), surveys will be added. Call ``config.get_config()`` to get the updated configs.

``config.remove_survey(surveys)``
    ``surveys`` (list), surveys will be removed. Call ``config.get_config()`` to get the updated configs.

``config.save_config(conf)``
    ``conf`` (dict), save the current config to ``config.toml`` and ``config_(survey).toml``.

``config.check_config()``
    check if config has conflicts.


PreProcess
--------------

``preprocess = PreProcess(config)``
    ``config`` (dict), configuration.


Instance Methods
~~~~~~~~~~~~~~

``subhalos = preprocess.get_subhalos()``
    return subhalos including ``subhaloNum``, ``subhaloIDs`` and ``subhaloSFR``.

``preprocess.subhalo(subhaloID)``
    ``subhaloID`` (int), subhalo ID. Specify a subhalo.

``preprocess.prepare(data=None)``
    ``data`` (dict), For dynamical modification of parameters.
    Prepare stars, star-forming regions and dusts for current subhalo, and create SKIRT execution file.

``preprocess.inputs(data)``
    ``data`` (dict), Inputs for hydro-simulation except IllustrisTNG.
    data must include ``snapRedshift``, ``cosmology``, ``stellarMass``, ``subhaloID``, and ``boxLength`` used to select particles.

DataGeneration
------------

``dataGeneration = DataGeneration(config)``
    ``config`` (dict), configuration.


Instance Methods
~~~~~~~~~~~~~~

``dataGeneration.runSKIRT()``
    Perform SKIRT radiative transfer. 

PostProcess
----------

``postprocess = PostProcess(subhaloID, config)``
    ``subhaloID`` (int), subhalo ID.
    ``config`` (dict), configuration.


Instance Methods
~~~~~~~~~~~~~~

``postprocess.runPostprocess(showImage=False)``
    ``showImage`` (bool), show image. Only work in jupyter environment.
    Postprocess to perform mock observation.

.. autoclass:: galaxyGenius.preprocess.PreProcess
   :members:
   :undoc-members:
   :show-inheritance:

   