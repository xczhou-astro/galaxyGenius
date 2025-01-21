Usage
=====

Basic Usage
----------

Here's a basic example of how to use GalaxyGenius:

.. code-block:: python

    from galaxyGenius import PreProcess
    
    # Initialize configuration
    config = {
        'dataDir': '/path/to/data',
        'workingDir': '/path/to/working/dir',
        'snapNum': 99,
        'simulation': 'TNG100-1',
        'minStellarMass': 1e10,
        'maxStellarMass': 1e11
    }
    
    # Create preprocessor instance
    preprocessor = PreProcess(config)
    
    # Get subhalos
    subhalos = preprocessor.get_subhalos()
    
    # Process specific subhalo
    preprocessor.subhalo(subhaloID=123)
    
    # Prepare data for analysis
    preprocessor.prepare()

Advanced Usage
-------------

For more advanced usage and configuration options, please refer to the API documentation. 