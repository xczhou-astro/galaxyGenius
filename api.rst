API Reference
=============

This page provides detailed information about the galaxyGenius API.

Core Classes
-----------

Configuration
~~~~~~~~~~~~

.. autoclass:: galaxyGenius.config.Configuration
   :members:
   :undoc-members:
   :show-inheritance:

   The Configuration class handles all configuration settings for galaxyGenius, including:
   
   - Main configuration settings
   - Survey-specific configurations
   - Configuration validation and management

PreProcess
~~~~~~~~~

.. autoclass:: galaxyGenius.preprocess.PreProcess
   :members:
   :undoc-members:
   :show-inheritance:

   The PreProcess class handles the preparation of data for SKIRT simulation:
   
   - Reading and processing subhalo data
   - Preparing particle data
   - Creating SKIRT input files
   - Managing simulation parameters

DataGeneration
~~~~~~~~~~~~~

.. autoclass:: galaxyGenius.generation.DataGeneration
   :members:
   :undoc-members:
   :show-inheritance:

   The DataGeneration class manages the SKIRT radiative transfer simulation process:
   
   - Running SKIRT simulations
   - Managing input/output files
   - Handling data cube generation

PostProcess
~~~~~~~~~~

.. autoclass:: galaxyGenius.postprocess.PostProcess
   :members:
   :undoc-members:
   :show-inheritance:

   The PostProcess class handles the post-processing of simulation results:
   
   - Processing data cubes
   - Generating bandpass images
   - Creating SEDs
   - Visualizing results

Utility Functions
---------------

.. automodule:: galaxyGenius.utils
   :members:
   :undoc-members:
   :show-inheritance:

   Utility functions for various operations including:
   
   - Data conversion and manipulation
   - File handling
   - Mathematical operations
   - Visualization helpers 