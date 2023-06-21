.. PySorb documentation master file, created by
   sphinx-quickstart on Wed Jun 21 20:04:41 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PySorb's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. currentmodule:: pysorb

===========
pysorb Package
===========

C/C++ methods to run the PTA model.

EmpiricalIsotherms
------------------
.. autoclass:: EmpiricalIsotherms
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: EmpiricalIsotherms.get_pure_loading
      :noindex:
   
   .. automethod:: EmpiricalIsotherms.get_pure_loadings
      :noindex:
   
   .. automethod:: EmpiricalIsotherms.get_deviation
      :noindex:
   
   .. automethod:: EmpiricalIsotherms.get_mixture_loading
      :noindex:

VacancySolutionMethod
---------------------
.. autoclass:: VacancySolutionMethod
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: VacancySolutionMethod.get_pure_loading
      :noindex:
   
   .. automethod:: VacancySolutionMethod.get_pure_loadings
      :noindex:
   
   .. automethod:: VacancySolutionMethod.get_deviation
      :noindex:
   
   .. automethod:: VacancySolutionMethod.get_mixture_loading
      :noindex:

PotentialTheoryModels
---------------------
.. autoclass:: PotentialTheoryModels
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: PotentialTheoryModels.get_pure_loading
      :noindex:
   
   .. automethod:: PotentialTheoryModels.get_pure_loadings
      :noindex:
   
   .. automethod:: PotentialTheoryModels.get_deviation
      :noindex:
   
   .. automethod:: PotentialTheoryModels.get_mixture_loading
      :noindex:

   ============
   Example Usage
   ============

   This section provides some example usages of the `PotentialTheoryModels` class.

   Creating a Pure Potential Theory Adsorption Model
   -------------------------------------------------

   The following example demonstrates how to create a pure potential theory adsorption model using different adsorption potentials.

   .. code-block:: python

      from pysorb import PotentialTheoryModels, Fluid, Adsorbent
      import pytest
      from typing import List
      import numpy

      DRA_POTENTIAL = 'DRA'
      LEE_POTENTIAL = 'LEE'
      STEELE_POTENTIAL = 'STEELE'

      @pytest.fixture
      def setup_fluid() -> Fluid:
         return Fluid('CO2', 73.773e5, 304.13, 0.22394)

      @pytest.mark.parametrize(
         'adsorption_potential', (DRA_POTENTIAL, LEE_POTENTIAL, STEELE_POTENTIAL)
      )
      def test_create_pure_pta(adsorption_potential:str, setup_fluid:Fluid)->None:
         solver = PotentialTheoryModels(adsorption_potential, 'pr77', 'excess', 555, setup_fluid)
         assert solver is not None

      # Rest of the test cases...

   Calculating Loading in a Pure Potential Theory Adsorption Model
   --------------------------------------------------------------

   The following example demonstrates how to calculate the loading in a pure potential theory adsorption model.

   .. code-block:: python

      # Setup fluid and adsorbent objects
      adsorbent = Adsorbent("Z01x", 3.35, 0.382)
      setup_fluid.lennard_jonnes_diameter = 3.941

      # Create the PotentialTheoryModels instance
      pure_pta = PotentialTheoryModels(adsorption_potential, 'pr77', 'excess', 555, setup_fluid, adsorbent)

      # Calculate the loading
      loading = pure_pta.get_loading(1e6, 305, parameters)

      # Assert the calculated loading
      assert loading == pytest.approx(expected_calculated_loading)

   Calculating Loadings in a Pure Potential Theory Adsorption Model
   --------------------------------------------------------------

   The following example demonstrates how to calculate multiple loadings in a pure potential theory adsorption model.

   .. code-block:: python

      import numpy as np

      # Create the PotentialTheoryModels instance
      pure_pta = PotentialTheoryModels(DRA_POTENTIAL, 'pr77', 'excess', 555, setup_fluid)

      # Define the parameters
      DRA_PARAMETERS = [7880.19, 0.29, 2.]

      # Define the pressures
      pressures = np.arange(1e6, 11e6, 1e6)

      # Calculate the loadings
      loadings = pure_pta.get_loadings(pressures, 305, DRA_PARAMETERS)

      # Assert the length of loadings and pressures
      assert len(loadings) == len(pressures)

      # Define the expected loadings
      expected_loadings = [5.758932400768459,
               6.895706537068475,
               7.304286537152829,
               7.42048187387874,
               7.368745932028131,
               7.173182017172223,
               6.720995771700607,
               4.603105939555204,
               4.11310650481496,
               3.8378197855658422]

      # Assert each calculated loading with the expected value
      for calc, exp in zip(loadings, expected_loadings):
         assert pytest.approx(calc) == exp

   Calculating Deviation in a Pure Potential Theory Adsorption Model
   ----------------------------------------------------------------

   The following example demonstrates how to calculate the deviation in a pure potential theory adsorption model.

   .. code-block:: python

      import numpy as np

      # Create the PotentialTheoryModels instance
      pure_pta = PotentialTheoryModels(DRA_POTENTIAL, 'pr77', 'excess', 555, setup_fluid)

      # Define the parameters
      DRA_PARAMETERS = [7880.19, 0.29, 2.]

      # Define the pressures
      pressures = np.arange(1e6, 11e6, 1e6)

      # Define the experimental loadings
      experimental_loadings = [5.758932400768459,
               6.895706537068475,
               7.304286537152829,
               7.42048187387874,
               7.368745932028131,
               7.173182017172223,
               6.720995771700607,
               4.603105939555204,
               4.11310650481496,
               3.8378197855658422]

      # Round the experimental loadings to 2 decimal places
      experimental_loadings = [round(loading, 2) for loading in experimental_loadings]

      # Calculate the deviation
      deviation = pure_pta.get_deviation(pressures, experimental_loadings, 305, DRA_PARAMETERS, type_of_deviation)

      # Assert the calculated deviation
      assert pytest.approx(deviation) == expected
