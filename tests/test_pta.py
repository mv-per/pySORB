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

@pytest.mark.parametrize(
    'adsorption_potential, parameters, expected_calculated_loading',
    [
        (DRA_POTENTIAL, [7880.19, 0.29, 2.], 5.7589324),
        (LEE_POTENTIAL, [125.63, 12.26, 765.70], 6.113060),
        (STEELE_POTENTIAL, [109.32, 13.34, 611.88], 7.50000),
    ]
)
def test_pure_pta_calculate_loading_pr77(setup_fluid:Fluid, adsorption_potential:str, parameters:List[float], expected_calculated_loading:float)->None:

    adsorbent = Adsorbent("Z01x", 3.35, 0.382)
    setup_fluid.lennard_jonnes_diameter = 3.941

    pure_pta = PotentialTheoryModels(adsorption_potential, 'pr77', 'excess', 555, setup_fluid, adsorbent)

    loading = pure_pta.get_loading(1e6, 305, parameters)

    assert loading == pytest.approx(expected_calculated_loading)

@pytest.mark.parametrize(
    'adsorption_potential, parameters, expected_calculated_loading',
    [
        (DRA_POTENTIAL, [7880.19, 0.29, 2.], 5.0415179),
        (LEE_POTENTIAL, [125.63, 12.26, 765.70], 5.3310869),
        (STEELE_POTENTIAL, [109.32, 13.34, 611.88], 6.7061698),
    ]
)
def test_pure_pta_calculate_loading_srk(setup_fluid:Fluid, adsorption_potential:str, parameters:List[float], expected_calculated_loading:float)->None:


    adsorbent = Adsorbent("Z01x", 3.35, 0.382)

    setup_fluid.lennard_jonnes_diameter = 3.941
    pure_pta = PotentialTheoryModels(adsorption_potential, 'srk', 'excess', 555, setup_fluid, adsorbent)

    loading = pure_pta.get_loading(1e6, 305, parameters)

    assert loading == pytest.approx(expected_calculated_loading)

@pytest.mark.parametrize(
    'adsorption_potential, parameters',
    [
        (LEE_POTENTIAL, [125.63, 12.26, 765.70]),
        (STEELE_POTENTIAL, [109.32, 13.34, 611.88]),
    ]
)
def test_pure_pta_calculate_loading_error_without_adsorbent(setup_fluid:Fluid, adsorption_potential:str, parameters:List[float])->None:

    pure_pta = PotentialTheoryModels(adsorption_potential, 'pr77', 'excess', 555, setup_fluid)

    expected_error_message = """Adsorbent properties are needed for LJ-based potentials and is not defined."""
    with pytest.raises(ValueError) as e:
        pure_pta.get_loading(1e6, 305, parameters)

    assert expected_error_message in str(e)


def test_pure_pta_calculate_loadings(setup_fluid:Fluid)->None:
    import numpy as np


    pure_pta = PotentialTheoryModels(DRA_POTENTIAL, 'pr77', 'excess', 555, setup_fluid)

    DRA_PARAMETERS = [7880.19, 0.29, 2.]

    pressures = np.arange(1e6, 11e6, 1e6)

    loadings = pure_pta.get_loadings(pressures, 305, DRA_PARAMETERS)

    assert len(loadings) == len(pressures)

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
    
    for calc, exp in zip(loadings, expected_loadings):
        assert pytest.approx(calc) == exp

@pytest.mark.parametrize(
    'type_of_deviation, expected',
    [
        ('RABS', 0.0422067131765064),
        ('EABS', 0.023953988843656404),
        ('ARE', -0.011296498267532325),
    ]
)
def test_pure_pta_calculate_deviation_of_range(setup_fluid:Fluid, type_of_deviation:str, expected:float)->None:
    import numpy as np

    pure_pta = PotentialTheoryModels(DRA_POTENTIAL, 'pr77', 'excess', 555, setup_fluid)

    DRA_PARAMETERS = [7880.19, 0.29, 2.]

    pressures = np.arange(1e6, 11e6, 1e6)

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

    experimental_loadings = [round(loading, 2) for loading in experimental_loadings]

    deviation = pure_pta.get_deviation(pressures, experimental_loadings, 305, DRA_PARAMETERS, type_of_deviation)

    assert pytest.approx(deviation) == expected