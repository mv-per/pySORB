# Classic Isotherms

```python

def test_get_mixture_loading_two_components(setup_fluids:List[Fluid], potential:str, CO2_params:List[float], CH4_params:List[float], expected:List[float])->None:

    bulk_composition = [.25, .75]
    mixture_params = [CO2_params, CH4_params]
    fluids = setup_fluids[:2]
    adsorbent = Adsorbent("Z01x", 3.35, 0.382)

    model = PotentialTheoryModels(potential, 'pr77', 'excess', 555, fluids,adsorbent)


    calculated = model.get_mixture_loading( 1e6, 310.2, bulk_composition,mixture_params)

    calculated = [round(val, 5) for val in calculated]

    assert len(calculated) == len(bulk_composition)
    assert calculated == expected


```
