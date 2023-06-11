# ClassicIsotherms

ClassicIsotherms is a Python module that provides implementations of various classic adsorption isotherm equations. It allows you to calculate the loading and deviations for different isotherm models based on pressure, temperature, and parameters.

## Installation

You can install ClassicIsotherms using pip:

```shell
    pip install classic-isotherms
```

## Usage

```python
    from ClassicIsotherms import ClassicIsotherms

    # Create an instance of ClassicIsotherm
    isotherm = ClassicIsotherms("freundlich")

    # Calculate loading for a given pressure, temperature, and parameters
    pressure = 1.5
    temperature = 298.15
    parameters = [2.0, 1.5]
    loading = isotherm.get_loading(pressure, temperature, parameters)
    print("Loading:", loading)

    # Calculate loadings for multiple pressures
    pressures = [1.0, 2.0, 3.0]
    loadings = isotherm.get_loadings(pressures, temperature, parameters)
    print("Loadings:", loadings)

    # Calculate deviation using a specific deviation function
    experimental_loadings = [0.5, 1.0, 1.5]
    deviation_equation = "SSE"
    deviation = isotherm.get_deviation(pressures, experimental_loadings, temperature, parameters, deviation_equation)
    print("Deviation:", deviation)
```

# Supported Isotherm Models

The following classic adsorption isotherm models are currently supported:

- Freundlich (`'freundlich'`)
- Langmuir (`'langmuir'`)
- Redlich-Peterson (`'redlich-perterson'`)
- Sips (`'sips'`)
- Toth (`'toth'`)
- Unilan (`'unilan'`)
- Keller-Staudt-Toth (`'keller-staudt-toth'`)

# Supported Deviation Methods

The following deviation methods are currently supported:

- Sum of Squared Errors (`'SSE'`)
- Average Relative Error (`'ARE'`)
- Absolute Error (`'EABS'`)
- Hybrid Deviation (`'HYBRID'`)
- Mean Percentage Standard Deviation (`'MPSD'`)
- Sum of Relative Errors (`'SRE'`)
- Chi-square Deviation ('`CHI_2'`)
- Residual Sum of Squares Deviation (`'R_S'`)

# Contributing

Contributions are welcome! If you have suggestions, bug reports, or feature requests, please open an issue on the [GitHub repository](https://github.com/mv-per/ClassicIsotherms).

# License

ClassicIsotherms is licensed under the MIT License. See the [LICENSE](https://github.com/mv-per/ClassicIsotherms/blob/main/LICENSE) file for more information.
