[![Documentation Status](https://readthedocs.org/projects/pysorb/badge/?version=latest)](https://pysorb.readthedocs.io/en/latest/?badge=latest)

# pysorb

PySORB is a Python module that provides implementations of various adsorption isotherm models for data fitting, analysis and simulation.

## Installation

You can install pysorb using pip:

```shell
    pip install pysorb
```

## Usage

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
- Relative Absolute Error (`'RABS'`)
- Hybrid Deviation (`'HYBRID'`)
- Mean Percentage Standard Deviation (`'MPSD'`)
- Sum of Relative Errors (`'SRE'`)
- Chi-square Deviation ('`CHI_2'`)
- Residual Sum of Squares Deviation (`'R_S'`)

# Contributing

Contributions are welcome! If you have suggestions, bug reports, or feature requests, please open an issue on the [GitHub repository](https://github.com/mv-per/pysorb).

# License

PySorb is licensed under the GNU GENERAL License. See the [LICENSE](https://github.com/mv-per/pysorb/blob/main/LICENSE) file for more information.
