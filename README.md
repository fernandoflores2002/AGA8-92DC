# AGA8-92DC Equation of State Implementation

## Overview

This repository contains a Python implementation of the AGA8-92DC equation of state (EOS) used for calculating the compressibility factor \( Z \) of natural gas mixtures. The code processes data from an Excel file to compute various parameters and applies the AGA8-92DC model to predict the behavior of natural gas.

## Version

**0.1** - Initial release with basic functionality for calculating the compressibility factor \( Z \) using the AGA8-92DC EOS.

## Files

- `aga92.py`: Python script implementing the AGA8-92DC equation of state.
- `coef_aga.xlsx`: Excel file containing coefficients and parameters required for the calculations.

## Installation

To use this code, you need to have Python installed along with the following libraries:

- `pandas`
- `numpy`
- `scipy`

You can install these dependencies using `pip`:

```bash
pip install pandas numpy scipy

## Usage

1. Place the coef_aga.xlsx file in the same directory AGA8-92DC.py.

2. Run the Python script:

```python 
python AGA8-92DC.py

```

3. Check the output for the calculated compressibility factor Z

## Example

The provided script includes a sample calculation with the following parameters:

```python 
P = 6657160
T = 225
R = 8.314  # Gas constant
x = np.array([0.90644, 0.04553, 0.00833, 0.001, 0.00156, 0.0003, 0.00045, 0.0004, 0.03134, 0.00466])
```

You can modify these values to suit your specific needs.

## Contact

For any questions or issues, please contact me at cesar.flores.b@uni.pe.