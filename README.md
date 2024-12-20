# SpectralCubePlus

SpectralCubePlus is an extension of the popular [SpectralCube](https://spectral-cube.readthedocs.io/) library, designed to provide additional functionality for handling spectral data cubes. This package includes methods for automated RMS calculation, error propagation for moment maps, and advanced masking operations.

## Features

- **Automated RMS Calculation:** Easily compute the RMS noise of your spectral data cube using `get_rms` or `get_rms_auto`.
- **Moment Maps with Error Propagation:** Generate moment maps (e.g., moment-0) and calculate associated errors using `moment0err`.
- **Enhanced Masking Tools:** Create and apply custom masks to your data cubes with `get_expmask`.
- **Seamless Integration with SpectralCube:** Inherits all functionality from the SpectralCube library while adding new features.

## Installation
## STILL TO DO.... 

To install SpectralCubePlus, clone the repository and install the package using pip:

```bash
$ git clone https://github.com/yourusername/spectral_cube_plus.git
$ cd spectral_cube_plus
$ pip install .
```

## Usage

### Reading a Data Cube

```python
from spectral_cube_plus import SpectralCubePlus

# Read a spectral data cube
cube = SpectralCubePlus.read("path/to/your/file.fits")
```

### Calculate RMS Noise

```python
# Automated RMS calculation
cube.get_rms_auto()
print(cube.rms)
```

### Generate Moment Maps and Error Propagation

```python
# Calculate the moment-0 map
moment0 = cube.moment0()

# Calculate the error for the moment-0 map
moment0_error = cube.moment0err()
```

### Advanced Masking

```python
# Generate an expanded mask based on RMS
expmask = cube.get_expmask(threshold=3)

# Apply the mask to the data cube
masked_cube = cube.with_mask(expmask)
```

## Requirements

- Python 3.7+
- [SpectralCube](https://spectral-cube.readthedocs.io/)
- [Astropy](https://www.astropy.org/)
- [Numpy](https://numpy.org/)

## Contributing

Contributions are welcome! Please follow these steps to contribute:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Write your code and add tests.
4. Submit a pull request.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [SpectralCube](https://spectral-cube.readthedocs.io/) for providing the foundation for this package.
- [Astropy](https://www.astropy.org/) for its invaluable tools for astronomy.

## Contact

For questions or suggestions, please open an issue on the GitHub repository or contact Ashley Barnes