# quickLines
This package is meant to extract spectral information by fitting galactic emission lines given certain 1D spectra data. As this far, the package only permits to fit gaussian functions. 

With this package you can;
  - Obtain the flux of an emission line
  - Obtain the flux of the continuum of an emission line
  - Obtain the velocity dispersion of an emission line

# Installation
The installation for this package is easy! just
```
pip install quicklines
```

# How to use
import the module as 
```
import quicklines.quicklines as quicklines
```
Then you can instance the Galaxy class like
```
superDuperCoolestGalaxy = quicklines.Galaxy(id=123456)
```
If the id doesn't exist in the zCOSMOS catalogue, you'll need to pass to the Galaxy class the array for the wavelength, flux and error in this order.

# Questions?
Do you have any questions or suggestions? Please send an email to
akhostov [at] gmail [dot] com or open an
`issue <https://github.com/akhostov/quickLines/issues>`

# License
The code is under the license of **MIT**
