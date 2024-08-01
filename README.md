![alt text](quicklines/quicklines_logo.png)

Developed as part of the 2024 Code/Astro Software Engineering for Astronomers Workshop, quickLines is designed and packaged to be a simple, quick way of extracting emission line properties on-the-fly given a 1D spectra of a galaxy. quickLines is meant to be a great companion in observing runs or exploring spectroscopic data sets with efficiency where we can essentially make on-the-fly calculations of emission line features before devling into more detailed calculations.

The program quickly measures:
  - Emission Line Flux
  - Continuum Flux Density
  - Velocity Dispersion

The package takes in a given 1D spectra (wavelength, flux, and error spectrum), redshift, and a user-inputted line-of-interest (e.g., [OIII]5007A). quickLines then quickly fits a single gaussian line profile about the wavelength of interest (hence the "quick" in the name) and derives the above mentioned properties.

At the moment, this package can only fit single gaussian functions but can be expanded to include multi-component gaussians (e.g., narrow + broad emission line components, fitting doublets such as CIII]1907,1909A and [OII]3726,3728A) and also uses the ```scipy.optimize.curve_fit``` in the fitting process but can be expanded to incorporate other fitting approaches (e.g., Nested Sampling) in the case of complex parameter spaces (e.g., fitting multiple features).

# Installation
The installation for this package is easy! just
```
pip install quicklines
```
or you can also run
```
pip install git+https://github.com/akhostov/quickLines.git
```
or clone the repository.

# Usage
Import the module as 
```
import quicklines.quicklines as quicklines
```

Then initialize the galaxy class which will take in your spectra. For example
```
# Load in the Astropy Module that Reads a Fits File
from astropy.io import fits

# Load in my 1D spectra for the galaxy of interest
spec1D = fits.open("my_spectra.fits")[1].data

# Assign ID and Redshift
id = 123
zSpec = 0.84

# Extract the Wavalength, Flux Spectrum, and Error Spectrum (all observer-frame wavelengths and fluxes)
wave = spec1D["WAVELENGTH"]
flux = spec1D["FLUX"]
err = spec1D["ERROR"]

# Initialize the Galaxy Class
my_galaxy = quicklines.Galaxy(id=id, zSpec = zSpec, wave=wave, flux=flux, err=err)
```

If your source of interest is a zCOSMOS 20K survey source, then you can simply initialize the galaxy class as:
```
zCOS2O_id = 805732
my_galaxy = quicklines.Galaxy(id=zCOS20_id)
```
as long as the zCOSMOS spectra is placed within the ```examples``` folder and has the format: ```805732_1d.fits``` The package is hardcoded such that if a zCOSMOS ID is provided, then it will load in the 1D spectra and extract the spectroscopic redshift via the zCOSMOS 20K catalog found within the ESO Science Archive.

After initializing the Galaxy class, you will be presented with a statement similar to this:
```Spectra covers rest-frame wavelengths between 3298 and 5808 Angstroms```
which will notify you that you can investigate any line for which the rest-frame wavelength falls within this range. Outside of this range will raise an error.

At this point, you can now enter a line of interest (rest-frame wavelengths) as such
```
Hbeta_line = my_galaxy.run_line(4861.)
```
where in this example we are interest in the Hydrogen Balmer Line (4 --> 2 transition) about 4861 Angstroms. quickLines will immediately convert this to observer-frame wavelength and automatically do the emission line profile fitting and print out the following:
```
Line Flux: 8.496e-18 +- 1.136e-18 erg/s/cm2
Central Wavelength: 8116.8 +- 0.4 Angstrom
Sigma: 3.90 +- 0.39 Angstrom
Continuum Flux Density: 8.306e-19 +- 8.153e-20 erg/s/cm2/A
S/N:7.48
Reduced Chi-Square: 1.54

Old Redshift: 0.6691
Refined Redshift: 0.6698 +- 0.0001
```
which provides a quick, on-th-fly information about your source. It also will open a matplotlib GUI display showing you the observed spectra along with the modeled Gaussian emission line profile.

You can also store the resulting properties by using the following functions:
```
 # Returns the Line Flux and its associated error
lineflux, lineflux_err = Hbeta_line.getLineFlux(include_err=True)

# Returns the Continuum Flux Density about the emission line wavelength
cont_flam = Hbeta_line.getContinuumFluxDensity()

# Returns the Velocity Dispersion in units of inputted wavelength (e.g., Angstrom)
velDisp = Hbeta_line.getVelocityDisp()

# Returns the Velocity Dispersion in units of km/s and includes the associated error
velDisp,velDisp_err = Hbeta_line.getVelocityDisp(units="km/s", include_err=True)
```

# Questions?
Do you have any questions or suggestions? Please send an email to
akhostov [at] gmail [dot] com or open an
`issue <https://github.com/akhostov/quickLines/issues>`

# License
The code is under the license of **MIT**
