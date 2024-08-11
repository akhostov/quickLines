Usage
============

Import the module as 

.. code-block::  python
    
    import quicklines.quicklines as quicklines


Then initialize the galaxy class which will take in your spectra. For example

.. code-block:: python

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


If your source of interest is a zCOSMOS 20K survey source, then you can simply initialize the galaxy class as:

.. code-block:: python

    zCOS2O_id = 805732
    my_galaxy = quicklines.Galaxy(id=zCOS20_id)

as long as the zCOSMOS spectra is placed within the ```examples``` folder and has the format: ```805732_1d.fits``` The package is hardcoded such that if a zCOSMOS ID is provided, then it will load in the 1D spectra and extract the spectroscopic redshift via the zCOSMOS 20K catalog found within the ESO Science Archive.

After initializing the Galaxy class, you will be presented with a statement similar to this:

.. code-block:: python

    Spectra covers rest-frame wavelengths between 3298 and 5808 Angstroms
    
which will notify you that you can investigate any line for which the rest-frame wavelength falls within this range. Outside of this range will raise an error.

At this point, you can now enter a line of interest (rest-frame wavelengths) as such

.. code-block:: python

    Hbeta_line = my_galaxy.run_line(4861.)

where in this example we are interest in the Hydrogen Balmer Line (4 --> 2 transition) about 4861 Angstroms. quickLines will immediately convert this to observer-frame wavelength and automatically do the emission line profile fitting and print out the following:

.. code-block:: python

    Line Flux: 8.496e-18 +- 1.136e-18 erg/s/cm2
    Central Wavelength: 8116.8 +- 0.4 Angstrom
    Sigma: 3.90 +- 0.39 Angstrom
    Continuum Flux Density: 8.306e-19 +- 8.153e-20 erg/s/cm2/A
    S/N:7.48
    Reduced Chi-Square: 1.54

    Old Redshift: 0.6691
    Refined Redshift: 0.6698 +- 0.0001
    
which provides a quick, on-th-fly information about your source. It also will open a matplotlib GUI display showing you the observed spectra along with the modeled Gaussian emission line profile.

You can also store the resulting properties by using the following functions:

.. code-block:: python

    # Returns the Line Flux and its associated error
    lineflux, lineflux_err = Hbeta_line.getLineFlux(include_err=True)

    # Returns the Continuum Flux Density about the emission line wavelength
    cont_flam = Hbeta_line.getContinuumFluxDensity()

    # Returns the Velocity Dispersion in units of inputted wavelength (e.g., Angstrom)
    velDisp = Hbeta_line.getVelocityDisp()

    # Returns the Velocity Dispersion in units of km/s and includes the associated error
    velDisp,velDisp_err = Hbeta_line.getVelocityDisp(units="km/s", include_err=True)

