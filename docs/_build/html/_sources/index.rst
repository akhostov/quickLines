.. QuickLines documentation master file, created by
   sphinx-quickstart on Thu Jul 18 16:35:57 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/quicklines_logo.png
   :alt: quickLines Logo
   :width: 800px
   :align: center


Welcome to quickLines
=====================
Developed as part of the 2024 Code/Astro Software Engineering for Astronomers Workshop, quickLines is designed and packaged to be a simple, quick way of extracting emission line properties on-the-fly given a 1D spectra of a galaxy. quickLines is meant to be a great companion in observing runs or exploring spectroscopic data sets with efficiency where we can essentially make on-the-fly calculations of emission line features before delving into more detailed calculations.

The program quickly measures:

- Emission Line Flux
- Continuum Flux Density
- Velocity Dispersion

The package takes in a given 1D spectra (wavelength, flux, and error spectrum), redshift, and a user-inputted line-of-interest (e.g., [OIII]5007A). quickLines then quickly fits a single gaussian line profile about the wavelength of interest (hence the "quick" in the name) and derives the above mentioned properties.

At the moment, this package can only fit single gaussian functions but can be expanded to include multi-component gaussians (e.g., narrow + broad emission line components, fitting doublets such as CIII]1907,1909A and [OII]3726,3728A) and also uses the scipy.optimize.curve_fit in the fitting process but can be expanded to incorporate other fitting approaches (e.g., Nested Sampling) in the case of complex parameter spaces (e.g., fitting multiple features).

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   install.rst
   usage.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
