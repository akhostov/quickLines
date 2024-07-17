from astropy.io import fits
import astropy.constants as const
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Placeholder for the Gaussian Profile
def gaussian(x,lineFlux,cent_wave,sigma,cont):
    return None

def reduced_Chi2(x,model,obs,err):
    return 0.5*np.sum( ((model - obs)/err)**2. )


class Galaxy():
    
    # Initialize the Galaxy Class
    def __init__(self, id,zSpec=None):
        
        # Initialze the id which is the main user input
        # spectroscopic redshift is dervied from cross matching with the zCOSMOS 20K catalog
        # A Child Class Line will inherit all tha Galaxy Class attributes
        # The 1D spectra is retrieved from the examples folder
        self.id = id
        self.zSpec = self.find_z_spec()
        self.line = self.Line
        self.wave,self.flux,self.err = self.retrieve_1dspec()

        # Rest-Frame Wavelength Range Observed
        print(f"Spectra covers rest-frame wavelengths between {np.min(self.wave)} and {np.max(self.wave)} Angstroms")

    # Match with the zCOSMOS catalog to extract the spectroscopic redshift
    def find_z_spec(self):
        data = fits.open('zCOSMOS_20K.fits')[1].data
        this_one = data['object_id']==self.id

        if True in this_one:
            return print(data['Redshift'][this_one])
        else:
            raise ValueError("ID is not in catalog")

    # Get the 1D spectra and load it
    def retrieve_1dspec(self):
        spec_1d = fits.open(f'examples/{self.id}_1d.fits')
        wave = spec_1d['WAVE']
        flux = spec_1d['FLUX_REDUCED']
        err = spec_1d['ERR']
        return (wave,flux,err)


    # This class defines all the emission line profile measurements
    class Line():

        # Initialize the attributes
        def __init__(self,linewave):
            self.linewave = linewave

            # TODO: Need to include a check
            # if linewave < minimum covered wavelength (rest-frame) OR linewave > max wavelength (rest-frame):
            # then raise error and exit 

            self.best_model = self.fit()

        # Fit a Simple Gaussian Profile to the emission line
        def fit(self):

            # Fit the Model


            # Plot the Model against Observations


            # Print out the Reduced Chi^2


            return None # Placeholder

        def getLineFlux(self):
            return self.best_model[0]
        
        def getContinuumFluxDensity(self):
            return self.best_model[3]

        def getVelocityDisp(self,units="Angstrom"):

            if units == "Angstrom":
                return self.best_model[2]
            
            if units == "km/s":
                return self.best_model[2]/self.best_model[1]*const.c.to('km/s').value
            
            # In case user uses the wrong units
            if ("Angstrom" not in units) or ("km/s" not in units):
                raise ValueError("Available Units are Angstrom and km/s")


            

pedro = Galaxy(701230)
hb = pedro.Line(4861.)
pedro.zSpec
print(hb.linewave)