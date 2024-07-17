from astropy.io import fits
import astropy.constants as const
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Placeholder for the Gaussian Profile
def gaussian(x,lineFlux,cent_wave,sigma,cont):
    return (1/sigma*np.sqrt(2*np.pi))*lineFlux * np.exp(-(x - cent_wave)**2 / (2 * sigma**2)) + cont

def reduced_Chi2(model,obs,err):
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
        self.wave,self.flux,self.err = self.retrieve_1dspec()

        # Rest-Frame Wavelength Range Observed
        print(f"Spectra covers rest-frame wavelengths between {np.min(self.wave)} and {np.max(self.wave)} Angstroms")

    # Match with the zCOSMOS catalog to extract the spectroscopic redshift
    def find_z_spec(self):
        try:
            data = fits.open('data/zCOSMOS_20K.fits')[1].data
        except:
            data = fits.open('../data/zCOSMOS_20K.fits')[1].data
        
        this_one = data['object_id']==self.id

        if True in this_one:
            return data['Redshift'][this_one][0]
        else:
            raise ValueError("ID is not in catalog")

    # Get the 1D spectra and load it
    def retrieve_1dspec(self):
        try:
            spec_1d = fits.open(f'examples/{self.id}_1d.fits')[1].data
        except:
            spec_1d = fits.open(f'../examples/{self.id}_1d.fits')[1].data
           
        wave = spec_1d['WAVE'][0]
        flux = spec_1d['FLUX_REDUCED'][0]
        #err = spec_1d['ERR'][0]
        err = np.abs(0.2*flux) # TODO: zCOSMOS has 0 for error but says there is 20% uncertainty in the fluxcalibration
        return (wave,flux,err)

    def run_line(self, linewave,**kwargs):
        return self.Line(self, linewave,**kwargs)
    
    # This class defines all the emission line profile measurements
    class Line():

        # Initialize the attributes
        def __init__(self,galaxy,linewave,window=40,plot=True):
            self.galaxy = galaxy
            self.linewave = linewave
            self.window = window
            self.plot = plot

            # This will trigger a check if wavelength inputted is in the 1D spectral coverage
            self.check_input_in_wave_coverage()


            # Fitting the Line
            self.best_param,self.best_perr = self.fit()

        def check_input_in_wave_coverage(self):
            obs_linewave = self.linewave*(1. + self.galaxy.zSpec)
            if (obs_linewave < np.min(self.galaxy.wave)) or \
                (obs_linewave > np.max(self.galaxy.wave)):
                raise ValueError("Inputted Wavelength is outside the coverage of the 1D Spectra!")

        # Fit a Simple Gaussian Profile to the emission line
        def fit(self):

            # Convert to Observed Wavelength
            obs_linewave = self.linewave*(1. + self.galaxy.zSpec)

            # Limit Fit
            keep = np.where((self.galaxy.wave > obs_linewave-self.window/2.) & (self.galaxy.wave < obs_linewave+self.window/2.))

            # Fit the Model
            params,pcov = curve_fit(gaussian,self.galaxy.wave[keep],
                                            self.galaxy.flux[keep],
                                            p0=[1e-17,obs_linewave,1.,0.],
                                            sigma=self.galaxy.err[keep])

            # Get Errors
            perr = np.sqrt(np.diag(pcov))

            # Get Bestfit Model
            bestfit_model = gaussian(self.galaxy.wave[keep],*params)

            # Print Parameters
            print("")
            print(f"Line Flux: {params[0]:.3e} +- {perr[0]:.3e} erg/s/cm2")
            print(f"Central Wavelength: {params[1]:.1f} +- {perr[1]:.1f} Angstrom")
            print(f"Sigma: {params[2]:.2f} +- {perr[2]:.2f} Angstrom")
            print(f"Continuum Flux Density: {params[3]:.3e} +- {perr[3]:.3e} erg/s/cm2/A")
            print(f'S/N:{params[0]/perr[0]:.2f}')
         
            # Get and Print Reduced Chi2
            dof = len(self.galaxy.flux[keep]) - len(params) 
            red_chi2 = reduced_Chi2(bestfit_model,self.galaxy.flux[keep],self.galaxy.err[keep])/dof
            print(f"Reduced Chi-Square: {red_chi2:0.2f}")

            # Get Refined Redshift
            new_specz = params[1]/self.linewave - 1.
            new_specz_err = perr[1]/self.linewave

            print("")
            print(f"Old Redshift: {self.galaxy.zSpec:0.4f}")
            print(f"Refined Redshift: {new_specz:0.4f} +- {new_specz_err:0.4f}")


            # Plot the Model against Observations
            if self.plot == True:
                plt.plot(self.galaxy.wave[keep],self.galaxy.flux[keep],label="Observed")
                plt.plot(self.galaxy.wave[keep],bestfit_model,ls='--',label="Model")
                plt.legend(loc="upper right",ncol=1,numpoints=1,fontsize=8,frameon=False)
                plt.show()

            return (params,perr) 

        def getLineFlux(self,include_err=False):
            if include_err == False:
                return self.best_param[0]
            if include_err == True:
                return (self.best_param[0],self.best_perr[0])


        def getContinuumFluxDensity(self,include_err=False):
            if include_err == False:
                return self.best_param[3]
            if include_err == True:
                return (self.best_param[3],self.best_perr[3])


        def getVelocityDisp(self,units="Angstrom",include_err=False):

            if units == "Angstrom":
                if include_err == True:
                    return (self.best_param[2],self.best_perr[2])
                if include_err == False:
                    return self.best_param[2]
            
            if units == "km/s":
                if include_err == True:
                    velDisp = self.best_param[2]/self.best_param[1]*const.c.to('km/s').value
                    velDisp_err = velDisp*np.sqrt( (self.best_perr[2]/self.best_param[2])**2. + (self.best_perr[1]/self.best_param[1])**2. )
                    return (velDisp,velDisp_err)
                
                if include_err == False:
                    return self.best_param[2]/self.best_param[1]*const.c.to('km/s').value
            
            # In case user uses the wrong units
            if ("Angstrom" not in units) or ("km/s" not in units):
                raise ValueError("Available Units are Angstrom and km/s")
    
def main():
    pedro = Galaxy(701230)
    hb = pedro.run_line(4959.)

    print(hb.getLineFlux())
    print(hb.getContinuumFluxDensity())
    print(hb.getVelocityDisp())
    print(hb.getVelocityDisp(units="km/s",include_err=True))

            
if __name__ == "__main__":
    main()