from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit



class Galaxy():
    """A simple attempt to model a dog."""
    def __init__(self, id,zSpec=None):
        """Initialize name and age attributes."""
        self.id = id
        self.zSpec = self.find_z_spec()
        self.line = self.Line
        self.wave,self.flux,self.err = self.retrieve_1dspec()

    def find_z_spec(self):
        data = fits.open('zCOSMOS_20K.fits')[1].data
        this_one = data['object_id']==self.id

        if True in this_one:
            return print(data['Redshift'][this_one])
        else:
            raise ValueError("ID is not in catalog")

    def retrieve_1dspec(self):
        spec_1d = fits.open(f'examples/{self.id}_1d.fits')
        wave = spec_1d['WAVE']
        flux = spec_1d['FLUX_REDUCED']
        err = spec_1d['ERR']
        return (wave,flux,err)

    class Line():
        def __init__(self,linewave):
            self.linewave = linewave
    
        def fit(self):

            

pedro = Galaxy(701230)
hb = pedro.Line(4861.)
pedro.zSpec
print(hb.linewave)