import pytest
import numpy as np
from quicklines.quicklines import gaussian
from quicklines.quicklines import Galaxy

# UNIT TESTS
# Test the Fitting Routine
# Mock Gaussian Profile -- (Mock Line Flux (1e-17), Observer Frame Central Wavelength (8000A), Sigma = 2A, Cont = 1e-20)
# Run the Gaussian Fit within the Line Class
# Compare the best-fit parameters versus the Mock parameters ()

def test_fits():
    wave = np.linspace(7950, 8950, 1000)
    lineflux = 1e-17
    cent_wave = 8000
    sigma = 2
    cont = 1e-20

    flux = gaussian(wave, lineflux, cent_wave, sigma, cont)
    error_flux = flux * 0.2

    my_galaxy = Galaxy(id=12346, zSpec=0.6, wave=wave, flux=flux, err=error_flux)

    mock_line = my_galaxy.run_line(5000,window=100,plot=False)

    mock_line.fit()

# Test getLineFlux()
# if we get the best-fit line flux returned in two cases: include_err = False, include_err = True

# Test getContinuumFluxDensity

# Test getVelocityDisp






# END-TO-END TEST
# Guideline
def test_all_zCOSMOS_spectra():
    pedro = Galaxy(803032)
    hb = pedro.run_line(4959.)

    print(hb.getLineFlux())
    print(hb.getContinuumFluxDensity())
    print(hb.getVelocityDisp())
    print(hb.getVelocityDisp(units="km/s",include_err=True))


# END-TO-END TEST WHERE USER PROVIDES THE SPECTRA
def test_all_zCOSMOS_spectra():
    pedro = Galaxy(222,wave=,flux=,err=,zSpec=)
    hb = pedro.run_line(4959.)

    print(hb.getLineFlux())
    print(hb.getContinuumFluxDensity())
    print(hb.getVelocityDisp())
    print(hb.getVelocityDisp(units="km/s",include_err=True))

