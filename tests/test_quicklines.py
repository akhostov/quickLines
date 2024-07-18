import pytest
import numpy as np
import astropy.constants as const
from astropy.io import fits

import os
import sys
sys.path.insert(0, os.path.abspath('../'))

from quicklines.quicklines import gaussian
from quicklines.quicklines import Galaxy

# UNIT TESTS
# Test the Fitting Routine
# Mock Gaussian Profile -- (Mock Line Flux (1e-17), Observer Frame Central Wavelength (8000A), Sigma = 2A, Cont = 1e-20)
# Run the Gaussian Fit within the Line Class
# Compare the best-fit parameters versus the Mock parameters ()

def test_fit():
    wave = np.linspace(7950, 8950, 1000)
    lineflux = 1e-17
    cent_wave = 8000
    sigma = 2
    cont = 1e-20

    flux = gaussian(wave, lineflux, cent_wave, sigma, cont)
    err_flux = flux * 0.2

    my_galaxy = Galaxy(id=12346, zSpec=0.6, wave=wave, flux=flux, err=err_flux)

    mock_line = my_galaxy.run_line(5000, window=100, plot=False)

    params, perr = mock_line.fit()
    assert lineflux == pytest.approx(params[0])
    assert cent_wave == pytest.approx(params[1])
    assert sigma == pytest.approx(params[2])
    assert cont == pytest.approx(params[3])


# Test getLineFlux()
# if we get the best-fit line flux returned in two cases: include_err = False, include_err = True
def test_getLineFlux():
    wave = np.linspace(7950, 8950, 1000)
    lineflux = 1e-17
    cent_wave = 8000
    sigma = 2
    cont = 1e-20

    flux = gaussian(wave, lineflux, cent_wave, sigma, cont)
    err_flux = flux * 0.2

    my_galaxy = Galaxy(id=12346, zSpec=0.6, wave=wave, flux=flux, err=err_flux)

    mock_line = my_galaxy.run_line(5000, window=100, plot=False)

    assert lineflux == pytest.approx(mock_line.getLineFlux())

    assert lineflux == pytest.approx(mock_line.getLineFlux(include_err=True)[0])



# Test getContinuumFluxDensity
def test_getContinuumFluxDensity():
    wave = np.linspace(7950, 8950, 1000)
    lineflux = 1e-17
    cent_wave = 8000
    sigma = 2
    cont = 1e-20

    flux = gaussian(wave, lineflux, cent_wave, sigma, cont)
    err_flux = flux * 0.2

    my_galaxy = Galaxy(id=12346, zSpec=0.6, wave=wave, flux=flux, err=err_flux)

    mock_line = my_galaxy.run_line(5000, window=100, plot=False)

    assert cont == pytest.approx(mock_line.getContinuumFluxDensity())

    assert cont == pytest.approx(mock_line.getContinuumFluxDensity(include_err=True)[0])

# Test getVelocityDisp
def test_getVelocityDisp():
    wave = np.linspace(7950, 8950, 1000)
    lineflux = 1e-17
    cent_wave = 8000
    sigma = 2
    cont = 1e-20

    flux = gaussian(wave, lineflux, cent_wave, sigma, cont)
    err_flux = flux * 0.2

    my_galaxy = Galaxy(id=12346, zSpec=0.6, wave=wave, flux=flux, err=err_flux)

    mock_line = my_galaxy.run_line(5000, window=100, plot=False)

    assert sigma == pytest.approx(mock_line.getVelocityDisp())
    assert sigma == pytest.approx(mock_line.getVelocityDisp(include_err=True)[0])

    sigma_kms = sigma / cent_wave * const.c.to('km/s').value

    assert sigma_kms == pytest.approx(mock_line.getVelocityDisp(units='km/s'))
    assert sigma_kms == pytest.approx(mock_line.getVelocityDisp(units='km/s', include_err=True)[0])


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
def test_all_user_spectra():
    data = fits.open("../examples/701230_1d.fits")[1].data

    wave = data["WAVE"][0]
    flux = data["FLUX_REDUCED"][0]
    err = flux*0.2

    pedro = Galaxy("123",wave=wave,flux=flux,err=err,zSpec=0.6691)

    hb = pedro.run_line(4959.)

    print(hb.getLineFlux())
    print(hb.getContinuumFluxDensity())
    print(hb.getVelocityDisp())
    print(hb.getVelocityDisp(units="km/s",include_err=True))

if __name__ == '__main__':
    test_fit()
    test_getLineFlux()
    test_getContinuumFluxDensity()
    test_getVelocityDisp()

    test_all_zCOSMOS_spectra()
    test_all_user_spectra()
