from typing import Dict

try:
    from enum import StrEnum
except ImportError:
    from backports.strenum import StrEnum


class LightCurveColumnName(StrEnum):
    HJD = 'HJD'
    COR_FLUX = 'cor_flux'
    FLUX = 'flux'
    FLUX_ERR = 'flux_err'
    OBS_ID = 'obsID'
    JD = 'JD'
    FWHM = 'fwhm'
    SKY = 'sky'
    AIRMASS = 'airmass'
    NSTAR = 'nstar'
    SCALE = 'scale'
    EXPTIME = 'exptime'
    SKYDIFF = 'skydiff'
    CHISQ = 'chisq'
    NPIX = 'npix'
    AIRMASS_1 = 'airmass1'
    ANG_1 = 'ang1'
    INCLUDED = 'included'


phot_all_column_names = [LightCurveColumnName.HJD, LightCurveColumnName.FLUX, LightCurveColumnName.FLUX_ERR,
                         LightCurveColumnName.OBS_ID, LightCurveColumnName.JD, LightCurveColumnName.FWHM,
                         LightCurveColumnName.SKY, LightCurveColumnName.AIRMASS, LightCurveColumnName.NSTAR,
                         LightCurveColumnName.SCALE, LightCurveColumnName.EXPTIME, LightCurveColumnName.SKYDIFF,
                         LightCurveColumnName.CHISQ, LightCurveColumnName.NPIX]

phot_column_names = [LightCurveColumnName.HJD, LightCurveColumnName.FLUX, LightCurveColumnName.FLUX_ERR,
                     LightCurveColumnName.OBS_ID, LightCurveColumnName.JD, LightCurveColumnName.FWHM,
                     LightCurveColumnName.SKY, LightCurveColumnName.AIRMASS, LightCurveColumnName.NSTAR,
                     LightCurveColumnName.SCALE, LightCurveColumnName.EXPTIME, LightCurveColumnName.SKYDIFF,
                     LightCurveColumnName.CHISQ, LightCurveColumnName.NPIX]

phot_cor_column_names = [LightCurveColumnName.HJD, LightCurveColumnName.COR_FLUX, LightCurveColumnName.FLUX_ERR,
                         LightCurveColumnName.OBS_ID, LightCurveColumnName.JD, LightCurveColumnName.FWHM,
                         LightCurveColumnName.SKY, LightCurveColumnName.AIRMASS, LightCurveColumnName.NSTAR,
                         LightCurveColumnName.SCALE, LightCurveColumnName.EXPTIME, LightCurveColumnName.SKYDIFF,
                         LightCurveColumnName.CHISQ, LightCurveColumnName.NPIX, LightCurveColumnName.AIRMASS_1,
                         LightCurveColumnName.ANG_1]

merged_column_names = [LightCurveColumnName.HJD, LightCurveColumnName.FLUX, LightCurveColumnName.COR_FLUX,
                       LightCurveColumnName.FLUX_ERR, LightCurveColumnName.OBS_ID, LightCurveColumnName.JD,
                       LightCurveColumnName.FWHM, LightCurveColumnName.SKY, LightCurveColumnName.AIRMASS,
                       LightCurveColumnName.NSTAR, LightCurveColumnName.SCALE, LightCurveColumnName.EXPTIME,
                       LightCurveColumnName.SKYDIFF, LightCurveColumnName.CHISQ, LightCurveColumnName.NPIX,
                       LightCurveColumnName.AIRMASS_1, LightCurveColumnName.ANG_1, LightCurveColumnName.INCLUDED]

light_curve_text_format_dictionary = {LightCurveColumnName.HJD: '.6f',
                                      LightCurveColumnName.FLUX: '.6f',
                                      LightCurveColumnName.COR_FLUX: '.6f',
                                      LightCurveColumnName.FLUX_ERR: '.6f',
                                      LightCurveColumnName.OBS_ID: '.0f',
                                      LightCurveColumnName.JD: '.6f',
                                      LightCurveColumnName.FWHM: '.5f',
                                      LightCurveColumnName.SKY: '.2f',
                                      LightCurveColumnName.AIRMASS: '.5f',
                                      LightCurveColumnName.NSTAR: '.0f',
                                      LightCurveColumnName.SCALE: '.6f',
                                      LightCurveColumnName.EXPTIME: '.0f',
                                      LightCurveColumnName.SKYDIFF: '.4f',
                                      LightCurveColumnName.CHISQ: '.6f',
                                      LightCurveColumnName.NPIX: '.0f',
                                      LightCurveColumnName.AIRMASS_1: '.4f',
                                      LightCurveColumnName.ANG_1: '.4f'}

light_curve_type_dictionary: Dict[LightCurveColumnName, str] = {
    LightCurveColumnName.HJD: 'float64',
    LightCurveColumnName.FLUX: 'float64',
    LightCurveColumnName.COR_FLUX: 'float64',
    LightCurveColumnName.FLUX_ERR: 'float64',
    LightCurveColumnName.OBS_ID: 'int32',
    LightCurveColumnName.JD: 'float64',
    LightCurveColumnName.FWHM: 'float32',
    LightCurveColumnName.SKY: 'float32',
    LightCurveColumnName.AIRMASS: 'float32',
    LightCurveColumnName.NSTAR: 'int32',
    LightCurveColumnName.SCALE: 'float32',
    LightCurveColumnName.EXPTIME: 'long',
    LightCurveColumnName.SKYDIFF: 'float32',
    LightCurveColumnName.CHISQ: 'float32',
    LightCurveColumnName.NPIX: 'int32',
    LightCurveColumnName.AIRMASS_1: 'float32',
    LightCurveColumnName.ANG_1: 'float32',
    LightCurveColumnName.INCLUDED: 'bool'
}

light_curve_description_dictionary: Dict[LightCurveColumnName, str] = {
    LightCurveColumnName.HJD: 'Heliocentric Julian Date of observational epoch.',
    LightCurveColumnName.FLUX: 'Number of counts measured for observational epoch: for a MOA magnitude of 20, this is 691.8 for magnitude of 20 for chip2 and 1,445 for all other chips.',
    LightCurveColumnName.COR_FLUX: 'Detrended flux.',
    LightCurveColumnName.FLUX_ERR: 'Uncertainty in measured flux for observational epoch.',
    LightCurveColumnName.OBS_ID: 'Run number, i.e., B#.',
    LightCurveColumnName.JD: 'Julian Date of observational epoch.',
    LightCurveColumnName.FWHM: 'Seeing estimated from the full-width half-maximum (FWHM) in arcseconds.',
    LightCurveColumnName.SKY: 'Sky level in counts on raw images.',
    LightCurveColumnName.AIRMASS: 'Airmass of observational epoch.',
    LightCurveColumnName.NSTAR: 'Number of resolved bright stars by "probe" (it is small if there are problems on the image).',
    LightCurveColumnName.SCALE: 'Numerical factor used to scale image relative to the reference image.',
    LightCurveColumnName.EXPTIME: 'Exposure time of the observational epoch in seconds.',
    LightCurveColumnName.SKYDIFF: 'Sky level in counts on the difference image (a residual of the Difference Imaging Analysis).',
    LightCurveColumnName.CHISQ: 'Chi-square of PSF fitting on the difference images. Two parameters are used in PSF fitting (χ2/dof is (/npix-2)).',
    LightCurveColumnName.NPIX: 'Number of pixels used in PSF fitting.',
    LightCurveColumnName.AIRMASS_1: 'Airmass calculated in the detrending code.',
    LightCurveColumnName.ANG_1: 'Angle relative to zenith calculated in the detrending code.',
    LightCurveColumnName.INCLUDED: 'Data included in fitting processes. Some data excluded due to bad quality flags.'
}
