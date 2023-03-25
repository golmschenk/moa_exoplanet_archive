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

merged_column_formats = {LightCurveColumnName.HJD: '.6f',
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
