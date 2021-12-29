try:
    from enum import StrEnum
except ImportError:
    from backports.strenum import StrEnum


class ColumnName(StrEnum):
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
