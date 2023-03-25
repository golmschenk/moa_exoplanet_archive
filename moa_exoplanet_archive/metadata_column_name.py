from dataclasses import dataclass
from typing import Optional, List

from astropy.coordinates import Angle
try:
    from enum import StrEnum
except ImportError:
    from backports.strenum import StrEnum


@dataclass
class AlertMetadata:
    separation_to_alert_position__pixels: float
    alert_id: str
    alert_x__pixels: float
    alert_y__pixels: float


@dataclass
class CandidateMetadata:
    number_of_data_points: int
    number_of_frames_object_is_detected: int
    max_significance: float
    sum_of_continuous_significance: float
    chi_squared_outside_search_window: float
    dophot_object_reference_image_separation__pixels: float
    dophot_id: int
    dophot_type: int
    dophot_magnitude: float
    dophot_magnitude_error: float
    pspl_t0: float
    pspl_tE: float
    pspl_u0: float
    pspl_source_flux: float
    pspl_blending_flux: float
    pspl_t0_parabolic_error: float
    pspl_tE_parabolic_error: float
    pspl_tE_error_lower_limit: float
    pspl_tE_error_upper_limit: float
    pspl_u0_parabolic_error: float
    pspl_u0_error_lower_limit: float
    pspl_u0_error_upper_limit: float
    pspl_source_flux_error: float
    pspl_blending_flux_error: float
    pspl_chi_squared: float
    fspl_t0: float
    fspl_tE: float
    fspl_u0: float
    fspl_rho: float
    fspl_source_flux: float
    fspl_blending_flux: float
    fspl_t0_parabolic_error: float
    fspl_tE_parabolic_error: float
    fspl_tE_error_lower_limit: float
    fspl_tE_error_upper_limit: float
    fspl_u0_parabolic_error: float
    fspl_u0_error_lower_limit: float
    fspl_u0_error_upper_limit: float
    fspl_rho_parabolic_error: float
    fspl_rho_error_lower_limit: float
    fspl_rho_error_upper_limit: float
    fspl_source_flux_error: float
    fspl_blending_flux_error: float
    fspl_chi_squared: float


@dataclass
class TargetMetadata:
    field: int
    chip: int
    subframe: int
    id: int
    tag: str
    x__pixels: float
    y__pixels: float
    ra_j2000: Angle
    dec_j2000: Angle
    candidate_metadata: Optional[CandidateMetadata] = None
    alert_metadata_list: Optional[List[AlertMetadata]] = None


class MetadataColumnName(StrEnum):
    FIELD = 'field'
    CHIP = 'chip'
    SUBFRAME = 'subframe'
    ID = 'id'
    TAG = 'tag'
    X = 'x'
    Y = 'y'
    RA_J2000 = 'ra_j2000'
    DEC_J2000 = 'dec_j2000'
    NUMBER_OF_DATA_POINTS = 'number_of_data_points'
    NUMBER_OF_FRAMES_OBJECT_IS_DETECTED = 'number_of_frames_object_is_detected'
    MAX_SIGNIFICANCE = 'max_significance'
    SUM_OF_CONTINUOUS_SIGNIFICANCE = 'sum_of_continuous_significance'
    CHI_SQUARED_OUTSIDE_SEARCH_WINDOW = 'chi_squared_outside_search_window'
    DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION = 'dophot_object_reference_image_separation'
    DOPHOT_ID = 'dophot_id'
    DOPHOT_TYPE = 'dophot_type'
    DOPHOT_MAGNITUDE = 'dophot_magnitude'
    DOPHOT_MAGNITUDE_ERROR = 'dophot_magnitude_error'
    PSPL_T0 = 'pspl_t0'
    PSPL_TE = 'pspl_tE'
    PSPL_U0 = 'pspl_u0'
    PSPL_SOURCE_FLUX = 'pspl_source_flux'
    PSPL_BLENDING_FLUX = 'pspl_blending_flux'
    PSPL_T0_PARABOLIC_ERROR = 'pspl_t0_parabolic_error'
    PSPL_TE_PARABOLIC_ERROR = 'pspl_tE_parabolic_error'
    PSPL_TE_ERROR_LOWER_LIMIT = 'pspl_tE_error_lower_limit'
    PSPL_TE_ERROR_UPPER_LIMIT = 'pspl_tE_error_upper_limit'
    PSPL_U0_PARABOLIC_ERROR = 'pspl_u0_parabolic_error'
    PSPL_U0_ERROR_LOWER_LIMIT = 'pspl_u0_error_lower_limit'
    PSPL_U0_ERROR_UPPER_LIMIT = 'pspl_u0_error_upper_limit'
    PSPL_SOURCE_FLUX_ERROR = 'pspl_source_flux_error'
    PSPL_BLENDING_FLUX_ERROR = 'pspl_blending_flux_error'
    PSPL_CHI_SQUARED = 'pspl_chi_squared'
    FSPL_T0 = 'fspl_t0'
    FSPL_TE = 'fspl_tE'
    FSPL_U0 = 'fspl_u0'
    FSPL_RHO = 'fspl_rho'
    FSPL_SOURCE_FLUX = 'fspl_source_flux'
    FSPL_BLENDING_FLUX = 'fspl_blending_flux'
    FSPL_T0_PARABOLIC_ERROR = 'fspl_t0_parabolic_error'
    FSPL_TE_PARABOLIC_ERROR = 'fspl_tE_parabolic_error'
    FSPL_TE_ERROR_LOWER_LIMIT = 'fspl_tE_error_lower_limit'
    FSPL_TE_ERROR_UPPER_LIMIT = 'fspl_tE_error_upper_limit'
    FSPL_U0_PARABOLIC_ERROR = 'fspl_u0_parabolic_error'
    FSPL_U0_ERROR_LOWER_LIMIT = 'fspl_u0_error_lower_limit'
    FSPL_U0_ERROR_UPPER_LIMIT = 'fspl_u0_error_upper_limit'
    FSPL_RHO_PARABOLIC_ERROR = 'fspl_rho_parabolic_error'
    FSPL_RHO_ERROR_LOWER_LIMIT = 'fspl_rho_error_lower_limit'
    FSPL_RHO_ERROR_UPPER_LIMIT = 'fspl_rho_error_upper_limit'
    FSPL_SOURCE_FLUX_ERROR = 'fspl_source_flux_error'
    FSPL_BLENDING_FLUX_ERROR = 'fspl_blending_flux_error'
    FSPL_CHI_SQUARED = 'fspl_chi_squared'
    SEPARATION_TO_ALERT_POSITION0 = 'separation_to_alert_position0'
    ALERT_ID0 = 'alert_id0'
    ALERT_X0 = 'alert_x0'
    ALERT_Y0 = 'alert_y0'
    SEPARATION_TO_ALERT_POSITION1 = 'separation_to_alert_position1'
    ALERT_ID1 = 'alert_id1'
    ALERT_X1 = 'alert_x1'
    ALERT_Y1 = 'alert_y1'


class TakahiroSumiCandlistRaDecFileColumnName(StrEnum):
    FIELD = 'field'
    CHIP = 'chip'
    NSUB = 'nsub'
    ID = 'ID'
    RA = 'RA'
    DEC = 'Dec'
    X = 'x'
    Y = 'y'
    NDATA = 'ndata'
    NDETECT = 'ndetect'
    SIGMA = 'sigma'
    SUMSIGMA = 'sumsigma'
    REDCHI2_OUT = 'redchi2_out'
    SEPMIN = 'sepmin'
    ID_DOPHOT = 'ID_dophot'
    TYPE = 'type'
    MAG = 'mag'
    MAGE = 'mage'
    T0 = 't0'
    TE = 'tE'
    UMIN = 'umin'
    FS = 'fs'
    FB = 'fb'
    T0E = 't0e'
    TEE = 'tEe'
    TEE1 = 'tEe1'
    TEE2 = 'tEe2'
    UMINE = 'umine'
    UMINE1 = 'umine1'
    UMINE2 = 'umine2'
    FSE = 'fse'
    FBE = 'fbe'
    CHI2 = 'chi2'
    T0FS = 't0FS'
    TEFS = 'tEFS'
    UMINFS = 'uminFS'
    RHOFS = 'rhoFS'
    FSFS = 'fsFS'
    FBFS = 'fbFS'
    T0EFS = 't0eFS'
    TEEFS = 'tEeFS'
    TEE1FS = 'tEe1FS'
    TEE2FS = 'tEe2FS'
    UMINEFS = 'umineFS'
    UMINE1FS = 'umine1FS'
    UMINE2FS = 'umine2FS'
    RHOEFS = 'rhoeFS'
    RHOE1FS = 'rhoe1FS'
    RHOE2FS = 'rhoe2FS'
    FSEFS = 'fseFS'
    FBEFS = 'fbeFS'
    CHI2FS = 'chi2FS'


class TakahiroSumiCandlistFileColumnName(StrEnum):
    FIELD = 'field'
    CHIP = 'chip'
    NSUB = 'nsub'
    ID = 'ID'
    TAG = 'tag'
    X = 'x'
    Y = 'y'


class TakahiroSumiCandlistAlertIdFileColumnName(StrEnum):
    FIELD = 'field'
    CHIP = 'chip'
    SUBFRAME = 'subframe'
    ID = 'ID'
    TAG = 'tag'
    X__PIXELS = 'x__pixels'
    Y__PIXELS = 'y__pixels'
    RA_J2000 = 'RA_j2000'
    DEC_J2000 = 'Dec_j2000'
    SEPARATION_TO_ALERT_POSITION0__PIXELS = 'separation_to_alert_position0__pixels'
    ALERT_ID0 = 'alert_id0'
    ALERT_X0__PIXELS = 'alert_x0__pixels'
    ALERT_Y0__PIXELS = 'alert_y0__pixels'
    SEPARATION_TO_ALERT_POSITION1__PIXELS = 'separation_to_alert_position1__pixels'
    ALERT_ID1 = 'alert_id1'
    ALERT_X1__PIXELS = 'alert_x1__pixels'
    ALERT_Y1__PIXELS = 'alert_y1__pixels'
