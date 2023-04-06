from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict

import numpy as np
from astropy import units
from astropy.coordinates import Angle

try:
    from enum import StrEnum
except ImportError:
    from backports.strenum import StrEnum


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


metadata_column_description_dictionary = {
    MetadataColumnName.FIELD: 'Galactic bulge observing field number.',
    MetadataColumnName.CHIP: 'Chip number.',
    MetadataColumnName.SUBFRAME: 'Subframe number.',
    MetadataColumnName.ID: 'Target identifier.',
    MetadataColumnName.TAG: 'Target classification.',
    MetadataColumnName.X: 'Pixel position x coordinate.',
    MetadataColumnName.Y: 'Pixel position y coordinate.',
    MetadataColumnName.RA_J2000: 'Right ascension J2000.',
    MetadataColumnName.DEC_J2000: 'Declination J2000.',
    MetadataColumnName.NUMBER_OF_DATA_POINTS: 'Number of data points in light curve.',
    MetadataColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED: 'Number of frames in which the object is detected.',
    MetadataColumnName.MAX_SIGNIFICANCE: 'Significance at maximum significant point in light curve.',
    MetadataColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE: 'Sum of significance of continuous significant points in light curve.',
    MetadataColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW: 'Chi square outside of the search window in light curve.',
    MetadataColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION: 'Separation in pixels to the closest DoPHOT object in the reference image.',
    MetadataColumnName.DOPHOT_ID: 'DoPHOT identifier.',
    MetadataColumnName.DOPHOT_TYPE: 'DoPHOT type.',
    MetadataColumnName.DOPHOT_MAGNITUDE: 'DoPHOT magnitude.',
    MetadataColumnName.DOPHOT_MAGNITUDE_ERROR: 'DoPHOT magnitude error.',
    MetadataColumnName.PSPL_T0: 'The time of minimum angular separation between lens and source in the fitted PSPL model.',
    MetadataColumnName.PSPL_TE: 'The Einstein crossing time in the fitted PSPL model.',
    MetadataColumnName.PSPL_U0: 'The minimum angular separation between lens and source normalized by the angular Einstein radius in the fitted PSPL model.',
    MetadataColumnName.PSPL_SOURCE_FLUX: 'The source flux in the fitted PSPL model.',
    MetadataColumnName.PSPL_BLENDING_FLUX: 'The blending flux in the fitted PSPL model.',
    MetadataColumnName.PSPL_T0_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.PSPL_T0}`.',
    MetadataColumnName.PSPL_TE_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.PSPL_TE}`.',
    MetadataColumnName.PSPL_TE_ERROR_LOWER_LIMIT: f'Lower limit error in `{MetadataColumnName.PSPL_TE}`.',
    MetadataColumnName.PSPL_TE_ERROR_UPPER_LIMIT: f'Upper limit error in `{MetadataColumnName.PSPL_TE}`.',
    MetadataColumnName.PSPL_U0_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.PSPL_U0}`.',
    MetadataColumnName.PSPL_U0_ERROR_LOWER_LIMIT: f'Lower limit error in `{MetadataColumnName.PSPL_U0}`.',
    MetadataColumnName.PSPL_U0_ERROR_UPPER_LIMIT: f'Upper limit error in `{MetadataColumnName.PSPL_U0}`.',
    MetadataColumnName.PSPL_SOURCE_FLUX_ERROR: f'Error in `{MetadataColumnName.PSPL_SOURCE_FLUX}`.',
    MetadataColumnName.PSPL_BLENDING_FLUX_ERROR: f'Error in `{MetadataColumnName.PSPL_BLENDING_FLUX}`.',
    MetadataColumnName.PSPL_CHI_SQUARED: 'The chi squared statistic of the fitted PSPL model.',
    MetadataColumnName.FSPL_T0: 'The time of minimum angular separation between lens and source in the fitted FSPL model.',
    MetadataColumnName.FSPL_TE: 'The Einstein crossing time in the fitted FSPL model.',
    MetadataColumnName.FSPL_U0: 'The minimum angular separation between lens and source normalized by the angular Einstein radius in the fitted FSPL model.',
    MetadataColumnName.FSPL_RHO: 'The angular source size normalized by the angular Einstein radius in the fitted FSPL model.',
    MetadataColumnName.FSPL_SOURCE_FLUX: 'The source flux in the fitted FSPL model.',
    MetadataColumnName.FSPL_BLENDING_FLUX: 'The blending flux in the fitted FSPL model.',
    MetadataColumnName.FSPL_T0_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.FSPL_T0}`.',
    MetadataColumnName.FSPL_TE_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.FSPL_TE}`.',
    MetadataColumnName.FSPL_TE_ERROR_LOWER_LIMIT: f'Lower limit error in `{MetadataColumnName.FSPL_TE}`.',
    MetadataColumnName.FSPL_TE_ERROR_UPPER_LIMIT: f'Upper limit error in `{MetadataColumnName.FSPL_TE}`',
    MetadataColumnName.FSPL_U0_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.FSPL_U0}`.',
    MetadataColumnName.FSPL_U0_ERROR_LOWER_LIMIT: f'Lower limit error in `{MetadataColumnName.FSPL_U0}`.',
    MetadataColumnName.FSPL_U0_ERROR_UPPER_LIMIT: f'Upper limit error in `{MetadataColumnName.FSPL_U0}`.',
    MetadataColumnName.FSPL_RHO_PARABOLIC_ERROR: f'Parabolic error in `{MetadataColumnName.FSPL_RHO}`.',
    MetadataColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: f'Lower limit error in `{MetadataColumnName.FSPL_RHO}`.',
    MetadataColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: f'Upper limit error in `{MetadataColumnName.FSPL_RHO}`.',
    MetadataColumnName.FSPL_SOURCE_FLUX_ERROR: f'Error in `{MetadataColumnName.FSPL_SOURCE_FLUX}`.',
    MetadataColumnName.FSPL_BLENDING_FLUX_ERROR: f'Error in `{MetadataColumnName.FSPL_BLENDING_FLUX}`.',
    MetadataColumnName.FSPL_CHI_SQUARED: 'The chi squared statistic of the fitted FSPL model.',
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION0: 'Separation to alert 1 position.',
    MetadataColumnName.ALERT_ID0: 'Alert 1 identifier.',
    MetadataColumnName.ALERT_X0: 'Alert 1 x pixel position.',
    MetadataColumnName.ALERT_Y0: 'Alert 1 y pixel position.',
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION1: 'Separation to alert 2 position.',
    MetadataColumnName.ALERT_ID1: 'Alert 2 identifier.',
    MetadataColumnName.ALERT_X1: 'Alert 2 x pixel position.',
    MetadataColumnName.ALERT_Y1: 'Alert 2 y pixel position.',
}

for metadata_column_name in MetadataColumnName:
    assert metadata_column_description_dictionary[metadata_column_name] is not None

metadata_human_readable_short_name_dictionary = {
    MetadataColumnName.FIELD: 'Field',
    MetadataColumnName.CHIP: 'Chip',
    MetadataColumnName.SUBFRAME: 'Subframe',
    MetadataColumnName.ID: 'ID',
    MetadataColumnName.TAG: 'Tag',
    MetadataColumnName.X: 'x',
    MetadataColumnName.Y: 'y',
    MetadataColumnName.RA_J2000: 'RA',
    MetadataColumnName.DEC_J2000: 'Dec',
    MetadataColumnName.NUMBER_OF_DATA_POINTS: 'Number of data points',
    MetadataColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED: 'Frames detected',
    MetadataColumnName.MAX_SIGNIFICANCE: 'Max significance',
    MetadataColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE: 'Sum significance',
    MetadataColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW: 'Chi squared outside window',
    MetadataColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION: 'DoPHOT reference separation',
    MetadataColumnName.DOPHOT_ID: 'DoPHOT ID',
    MetadataColumnName.DOPHOT_TYPE: 'DoPHOT type',
    MetadataColumnName.DOPHOT_MAGNITUDE: 'DoPHOT magnitude',
    MetadataColumnName.DOPHOT_MAGNITUDE_ERROR: 'DoPHOT magnitude error',
    MetadataColumnName.PSPL_T0: 't_0 (PSPL)',
    MetadataColumnName.PSPL_TE: 't_E (PSPL)',
    MetadataColumnName.PSPL_U0: 'u_0 (PSPL)',
    MetadataColumnName.PSPL_SOURCE_FLUX: 'Source flux (PSPL)',
    MetadataColumnName.PSPL_BLENDING_FLUX: 'Blending flux (PSPL)',
    MetadataColumnName.PSPL_T0_PARABOLIC_ERROR: f't_0 error (PSPL)',
    MetadataColumnName.PSPL_TE_PARABOLIC_ERROR: f't_E error (PSPL)',
    MetadataColumnName.PSPL_TE_ERROR_LOWER_LIMIT: f't_E lower limit (PSPL)',
    MetadataColumnName.PSPL_TE_ERROR_UPPER_LIMIT: f't_E upper limit (PSPL)',
    MetadataColumnName.PSPL_U0_PARABOLIC_ERROR: f'u_0 error (PSPL)',
    MetadataColumnName.PSPL_U0_ERROR_LOWER_LIMIT: f'u_0 lower limit (PSPL)',
    MetadataColumnName.PSPL_U0_ERROR_UPPER_LIMIT: f'u_0 upper limit (PSPL)',
    MetadataColumnName.PSPL_SOURCE_FLUX_ERROR: f'Source flux error (PSPL)',
    MetadataColumnName.PSPL_BLENDING_FLUX_ERROR: f'Blending flux error (PSPL)',
    MetadataColumnName.PSPL_CHI_SQUARED: 'Chi squared (PSPL)',
    MetadataColumnName.FSPL_T0: 't_0 (FSPL)',
    MetadataColumnName.FSPL_TE: 't_E (FSPL)',
    MetadataColumnName.FSPL_U0: 'u_0 (FSPL)',
    MetadataColumnName.FSPL_RHO: 'rho (FSPL)',
    MetadataColumnName.FSPL_SOURCE_FLUX: 'Source flux (FSPL)',
    MetadataColumnName.FSPL_BLENDING_FLUX: 'Blending flux (FSPL)',
    MetadataColumnName.FSPL_T0_PARABOLIC_ERROR: f't_0 error (FSPL)',
    MetadataColumnName.FSPL_TE_PARABOLIC_ERROR: f't_E error (FSPL)',
    MetadataColumnName.FSPL_TE_ERROR_LOWER_LIMIT: f't_E lower limit (FSPL)',
    MetadataColumnName.FSPL_TE_ERROR_UPPER_LIMIT: f't_E upper limit (FSPL)',
    MetadataColumnName.FSPL_U0_PARABOLIC_ERROR: f'u_0 error (FSPL)',
    MetadataColumnName.FSPL_U0_ERROR_LOWER_LIMIT: f'u_0 lower limit (FSPL)',
    MetadataColumnName.FSPL_U0_ERROR_UPPER_LIMIT: f'u_0 upper limit (FSPL)',
    MetadataColumnName.FSPL_RHO_PARABOLIC_ERROR: f'u_0 error (FSPL)',
    MetadataColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: f'u_0 lower limit (FSPL)',
    MetadataColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: f'u_0 upper limit (FSPL)',
    MetadataColumnName.FSPL_SOURCE_FLUX_ERROR: f'Source flux error (FSPL)',
    MetadataColumnName.FSPL_BLENDING_FLUX_ERROR: f'Blending flux error (FSPL)',
    MetadataColumnName.FSPL_CHI_SQUARED: 'Chi squared (FSPL)',
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION0: 'Alert 1 separation',
    MetadataColumnName.ALERT_ID0: 'Alert 1 ID',
    MetadataColumnName.ALERT_X0: 'Alert 1 x',
    MetadataColumnName.ALERT_Y0: 'Alert 1 y',
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION1: 'Alert 2 separation',
    MetadataColumnName.ALERT_ID1: 'Alert 2 ID',
    MetadataColumnName.ALERT_X1: 'Alert 2 x',
    MetadataColumnName.ALERT_Y1: 'Alert 2 y',
}

for metadata_column_name in MetadataColumnName:
    assert metadata_human_readable_short_name_dictionary[metadata_column_name] is not None


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

metadata_units_dictionary: Dict[str, units.Unit] = {
    MetadataColumnName.FIELD: units.dimensionless_unscaled,
    MetadataColumnName.CHIP: units.dimensionless_unscaled,
    MetadataColumnName.SUBFRAME: units.dimensionless_unscaled,
    MetadataColumnName.ID: units.dimensionless_unscaled,
    MetadataColumnName.TAG: units.dimensionless_unscaled,
    MetadataColumnName.X: units.pixel,
    MetadataColumnName.Y: units.pixel,
    MetadataColumnName.RA_J2000: units.degree,
    MetadataColumnName.DEC_J2000: units.degree,
    MetadataColumnName.NUMBER_OF_DATA_POINTS: units.dimensionless_unscaled,
    MetadataColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED: units.dimensionless_unscaled,
    MetadataColumnName.MAX_SIGNIFICANCE: units.dimensionless_unscaled,
    MetadataColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE: units.dimensionless_unscaled,
    MetadataColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW: units.dimensionless_unscaled,
    MetadataColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION: units.pixel,
    MetadataColumnName.DOPHOT_ID: units.dimensionless_unscaled,
    MetadataColumnName.DOPHOT_TYPE: units.dimensionless_unscaled,
    MetadataColumnName.DOPHOT_MAGNITUDE: units.dimensionless_unscaled,
    MetadataColumnName.DOPHOT_MAGNITUDE_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_T0: units.day,
    MetadataColumnName.PSPL_TE: units.day,
    MetadataColumnName.PSPL_U0: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_SOURCE_FLUX: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_BLENDING_FLUX: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_T0_PARABOLIC_ERROR: units.day,
    MetadataColumnName.PSPL_TE_PARABOLIC_ERROR: units.day,
    MetadataColumnName.PSPL_TE_ERROR_LOWER_LIMIT: units.day,
    MetadataColumnName.PSPL_TE_ERROR_UPPER_LIMIT: units.day,
    MetadataColumnName.PSPL_U0_PARABOLIC_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_U0_ERROR_LOWER_LIMIT: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_U0_ERROR_UPPER_LIMIT: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_SOURCE_FLUX_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_BLENDING_FLUX_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.PSPL_CHI_SQUARED: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_T0: units.day,
    MetadataColumnName.FSPL_TE: units.day,
    MetadataColumnName.FSPL_U0: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_RHO: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_SOURCE_FLUX: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_BLENDING_FLUX: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_T0_PARABOLIC_ERROR: units.day,
    MetadataColumnName.FSPL_TE_PARABOLIC_ERROR: units.day,
    MetadataColumnName.FSPL_TE_ERROR_LOWER_LIMIT: units.day,
    MetadataColumnName.FSPL_TE_ERROR_UPPER_LIMIT: units.day,
    MetadataColumnName.FSPL_U0_PARABOLIC_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_U0_ERROR_LOWER_LIMIT: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_U0_ERROR_UPPER_LIMIT: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_RHO_PARABOLIC_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_SOURCE_FLUX_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_BLENDING_FLUX_ERROR: units.dimensionless_unscaled,
    MetadataColumnName.FSPL_CHI_SQUARED: units.dimensionless_unscaled,
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION0: units.pixel,
    MetadataColumnName.ALERT_ID0: units.dimensionless_unscaled,
    MetadataColumnName.ALERT_X0: units.pixel,
    MetadataColumnName.ALERT_Y0: units.pixel,
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION1: units.pixel,
    MetadataColumnName.ALERT_ID1: units.dimensionless_unscaled,
    MetadataColumnName.ALERT_X1: units.pixel,
    MetadataColumnName.ALERT_Y1: units.pixel,
}

metadata_types_dictionary: Dict[str, np.dtype] = {
    MetadataColumnName.FIELD: np.int32,
    MetadataColumnName.CHIP: np.int32,
    MetadataColumnName.SUBFRAME: np.int32,
    MetadataColumnName.ID: np.int32,
    MetadataColumnName.TAG: str,
    MetadataColumnName.X: np.float64,
    MetadataColumnName.Y: np.float64,
    MetadataColumnName.RA_J2000: np.float64,
    MetadataColumnName.DEC_J2000: np.float64,
    MetadataColumnName.NUMBER_OF_DATA_POINTS: np.int32,
    MetadataColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED: np.int32,
    MetadataColumnName.MAX_SIGNIFICANCE: np.float64,
    MetadataColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE: np.float64,
    MetadataColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW: np.float64,
    MetadataColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION: np.float64,
    MetadataColumnName.DOPHOT_ID: np.int32,
    MetadataColumnName.DOPHOT_TYPE: np.int32,
    MetadataColumnName.DOPHOT_MAGNITUDE: np.float64,
    MetadataColumnName.DOPHOT_MAGNITUDE_ERROR: np.float64,
    MetadataColumnName.PSPL_T0: np.float64,
    MetadataColumnName.PSPL_TE: np.float64,
    MetadataColumnName.PSPL_U0: np.float64,
    MetadataColumnName.PSPL_SOURCE_FLUX: np.float64,
    MetadataColumnName.PSPL_BLENDING_FLUX: np.float64,
    MetadataColumnName.PSPL_T0_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.PSPL_TE_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.PSPL_TE_ERROR_LOWER_LIMIT: np.float64,
    MetadataColumnName.PSPL_TE_ERROR_UPPER_LIMIT: np.float64,
    MetadataColumnName.PSPL_U0_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.PSPL_U0_ERROR_LOWER_LIMIT: np.float64,
    MetadataColumnName.PSPL_U0_ERROR_UPPER_LIMIT: np.float64,
    MetadataColumnName.PSPL_SOURCE_FLUX_ERROR: np.float64,
    MetadataColumnName.PSPL_BLENDING_FLUX_ERROR: np.float64,
    MetadataColumnName.PSPL_CHI_SQUARED: np.float64,
    MetadataColumnName.FSPL_T0: np.float64,
    MetadataColumnName.FSPL_TE: np.float64,
    MetadataColumnName.FSPL_U0: np.float64,
    MetadataColumnName.FSPL_RHO: np.float64,
    MetadataColumnName.FSPL_SOURCE_FLUX: np.float64,
    MetadataColumnName.FSPL_BLENDING_FLUX: np.float64,
    MetadataColumnName.FSPL_T0_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.FSPL_TE_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.FSPL_TE_ERROR_LOWER_LIMIT: np.float64,
    MetadataColumnName.FSPL_TE_ERROR_UPPER_LIMIT: np.float64,
    MetadataColumnName.FSPL_U0_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.FSPL_U0_ERROR_LOWER_LIMIT: np.float64,
    MetadataColumnName.FSPL_U0_ERROR_UPPER_LIMIT: np.float64,
    MetadataColumnName.FSPL_RHO_PARABOLIC_ERROR: np.float64,
    MetadataColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: np.float64,
    MetadataColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: np.float64,
    MetadataColumnName.FSPL_SOURCE_FLUX_ERROR: np.float64,
    MetadataColumnName.FSPL_BLENDING_FLUX_ERROR: np.float64,
    MetadataColumnName.FSPL_CHI_SQUARED: np.float64,
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION0: np.float64,
    MetadataColumnName.ALERT_ID0: str,
    MetadataColumnName.ALERT_X0: np.float64,
    MetadataColumnName.ALERT_Y0: np.float64,
    MetadataColumnName.SEPARATION_TO_ALERT_POSITION1: np.float64,
    MetadataColumnName.ALERT_ID1: str,
    MetadataColumnName.ALERT_X1: np.float64,
    MetadataColumnName.ALERT_Y1: np.float64,
}


if __name__ == '__main__':
    with Path('column_descriptions.csv').open('w') as column_descriptions_file:
        column_descriptions_file.write('column_name,column_description\n')
        print()
        for column_name, column_description in metadata_column_description_dictionary.items():
            column_descriptions_file.write(f'{column_name},"{column_description}"\n')
            print(f'{column_name}: {column_description}')