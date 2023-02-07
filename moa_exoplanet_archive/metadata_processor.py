import astropy.io.ascii
import numpy as np
import re
import subprocess
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, List, Dict, Any

import pandas as pd
from astropy import units
from astropy.coordinates import Angle
from astropy.table import Table
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


class ExoplanetArchiveColumnName(StrEnum):
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


ExoplanetArchiveColumnNameToUnitDictionary: Dict[str, units.Unit] = {
    ExoplanetArchiveColumnName.FIELD: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.CHIP: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.SUBFRAME: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.ID: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.TAG: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.X: units.pixel,
    ExoplanetArchiveColumnName.Y: units.pixel,
    ExoplanetArchiveColumnName.RA_J2000: units.degree,
    ExoplanetArchiveColumnName.DEC_J2000: units.degree,
    ExoplanetArchiveColumnName.NUMBER_OF_DATA_POINTS: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.MAX_SIGNIFICANCE: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION: units.pixel,
    ExoplanetArchiveColumnName.DOPHOT_ID: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.DOPHOT_TYPE: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.DOPHOT_MAGNITUDE: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.DOPHOT_MAGNITUDE_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_T0: units.day,
    ExoplanetArchiveColumnName.PSPL_TE: units.day,
    ExoplanetArchiveColumnName.PSPL_U0: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_SOURCE_FLUX: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_BLENDING_FLUX: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_T0_PARABOLIC_ERROR: units.day,
    ExoplanetArchiveColumnName.PSPL_TE_PARABOLIC_ERROR: units.day,
    ExoplanetArchiveColumnName.PSPL_TE_ERROR_LOWER_LIMIT: units.day,
    ExoplanetArchiveColumnName.PSPL_TE_ERROR_UPPER_LIMIT: units.day,
    ExoplanetArchiveColumnName.PSPL_U0_PARABOLIC_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_U0_ERROR_LOWER_LIMIT: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_U0_ERROR_UPPER_LIMIT: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_SOURCE_FLUX_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_BLENDING_FLUX_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.PSPL_CHI_SQUARED: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_T0: units.day,
    ExoplanetArchiveColumnName.FSPL_TE: units.day,
    ExoplanetArchiveColumnName.FSPL_U0: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_RHO: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_SOURCE_FLUX: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_BLENDING_FLUX: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_T0_PARABOLIC_ERROR: units.day,
    ExoplanetArchiveColumnName.FSPL_TE_PARABOLIC_ERROR: units.day,
    ExoplanetArchiveColumnName.FSPL_TE_ERROR_LOWER_LIMIT: units.day,
    ExoplanetArchiveColumnName.FSPL_TE_ERROR_UPPER_LIMIT: units.day,
    ExoplanetArchiveColumnName.FSPL_U0_PARABOLIC_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_U0_ERROR_LOWER_LIMIT: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_U0_ERROR_UPPER_LIMIT: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_RHO_PARABOLIC_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_SOURCE_FLUX_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_BLENDING_FLUX_ERROR: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.FSPL_CHI_SQUARED: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.SEPARATION_TO_ALERT_POSITION0: units.pixel,
    ExoplanetArchiveColumnName.ALERT_ID0: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.ALERT_X0: units.pixel,
    ExoplanetArchiveColumnName.ALERT_Y0: units.pixel,
    ExoplanetArchiveColumnName.SEPARATION_TO_ALERT_POSITION1: units.pixel,
    ExoplanetArchiveColumnName.ALERT_ID1: units.dimensionless_unscaled,
    ExoplanetArchiveColumnName.ALERT_X1: units.pixel,
    ExoplanetArchiveColumnName.ALERT_Y1: units.pixel,
}

ExoplanetArchiveColumnNameToTypeDictionary: Dict[str, np.dtype] = {
    ExoplanetArchiveColumnName.FIELD: np.int32,
    ExoplanetArchiveColumnName.CHIP: np.int32,
    ExoplanetArchiveColumnName.SUBFRAME: np.int32,
    ExoplanetArchiveColumnName.ID: np.int32,
    ExoplanetArchiveColumnName.TAG: str,
    ExoplanetArchiveColumnName.X: np.float64,
    ExoplanetArchiveColumnName.Y: np.float64,
    ExoplanetArchiveColumnName.RA_J2000: np.float64,
    ExoplanetArchiveColumnName.DEC_J2000: np.float64,
    ExoplanetArchiveColumnName.NUMBER_OF_DATA_POINTS: np.int32,
    ExoplanetArchiveColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED: np.int32,
    ExoplanetArchiveColumnName.MAX_SIGNIFICANCE: np.float64,
    ExoplanetArchiveColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE: np.float64,
    ExoplanetArchiveColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW: np.float64,
    ExoplanetArchiveColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION: np.float64,
    ExoplanetArchiveColumnName.DOPHOT_ID: np.int32,
    ExoplanetArchiveColumnName.DOPHOT_TYPE: np.int32,
    ExoplanetArchiveColumnName.DOPHOT_MAGNITUDE: np.float64,
    ExoplanetArchiveColumnName.DOPHOT_MAGNITUDE_ERROR: np.float64,
    ExoplanetArchiveColumnName.PSPL_T0: np.float64,
    ExoplanetArchiveColumnName.PSPL_TE: np.float64,
    ExoplanetArchiveColumnName.PSPL_U0: np.float64,
    ExoplanetArchiveColumnName.PSPL_SOURCE_FLUX: np.float64,
    ExoplanetArchiveColumnName.PSPL_BLENDING_FLUX: np.float64,
    ExoplanetArchiveColumnName.PSPL_T0_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.PSPL_TE_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.PSPL_TE_ERROR_LOWER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.PSPL_TE_ERROR_UPPER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.PSPL_U0_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.PSPL_U0_ERROR_LOWER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.PSPL_U0_ERROR_UPPER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.PSPL_SOURCE_FLUX_ERROR: np.float64,
    ExoplanetArchiveColumnName.PSPL_BLENDING_FLUX_ERROR: np.float64,
    ExoplanetArchiveColumnName.PSPL_CHI_SQUARED: np.float64,
    ExoplanetArchiveColumnName.FSPL_T0: np.float64,
    ExoplanetArchiveColumnName.FSPL_TE: np.float64,
    ExoplanetArchiveColumnName.FSPL_U0: np.float64,
    ExoplanetArchiveColumnName.FSPL_RHO: np.float64,
    ExoplanetArchiveColumnName.FSPL_SOURCE_FLUX: np.float64,
    ExoplanetArchiveColumnName.FSPL_BLENDING_FLUX: np.float64,
    ExoplanetArchiveColumnName.FSPL_T0_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.FSPL_TE_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.FSPL_TE_ERROR_LOWER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.FSPL_TE_ERROR_UPPER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.FSPL_U0_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.FSPL_U0_ERROR_LOWER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.FSPL_U0_ERROR_UPPER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.FSPL_RHO_PARABOLIC_ERROR: np.float64,
    ExoplanetArchiveColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: np.float64,
    ExoplanetArchiveColumnName.FSPL_SOURCE_FLUX_ERROR: np.float64,
    ExoplanetArchiveColumnName.FSPL_BLENDING_FLUX_ERROR: np.float64,
    ExoplanetArchiveColumnName.FSPL_CHI_SQUARED: np.float64,
    ExoplanetArchiveColumnName.SEPARATION_TO_ALERT_POSITION0: np.float64,
    ExoplanetArchiveColumnName.ALERT_ID0: str,
    ExoplanetArchiveColumnName.ALERT_X0: np.float64,
    ExoplanetArchiveColumnName.ALERT_Y0: np.float64,
    ExoplanetArchiveColumnName.SEPARATION_TO_ALERT_POSITION1: np.float64,
    ExoplanetArchiveColumnName.ALERT_ID1: str,
    ExoplanetArchiveColumnName.ALERT_X1: np.float64,
    ExoplanetArchiveColumnName.ALERT_Y1: np.float64,
}


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


candlist_ra_dec_path = Path('metadata/candlist_RADec.dat.20230202.txt')
candlist_path = Path('metadata/candlist.txt')
candlist_alert_id_path = Path('metadata/candlist_AlertID.dat.20230202.txt')


class MetadataProcessor:
    def __init__(self):
        self.candlist_ra_dec_data_frame = pd.read_csv(candlist_ra_dec_path, delim_whitespace=True,
                                                      skipinitialspace=True, comment='#',
                                                      names=[name for name in TakahiroSumiCandlistRaDecFileColumnName])
        self.candlist_data_frame = pd.read_csv(candlist_path, delim_whitespace=True, skipinitialspace=True, comment='#',
                                               names=[name for name in TakahiroSumiCandlistFileColumnName],
                                               usecols=[0, 2, 3, 4, 5, 6, 7])
        self.candlist_alert_id_data_frame = pd.read_csv(
            candlist_alert_id_path, delim_whitespace=True, skipinitialspace=True, comment='#',
            names=[name for name in TakahiroSumiCandlistAlertIdFileColumnName],
        )

    def get_ra_and_dec_for_ccd_position(self, field: int, chip: int, x__pixels: float, y__pixels: float
                                        ) -> (Angle, Angle):
        ra_dec_perl_script_directory = Path('RADEC')
        result = subprocess.run(['./ccd2sky.pl', str(field), str(chip), str(x__pixels), str(y__pixels)],
                                cwd=ra_dec_perl_script_directory, capture_output=True)
        result_string = result.stdout.decode('utf-8')
        match = re.match(r'RA=([^\n]*)Dec=([^\n]*)\n', result_string)
        ra_hms_string = match.group(1).strip()
        dec_dms_string = match.group(2).strip()
        ra__hour_angles = Angle(ra_hms_string, unit=units.hourangle)
        ra__degrees = ra__hour_angles.to(units.deg)
        dec__degrees = Angle(dec_dms_string, unit=units.deg)
        return ra__degrees, dec__degrees

    def get_ra_and_dec_for_light_curve_path(self, light_curve_path: Path) -> (Angle, Angle):
        field, chip, subframe, id_ = self.extract_field_chip_subframe_and_id_from_light_curve_path(light_curve_path)
        light_curve_row = self.candlist_data_frame[
            (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.FIELD] == f'gb{field}') &
            (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.CHIP] == chip) &
            (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.NSUB] == subframe) &
            (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.ID] == id_)].iloc[0]
        x__pixels = light_curve_row[TakahiroSumiCandlistFileColumnName.X]
        y__pixels = light_curve_row[TakahiroSumiCandlistFileColumnName.Y]
        return self.get_ra_and_dec_for_ccd_position(field, chip, x__pixels, y__pixels)

    def process_metadata(self):
        light_curve_root_directory = Path('light_curve_sample_ipac_format')
        light_curve_glob = light_curve_root_directory.glob('**/*.ipac.gz')
        target_exoplanet_archive_dictionary_list: List[Dict[str, Any]] = []
        for light_curve_path in light_curve_glob:
            field, chip, subframe, id_ = self.extract_field_chip_subframe_and_id_from_light_curve_path(light_curve_path)
            light_curve_row = self.candlist_data_frame[
                (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.FIELD] == f'gb{field}') &
                (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.CHIP] == chip) &
                (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.NSUB] == subframe) &
                (self.candlist_data_frame[TakahiroSumiCandlistFileColumnName.ID] == id_)].iloc[0]
            x__pixels = light_curve_row[TakahiroSumiCandlistFileColumnName.X]
            y__pixels = light_curve_row[TakahiroSumiCandlistFileColumnName.Y]
            ra, dec = self.get_ra_and_dec_for_ccd_position(field=field, chip=chip, x__pixels=x__pixels,
                                                           y__pixels=y__pixels)
            tag = light_curve_row[TakahiroSumiCandlistFileColumnName.TAG]
            candidate_metadata = self.get_candidate_metadata_for_light_curve(light_curve_path)
            alert_metadata_list = self.get_alert_metadata_list_for_light_curve(light_curve_path)
            target_metadata = TargetMetadata(field=field, chip=chip, subframe=subframe, id=id_, tag=tag,
                                             x__pixels=x__pixels, y__pixels=y__pixels, ra_j2000=ra, dec_j2000=dec,
                                             candidate_metadata=candidate_metadata,
                                             alert_metadata_list=alert_metadata_list)
            target_exoplanet_archive_dictionary = self.target_metadata_to_exoplanet_archive_dictionary(target_metadata)
            target_exoplanet_archive_dictionary_list.append(target_exoplanet_archive_dictionary)
        # TODO: Remove this artificial removal of first junk tag.
        for target_exoplanet_archive_dictionary in target_exoplanet_archive_dictionary_list:
            if target_exoplanet_archive_dictionary[ExoplanetArchiveColumnName.TAG] == 'j':
                del target_exoplanet_archive_dictionary[ExoplanetArchiveColumnName.TAG]
                break
        exoplanet_archive_table = Table(data=target_exoplanet_archive_dictionary_list,
                                        names=[name for name in ExoplanetArchiveColumnName],
                                        dtype=ExoplanetArchiveColumnNameToTypeDictionary.values(),
                                        units=ExoplanetArchiveColumnNameToUnitDictionary,
                                        )
        output_metadata_path = Path('metadata.ipac')
        astropy.io.ascii.write(exoplanet_archive_table, output=output_metadata_path, format='ipac', overwrite=True)

    def extract_field_chip_subframe_and_id_from_light_curve_path(self, light_curve_path) -> (int, int, int, int):
        match = re.match(r'gb(\d+)-R-(\d+)-(\d+)-(\d+).ipac.gz', light_curve_path.name)
        if match is None:
            raise ValueError(f'Could not find light curve naming in {light_curve_path}')
        field = int(match.group(1))
        chip = int(match.group(2))
        subframe = int(match.group(3))
        id_ = int(match.group(4))
        return field, chip, subframe, id_

    def get_candidate_metadata_for_light_curve(self, light_curve_path: Path) -> Optional[CandidateMetadata]:
        field, chip, subframe, id_ = self.extract_field_chip_subframe_and_id_from_light_curve_path(light_curve_path)
        try:
            light_curve_row = self.candlist_ra_dec_data_frame[
                (self.candlist_ra_dec_data_frame[TakahiroSumiCandlistRaDecFileColumnName.FIELD] == f'gb{field}') &
                (self.candlist_ra_dec_data_frame[TakahiroSumiCandlistRaDecFileColumnName.CHIP] == chip) &
                (self.candlist_ra_dec_data_frame[TakahiroSumiCandlistRaDecFileColumnName.NSUB] == subframe) &
                (self.candlist_ra_dec_data_frame[TakahiroSumiCandlistRaDecFileColumnName.ID] == id_)].iloc[0]
        except IndexError:
            return None
        candidate_metadata = CandidateMetadata(
            number_of_data_points=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.NDATA],
            number_of_frames_object_is_detected=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.NDETECT],
            max_significance=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.SIGMA],
            sum_of_continuous_significance=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.SUMSIGMA],
            chi_squared_outside_search_window=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.REDCHI2_OUT],
            dophot_object_reference_image_separation__pixels=light_curve_row[
                TakahiroSumiCandlistRaDecFileColumnName.SEPMIN],
            dophot_id=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.ID_DOPHOT],
            dophot_type=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TYPE],
            dophot_magnitude=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.MAG],
            dophot_magnitude_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.MAGE],
            pspl_t0=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.T0],
            pspl_tE=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TE],
            pspl_u0=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMIN],
            pspl_source_flux=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FS],
            pspl_blending_flux=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FB],
            pspl_t0_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.T0E],
            pspl_tE_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEE],
            pspl_tE_error_lower_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEE1],
            pspl_tE_error_upper_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEE2],
            pspl_u0_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINE],
            pspl_u0_error_lower_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINE1],
            pspl_u0_error_upper_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINE2],
            pspl_source_flux_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FSE],
            pspl_blending_flux_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FBE],
            pspl_chi_squared=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.CHI2],
            fspl_t0=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.T0FS],
            fspl_tE=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEFS],
            fspl_u0=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINFS],
            fspl_rho=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.RHOFS],
            fspl_source_flux=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FSFS],
            fspl_blending_flux=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FBFS],
            fspl_t0_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.T0EFS],
            fspl_tE_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEEFS],
            fspl_tE_error_lower_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEE1FS],
            fspl_tE_error_upper_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.TEE2FS],
            fspl_u0_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINEFS],
            fspl_u0_error_lower_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINE1FS],
            fspl_u0_error_upper_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.UMINE2FS],
            fspl_rho_parabolic_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.RHOEFS],
            fspl_rho_error_lower_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.RHOE1FS],
            fspl_rho_error_upper_limit=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.RHOE2FS],
            fspl_source_flux_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FSEFS],
            fspl_blending_flux_error=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.FBEFS],
            fspl_chi_squared=light_curve_row[TakahiroSumiCandlistRaDecFileColumnName.CHI2FS],
        )
        return candidate_metadata

    def get_alert_metadata_list_for_light_curve(self, light_curve_path: Path) -> Optional[List[AlertMetadata]]:
        field, chip, subframe, id_ = self.extract_field_chip_subframe_and_id_from_light_curve_path(light_curve_path)
        try:
            light_curve_row = self.candlist_alert_id_data_frame[
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.FIELD] == f'gb{field}') &
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.CHIP] == chip) &
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.SUBFRAME] == subframe) &
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.ID] == id_)].iloc[0]
        except IndexError:
            return None
        alert_metadata_list: List[AlertMetadata] = []
        if pd.notna(light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_ID0]):
            alert_metadata0 = AlertMetadata(
                separation_to_alert_position__pixels=light_curve_row[
                    TakahiroSumiCandlistAlertIdFileColumnName.SEPARATION_TO_ALERT_POSITION0__PIXELS],
                alert_id=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_ID0],
                alert_x__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_X0__PIXELS],
                alert_y__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_Y0__PIXELS],
            )
            alert_metadata_list.append(alert_metadata0)
        if pd.notna(light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_ID1]):
            alert_metadata1 = AlertMetadata(
                separation_to_alert_position__pixels=light_curve_row[
                    TakahiroSumiCandlistAlertIdFileColumnName.SEPARATION_TO_ALERT_POSITION1__PIXELS],
                alert_id=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_ID1],
                alert_x__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_X1__PIXELS],
                alert_y__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_Y1__PIXELS],
            )
            alert_metadata_list.append(alert_metadata1)
        if len(alert_metadata_list) == 0:
            return None
        return alert_metadata_list

    def target_metadata_to_exoplanet_archive_dictionary(self, target_metadata: TargetMetadata) -> Dict[str, Any]:
        exoplanet_archive_dictionary = {}
        exoplanet_archive_dictionary.update({
            ExoplanetArchiveColumnName.FIELD: target_metadata.field,
            ExoplanetArchiveColumnName.CHIP: target_metadata.chip,
            ExoplanetArchiveColumnName.SUBFRAME: target_metadata.subframe,
            ExoplanetArchiveColumnName.ID: target_metadata.id,
            ExoplanetArchiveColumnName.TAG: target_metadata.tag,
            ExoplanetArchiveColumnName.X: target_metadata.x__pixels,
            ExoplanetArchiveColumnName.Y: target_metadata.y__pixels,
            ExoplanetArchiveColumnName.RA_J2000: target_metadata.ra_j2000.value,
            ExoplanetArchiveColumnName.DEC_J2000: target_metadata.dec_j2000.value,
        })
        candidate_metadata = target_metadata.candidate_metadata
        if candidate_metadata is not None:
            exoplanet_archive_dictionary.update({
                ExoplanetArchiveColumnName.NUMBER_OF_DATA_POINTS: candidate_metadata.number_of_data_points,
                ExoplanetArchiveColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED:
                    candidate_metadata.number_of_frames_object_is_detected,
                ExoplanetArchiveColumnName.MAX_SIGNIFICANCE: candidate_metadata.max_significance,
                ExoplanetArchiveColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE:
                    candidate_metadata.sum_of_continuous_significance,
                ExoplanetArchiveColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW:
                    candidate_metadata.chi_squared_outside_search_window,
                ExoplanetArchiveColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION:
                    candidate_metadata.dophot_object_reference_image_separation__pixels,
                ExoplanetArchiveColumnName.DOPHOT_ID: candidate_metadata.dophot_id,
                ExoplanetArchiveColumnName.DOPHOT_TYPE: candidate_metadata.dophot_type,
                ExoplanetArchiveColumnName.DOPHOT_MAGNITUDE: candidate_metadata.dophot_magnitude,
                ExoplanetArchiveColumnName.DOPHOT_MAGNITUDE_ERROR: candidate_metadata.dophot_magnitude_error,
                ExoplanetArchiveColumnName.PSPL_T0: candidate_metadata.pspl_t0,
                ExoplanetArchiveColumnName.PSPL_TE: candidate_metadata.pspl_tE,
                ExoplanetArchiveColumnName.PSPL_U0: candidate_metadata.pspl_u0,
                ExoplanetArchiveColumnName.PSPL_SOURCE_FLUX: candidate_metadata.pspl_source_flux,
                ExoplanetArchiveColumnName.PSPL_BLENDING_FLUX: candidate_metadata.pspl_blending_flux,
                ExoplanetArchiveColumnName.PSPL_T0_PARABOLIC_ERROR: candidate_metadata.pspl_t0_parabolic_error,
                ExoplanetArchiveColumnName.PSPL_TE_PARABOLIC_ERROR: candidate_metadata.pspl_tE_parabolic_error,
                ExoplanetArchiveColumnName.PSPL_TE_ERROR_LOWER_LIMIT: candidate_metadata.pspl_tE_error_lower_limit,
                ExoplanetArchiveColumnName.PSPL_TE_ERROR_UPPER_LIMIT: candidate_metadata.pspl_tE_error_upper_limit,
                ExoplanetArchiveColumnName.PSPL_U0_PARABOLIC_ERROR: candidate_metadata.pspl_u0_parabolic_error,
                ExoplanetArchiveColumnName.PSPL_U0_ERROR_LOWER_LIMIT: candidate_metadata.pspl_u0_error_lower_limit,
                ExoplanetArchiveColumnName.PSPL_U0_ERROR_UPPER_LIMIT: candidate_metadata.pspl_u0_error_upper_limit,
                ExoplanetArchiveColumnName.PSPL_SOURCE_FLUX_ERROR: candidate_metadata.pspl_source_flux_error,
                ExoplanetArchiveColumnName.PSPL_BLENDING_FLUX_ERROR: candidate_metadata.pspl_blending_flux_error,
                ExoplanetArchiveColumnName.PSPL_CHI_SQUARED: candidate_metadata.pspl_chi_squared,
                ExoplanetArchiveColumnName.FSPL_T0: candidate_metadata.fspl_t0,
                ExoplanetArchiveColumnName.FSPL_TE: candidate_metadata.fspl_tE,
                ExoplanetArchiveColumnName.FSPL_U0: candidate_metadata.fspl_u0,
                ExoplanetArchiveColumnName.FSPL_RHO: candidate_metadata.fspl_rho,
                ExoplanetArchiveColumnName.FSPL_SOURCE_FLUX: candidate_metadata.fspl_source_flux,
                ExoplanetArchiveColumnName.FSPL_BLENDING_FLUX: candidate_metadata.fspl_blending_flux,
                ExoplanetArchiveColumnName.FSPL_T0_PARABOLIC_ERROR: candidate_metadata.fspl_t0_parabolic_error,
                ExoplanetArchiveColumnName.FSPL_TE_PARABOLIC_ERROR: candidate_metadata.fspl_tE_parabolic_error,
                ExoplanetArchiveColumnName.FSPL_TE_ERROR_LOWER_LIMIT: candidate_metadata.fspl_tE_error_lower_limit,
                ExoplanetArchiveColumnName.FSPL_TE_ERROR_UPPER_LIMIT: candidate_metadata.fspl_tE_error_upper_limit,
                ExoplanetArchiveColumnName.FSPL_U0_PARABOLIC_ERROR: candidate_metadata.fspl_u0_parabolic_error,
                ExoplanetArchiveColumnName.FSPL_U0_ERROR_LOWER_LIMIT: candidate_metadata.fspl_u0_error_lower_limit,
                ExoplanetArchiveColumnName.FSPL_U0_ERROR_UPPER_LIMIT: candidate_metadata.fspl_u0_error_upper_limit,
                ExoplanetArchiveColumnName.FSPL_RHO_PARABOLIC_ERROR: candidate_metadata.fspl_rho_parabolic_error,
                ExoplanetArchiveColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: candidate_metadata.fspl_rho_error_lower_limit,
                ExoplanetArchiveColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: candidate_metadata.fspl_rho_error_upper_limit,
                ExoplanetArchiveColumnName.FSPL_SOURCE_FLUX_ERROR: candidate_metadata.fspl_source_flux_error,
                ExoplanetArchiveColumnName.FSPL_BLENDING_FLUX_ERROR: candidate_metadata.fspl_blending_flux_error,
                ExoplanetArchiveColumnName.FSPL_CHI_SQUARED: candidate_metadata.fspl_chi_squared,
            })
        alert_metadata_list = target_metadata.alert_metadata_list
        if alert_metadata_list is not None:
            alert_metadata0 = alert_metadata_list[0]
            exoplanet_archive_dictionary.update({
                ExoplanetArchiveColumnName.SEPARATION_TO_ALERT_POSITION0:
                    alert_metadata0.separation_to_alert_position__pixels,
                ExoplanetArchiveColumnName.ALERT_ID0: alert_metadata0.alert_id,
                ExoplanetArchiveColumnName.ALERT_X0: alert_metadata0.alert_x__pixels,
                ExoplanetArchiveColumnName.ALERT_Y0: alert_metadata0.alert_y__pixels,
            })
            if len(alert_metadata_list) > 1:
                alert_metadata1 = alert_metadata_list[1]
                exoplanet_archive_dictionary.update({
                    ExoplanetArchiveColumnName.SEPARATION_TO_ALERT_POSITION1:
                        alert_metadata1.separation_to_alert_position__pixels,
                    ExoplanetArchiveColumnName.ALERT_ID1: alert_metadata1.alert_id,
                    ExoplanetArchiveColumnName.ALERT_X1: alert_metadata1.alert_x__pixels,
                    ExoplanetArchiveColumnName.ALERT_Y1: alert_metadata1.alert_y__pixels,
                })
        return exoplanet_archive_dictionary


if __name__ == '__main__':
    MetadataProcessor().process_metadata()
