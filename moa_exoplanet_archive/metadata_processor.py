import numpy as np
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

import pandas as pd
from astropy import units
from astropy.coordinates import Angle
from backports.strenum import StrEnum


@dataclass
class AlertMetadata:
    separation_to_alert_position0__pixels: float
    alert_id0: str
    alert_x0__pixels: float
    alert_y0__pixels: float
    separation_to_alert_position1__pixels: float
    alert_id1: str
    alert_x1__pixels: float
    alert_y1__pixels: float


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
    alert_metadata: Optional[AlertMetadata] = None


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
        metadata_list: List[TargetMetadata] = []
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
            alert_metadata = self.get_alert_metadata_for_light_curve(light_curve_path)
            target_metadata = TargetMetadata(field=field, chip=chip, subframe=subframe, id=id_, tag=tag,
                                             x__pixels=x__pixels, y__pixels=y__pixels, ra_j2000=ra, dec_j2000=dec,
                                             candidate_metadata=candidate_metadata, alert_metadata=alert_metadata)
            metadata_list.append(target_metadata)
        pass

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

    def get_alert_metadata_for_light_curve(self, light_curve_path: Path) -> Optional[AlertMetadata]:
        field, chip, subframe, id_ = self.extract_field_chip_subframe_and_id_from_light_curve_path(light_curve_path)
        try:
            light_curve_row = self.candlist_alert_id_data_frame[
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.FIELD] == f'gb{field}') &
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.CHIP] == chip) &
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.SUBFRAME] == subframe) &
                (self.candlist_alert_id_data_frame[TakahiroSumiCandlistAlertIdFileColumnName.ID] == id_)].iloc[0]
        except IndexError:
            return None
        alert_metadata = AlertMetadata(
            separation_to_alert_position0__pixels=light_curve_row[
                TakahiroSumiCandlistAlertIdFileColumnName.SEPARATION_TO_ALERT_POSITION0__PIXELS],
            alert_id0=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_ID0],
            alert_x0__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_X0__PIXELS],
            alert_y0__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_Y0__PIXELS],
            separation_to_alert_position1__pixels=light_curve_row[
                TakahiroSumiCandlistAlertIdFileColumnName.SEPARATION_TO_ALERT_POSITION1__PIXELS],
            alert_id1=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_ID1],
            alert_x1__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_X1__PIXELS],
            alert_y1__pixels=light_curve_row[TakahiroSumiCandlistAlertIdFileColumnName.ALERT_Y1__PIXELS],
        )
        return alert_metadata


if __name__ == '__main__':
    MetadataProcessor().process_metadata()
