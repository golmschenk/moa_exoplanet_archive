from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

import pandas as pd
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
    ra_j2000: float
    dec_j2000: float
    candidate_metadata: Optional[CandidateMetadata] = None
    alert_metadata: Optional[List[AlertMetadata]] = None


class CandlistRadecFileColumnName(StrEnum):
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


class CandlistFileColumnName(StrEnum):
    FIELD = 'field'
    CHIP = 'chip'
    NSUB = 'nsub'
    ID = 'ID'
    TAG = 'tag'
    X = 'x'
    Y = 'y'

class CandlistAlertIdFileColumnName(StrEnum):
    FIELD = 'field'
    CHIP = 'chip'
    SUBFRAME = 'subframe'
    ID = 'ID'
    TAG = 'tag'
    X__PIXELS = 'x__pixels'
    Y__PIXELS = 'y__pixels'
    RA_J2000 = 'RA_j2000'
    DEC_J2000 = 'Dec_j2000'
    SEPARATION_TO_ALERT_POSITION__PIXELS = 'separation_to_alert_position__pixels'
    ALERT_ID = 'alert_id'
    ALERT_X__PIXELS = 'alert_x__pixels'
    ALERT_Y__PIXELS = 'alert_y__pixels'


candlist_ra_dec_path = Path('metadata/candlist_RADec.dat.20230202.txt')
candlist_path = Path('metadata/candlist.txt')
candlist_alert_id_path = Path('metadata/candlist_AlertID.dat.20230202.txt')


def process_metadata():
    candlist_radec_data_frame = pd.read_csv(candlist_ra_dec_path, delim_whitespace=True, skipinitialspace=True,
                                            names=[name.value for name in CandlistRadecFileColumnName], comment='#')
    candlist_data_frame = pd.read_csv(candlist_path, delim_whitespace=True, skipinitialspace=True, comment='#',
                                      names=[name.value for name in CandlistFileColumnName],
                                      usecols=[0, 2, 3, 4, 5, 6, 7])
    pass


if __name__ == '__main__':
    process_metadata()
