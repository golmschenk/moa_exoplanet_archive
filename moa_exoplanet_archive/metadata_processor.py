import astropy.io.ascii
import re
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Any

import pandas as pd
from astropy import units
from astropy.coordinates import Angle
from astropy.table import Table

from moa_exoplanet_archive.metadata_column_name import AlertMetadata, CandidateMetadata, TargetMetadata, \
    MetadataColumnName, TakahiroSumiCandlistRaDecFileColumnName, TakahiroSumiCandlistFileColumnName, \
    TakahiroSumiCandlistAlertIdFileColumnName, metadata_units_dictionary, \
    metadata_types_dictionary
from moa_exoplanet_archive.paths_and_names import local_cut0_object_directory

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
            x__pixels, y__pixels, ra, dec, tag = self.get_basic_target_data(field, chip, subframe, id_)
            candidate_metadata = self.get_candidate_metadata_for_light_curve_from_candlist(light_curve_path)
            alert_metadata_list = self.get_alert_metadata_list_for_light_curve(light_curve_path)
            target_metadata = TargetMetadata(field=field, chip=chip, subframe=subframe, id=id_, tag=tag,
                                             x__pixels=x__pixels, y__pixels=y__pixels, ra_j2000=ra, dec_j2000=dec,
                                             candidate_metadata=candidate_metadata,
                                             alert_metadata_list=alert_metadata_list)
            target_exoplanet_archive_dictionary = self.target_metadata_to_exoplanet_archive_dictionary(target_metadata)
            target_exoplanet_archive_dictionary_list.append(target_exoplanet_archive_dictionary)
        exoplanet_archive_table = Table(data=target_exoplanet_archive_dictionary_list,
                                        names=[name for name in MetadataColumnName],
                                        dtype=metadata_types_dictionary.values(),
                                        units=metadata_units_dictionary,
                                        )
        output_metadata_path = Path('metadata.ipac')
        astropy.io.ascii.write(exoplanet_archive_table, output=output_metadata_path, format='ipac', overwrite=True)

    def get_basic_target_data(self, field: int, chip: int, subframe: int, id_: int
                              ) -> (float, float, float, float, Optional[str]):
        x__pixels, y__pixels, ra, dec, tag = self.get_basic_target_data_from_object_file(
            field, chip, subframe, id_)
        result = self.get_basic_target_data_from_candlist(field, chip, subframe, id_)
        if result is not None:
            cl_x__pixels, cl_y__pixels, cl_ra, cl_dec, cl_tag = result
            assert cl_x__pixels == x__pixels
            assert cl_y__pixels == y__pixels
            assert cl_ra == ra
            assert cl_dec == dec
            tag = cl_tag
        return x__pixels, y__pixels, ra, dec, tag

    def get_basic_target_data_from_object_file(self, field, chip, subframe, id_):
        object_path = Path(local_cut0_object_directory).joinpath(f'cut0-gb{field}-{chip}-{subframe}.obj')
        object_data_frame = pd.read_csv(object_path, delim_whitespace=True, skipinitialspace=True, comment='#',
                                        names=['id', 'x', 'y'], usecols=[0, 1, 2])
        light_curve_row = object_data_frame[(self.candlist_data_frame['id'] == id_)].iloc[0]
        x__pixels = light_curve_row['x']
        y__pixels = light_curve_row['y']
        ra, dec = self.get_ra_and_dec_for_ccd_position(field=field, chip=chip, x__pixels=x__pixels, y__pixels=y__pixels)
        tag = None
        return x__pixels, y__pixels, ra, dec, tag

    def get_basic_target_data_from_candlist(self, field, chip, subframe, id_):
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
        return x__pixels, y__pixels, ra, dec, tag

    def extract_field_chip_subframe_and_id_from_light_curve_path(self, light_curve_path) -> (int, int, int, int):
        match = re.match(r'gb(\d+)-R-(\d+)-(\d+)-(\d+).ipac.gz', light_curve_path.name)
        if match is None:
            raise ValueError(f'Could not find light curve naming in {light_curve_path}')
        field = int(match.group(1))
        chip = int(match.group(2))
        subframe = int(match.group(3))
        id_ = int(match.group(4))
        return field, chip, subframe, id_

    def get_candidate_metadata_for_light_curve_from_candlist(self, light_curve_path: Path
                                                             ) -> Optional[CandidateMetadata]:
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
            MetadataColumnName.FIELD: target_metadata.field,
            MetadataColumnName.CHIP: target_metadata.chip,
            MetadataColumnName.SUBFRAME: target_metadata.subframe,
            MetadataColumnName.ID: target_metadata.id,
            MetadataColumnName.TAG: target_metadata.tag,
            MetadataColumnName.X: target_metadata.x__pixels,
            MetadataColumnName.Y: target_metadata.y__pixels,
            MetadataColumnName.RA_J2000: target_metadata.ra_j2000.value,
            MetadataColumnName.DEC_J2000: target_metadata.dec_j2000.value,
        })
        candidate_metadata = target_metadata.candidate_metadata
        if candidate_metadata is not None:
            exoplanet_archive_dictionary.update({
                MetadataColumnName.NUMBER_OF_DATA_POINTS: candidate_metadata.number_of_data_points,
                MetadataColumnName.NUMBER_OF_FRAMES_OBJECT_IS_DETECTED:
                    candidate_metadata.number_of_frames_object_is_detected,
                MetadataColumnName.MAX_SIGNIFICANCE: candidate_metadata.max_significance,
                MetadataColumnName.SUM_OF_CONTINUOUS_SIGNIFICANCE:
                    candidate_metadata.sum_of_continuous_significance,
                MetadataColumnName.CHI_SQUARED_OUTSIDE_SEARCH_WINDOW:
                    candidate_metadata.chi_squared_outside_search_window,
                MetadataColumnName.DOPHOT_OBJECT_REFERENCE_IMAGE_SEPARATION:
                    candidate_metadata.dophot_object_reference_image_separation__pixels,
                MetadataColumnName.DOPHOT_ID: candidate_metadata.dophot_id,
                MetadataColumnName.DOPHOT_TYPE: candidate_metadata.dophot_type,
                MetadataColumnName.DOPHOT_MAGNITUDE: candidate_metadata.dophot_magnitude,
                MetadataColumnName.DOPHOT_MAGNITUDE_ERROR: candidate_metadata.dophot_magnitude_error,
                MetadataColumnName.PSPL_T0: candidate_metadata.pspl_t0,
                MetadataColumnName.PSPL_TE: candidate_metadata.pspl_tE,
                MetadataColumnName.PSPL_U0: candidate_metadata.pspl_u0,
                MetadataColumnName.PSPL_SOURCE_FLUX: candidate_metadata.pspl_source_flux,
                MetadataColumnName.PSPL_BLENDING_FLUX: candidate_metadata.pspl_blending_flux,
                MetadataColumnName.PSPL_T0_PARABOLIC_ERROR: candidate_metadata.pspl_t0_parabolic_error,
                MetadataColumnName.PSPL_TE_PARABOLIC_ERROR: candidate_metadata.pspl_tE_parabolic_error,
                MetadataColumnName.PSPL_TE_ERROR_LOWER_LIMIT: candidate_metadata.pspl_tE_error_lower_limit,
                MetadataColumnName.PSPL_TE_ERROR_UPPER_LIMIT: candidate_metadata.pspl_tE_error_upper_limit,
                MetadataColumnName.PSPL_U0_PARABOLIC_ERROR: candidate_metadata.pspl_u0_parabolic_error,
                MetadataColumnName.PSPL_U0_ERROR_LOWER_LIMIT: candidate_metadata.pspl_u0_error_lower_limit,
                MetadataColumnName.PSPL_U0_ERROR_UPPER_LIMIT: candidate_metadata.pspl_u0_error_upper_limit,
                MetadataColumnName.PSPL_SOURCE_FLUX_ERROR: candidate_metadata.pspl_source_flux_error,
                MetadataColumnName.PSPL_BLENDING_FLUX_ERROR: candidate_metadata.pspl_blending_flux_error,
                MetadataColumnName.PSPL_CHI_SQUARED: candidate_metadata.pspl_chi_squared,
                MetadataColumnName.FSPL_T0: candidate_metadata.fspl_t0,
                MetadataColumnName.FSPL_TE: candidate_metadata.fspl_tE,
                MetadataColumnName.FSPL_U0: candidate_metadata.fspl_u0,
                MetadataColumnName.FSPL_RHO: candidate_metadata.fspl_rho,
                MetadataColumnName.FSPL_SOURCE_FLUX: candidate_metadata.fspl_source_flux,
                MetadataColumnName.FSPL_BLENDING_FLUX: candidate_metadata.fspl_blending_flux,
                MetadataColumnName.FSPL_T0_PARABOLIC_ERROR: candidate_metadata.fspl_t0_parabolic_error,
                MetadataColumnName.FSPL_TE_PARABOLIC_ERROR: candidate_metadata.fspl_tE_parabolic_error,
                MetadataColumnName.FSPL_TE_ERROR_LOWER_LIMIT: candidate_metadata.fspl_tE_error_lower_limit,
                MetadataColumnName.FSPL_TE_ERROR_UPPER_LIMIT: candidate_metadata.fspl_tE_error_upper_limit,
                MetadataColumnName.FSPL_U0_PARABOLIC_ERROR: candidate_metadata.fspl_u0_parabolic_error,
                MetadataColumnName.FSPL_U0_ERROR_LOWER_LIMIT: candidate_metadata.fspl_u0_error_lower_limit,
                MetadataColumnName.FSPL_U0_ERROR_UPPER_LIMIT: candidate_metadata.fspl_u0_error_upper_limit,
                MetadataColumnName.FSPL_RHO_PARABOLIC_ERROR: candidate_metadata.fspl_rho_parabolic_error,
                MetadataColumnName.FSPL_RHO_ERROR_LOWER_LIMIT: candidate_metadata.fspl_rho_error_lower_limit,
                MetadataColumnName.FSPL_RHO_ERROR_UPPER_LIMIT: candidate_metadata.fspl_rho_error_upper_limit,
                MetadataColumnName.FSPL_SOURCE_FLUX_ERROR: candidate_metadata.fspl_source_flux_error,
                MetadataColumnName.FSPL_BLENDING_FLUX_ERROR: candidate_metadata.fspl_blending_flux_error,
                MetadataColumnName.FSPL_CHI_SQUARED: candidate_metadata.fspl_chi_squared,
            })
        alert_metadata_list = target_metadata.alert_metadata_list
        if alert_metadata_list is not None:
            alert_metadata0 = alert_metadata_list[0]
            exoplanet_archive_dictionary.update({
                MetadataColumnName.SEPARATION_TO_ALERT_POSITION1:
                    alert_metadata0.separation_to_alert_position__pixels,
                MetadataColumnName.ALERT_ID1: alert_metadata0.alert_id,
                MetadataColumnName.ALERT_X1: alert_metadata0.alert_x__pixels,
                MetadataColumnName.ALERT_Y1: alert_metadata0.alert_y__pixels,
            })
            if len(alert_metadata_list) > 1:
                alert_metadata1 = alert_metadata_list[1]
                exoplanet_archive_dictionary.update({
                    MetadataColumnName.SEPARATION_TO_ALERT_POSITION2:
                        alert_metadata1.separation_to_alert_position__pixels,
                    MetadataColumnName.ALERT_ID2: alert_metadata1.alert_id,
                    MetadataColumnName.ALERT_X2: alert_metadata1.alert_x__pixels,
                    MetadataColumnName.ALERT_Y2: alert_metadata1.alert_y__pixels,
                })
        return exoplanet_archive_dictionary


if __name__ == '__main__':
    MetadataProcessor().process_metadata()
