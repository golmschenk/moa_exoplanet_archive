import gzip
import multiprocessing
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from tabulate import tabulate

from moa_nexsci_data_hosting.column_name import phot_all_column_names, phot_cor_column_names, ColumnName, \
    merged_column_names, merged_column_formats


def merge_split_version_files_into_single_file(dot_phot_dot_all_path: Path, destination_merged_path: Path):
    containing_directory_path = dot_phot_dot_all_path.parent
    assert dot_phot_dot_all_path.name.endswith('.phot.all.gz')
    assert dot_phot_dot_all_path.exists()
    dot_phot_name = dot_phot_dot_all_path.name.replace('.phot.all.gz', '.phot.gz')
    dot_phot_path = containing_directory_path.joinpath(dot_phot_name)
    assert dot_phot_path.name.endswith('.phot.gz')
    assert dot_phot_path.exists()
    dot_phot_dot_cor_name = dot_phot_dot_all_path.name.replace('.phot.all.gz', '.phot.cor.gz')
    dot_phot_dot_cor_path = containing_directory_path.joinpath(dot_phot_dot_cor_name)
    assert dot_phot_dot_cor_path.name.endswith('.phot.cor.gz')
    assert dot_phot_dot_cor_path.exists()
    dot_phot_dot_all_data_frame = pd.read_csv(dot_phot_dot_all_path, comment='#', names=phot_all_column_names,
                                              delim_whitespace=True, skipinitialspace=True)
    dot_phot_dot_cor_data_frame = pd.read_csv(dot_phot_dot_cor_path, comment='#', names=phot_cor_column_names,
                                              delim_whitespace=True, skipinitialspace=True)
    columns_to_merge_in = dot_phot_dot_cor_data_frame.columns.difference(dot_phot_dot_all_data_frame.columns).tolist()
    columns_to_merge_in.append(ColumnName.HJD)
    unordered_merged_data_frame = pd.merge(dot_phot_dot_all_data_frame,
                                           dot_phot_dot_cor_data_frame[columns_to_merge_in], on=ColumnName.HJD,
                                           how='outer')
    if unordered_merged_data_frame[ColumnName.AIRMASS_1].isna().all():
        unordered_merged_data_frame[ColumnName.COR_FLUX] = np.nan
    merged_data_frame = unordered_merged_data_frame[merged_column_names]
    table_string = tabulate(merged_data_frame, headers=merged_column_names, tablefmt='plain', missingval='na',
                            floatfmt=tuple(merged_column_formats.values()), showindex=False, numalign='right')
    table_string = table_string.replace('nan', ' na')

    with gzip.open(destination_merged_path, 'wt') as merged_file:
        merged_file.write(table_string)


def merge_three_version(phot_all_path, three_version_directory, merged_version_directory):
    merged_sub_path = phot_all_path.relative_to(three_version_directory)
    merged_path = merged_version_directory.joinpath(merged_sub_path)
    merged_path = merged_path.parent.joinpath(merged_path.name.replace('.phot.all.gz', '.txt.gz'))
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    merge_split_version_files_into_single_file(phot_all_path, merged_path)


def convert_directory_to_merged_version(three_version_directory: Path, merged_version_directory: Path):
    merge_three_version_with_roots = partial(merge_three_version,
                                             three_version_directory=three_version_directory,
                                             merged_version_directory=merged_version_directory)
    paths = three_version_directory.glob('**/*.phot.all.gz')
    with multiprocessing.get_context('spawn').Pool(processes=20, maxtasksperchild=200) as pool:
        for _ in pool.imap_unordered(merge_three_version_with_roots, paths):
            pass


if __name__ == '__main__':
    three_version_directory_ = Path('MOA2dia')
    merged_version_directory_ = Path('merged_MOA2dia_cor_flux_na')
    convert_directory_to_merged_version(three_version_directory_, merged_version_directory_)
