import gzip
from multiprocessing import Pool
from pathlib import Path

import tqdm as tqdm


def post_process():
    light_curve_root_directory = Path('/local/data/emu/share/exoplanet_archive_moa_dataset/exoplanet_archive_moa_dataset/lcurve/lcurve')
    light_curve_path_list = list(light_curve_root_directory.glob(r'**/*.ipac.gz'))
    with Pool(20) as pool:
        for _ in tqdm.tqdm(pool.imap_unordered(replace_null_string, light_curve_path_list), total=len(light_curve_path_list)):
            pass


def replace_null_string(light_curve_path):
    with gzip.open(light_curve_path, 'rt') as unprocessed_file:
        unprocessed_light_curve_string = unprocessed_file.read()
    processed_light_curve_string = unprocessed_light_curve_string.replace('null|', ' nan|')
    with gzip.open(light_curve_path, 'wt') as processed_file:
        processed_file.write(processed_light_curve_string)


if __name__ == '__main__':
    post_process()
