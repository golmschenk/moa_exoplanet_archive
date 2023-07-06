import gzip
from pathlib import Path


def post_process():
    light_curve_root_directory = Path('/local/data/emu/share/exoplanet_archive_moa_dataset/exoplanet_archive_moa_dataset/lcurve/lcurve')
    light_curve_path_glob = light_curve_root_directory.glob(r'**/*.ipac.gz')
    for light_curve_index, light_curve_path in enumerate(light_curve_path_glob):
        print(f'{light_curve_index}: {light_curve_path}')
        with gzip.open(light_curve_path, 'rt') as unprocessed_file:
            unprocessed_light_curve_string = unprocessed_file.read()

        processed_light_curve_string = unprocessed_light_curve_string.replace('null|', ' nan|')

        with gzip.open(light_curve_path, 'wt') as processed_file:
            processed_file.write(processed_light_curve_string)


if __name__ == '__main__':
    post_process()
