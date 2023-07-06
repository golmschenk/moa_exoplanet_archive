import re
from pathlib import Path

import pysftp as pysftp

from moa_exoplanet_archive.paths_and_names import osaka_entry_in_ssh_config, osaka_data_directories, \
    osaka_cut0_object_subdirectory
from secret import osaka_hostname, osaka_username, osaka_ssh_key_path


def transfer_all_cut0_object_files():
    def pass_function(_):
        pass

    file_list = []
    def transfer_cut0_object_file(path: Path):
        match = re.search(str(path), r'cut0-gb\d+-\d+\.obj')
        print(path)
        nonlocal file_list
        file_list.append(path)

    for data_path in osaka_data_directories:
        path_cut0_directory = data_path.joinpath(osaka_cut0_object_subdirectory)
        with pysftp.Connection(osaka_hostname, username=osaka_username, private_key=osaka_ssh_key_path) as sftp:
            sftp.walktree(remotepath=str(path_cut0_directory),
                          fcallback=transfer_cut0_object_file,
                          dcallback=pass_function, ucallback=pass_function)

    dups = set([x for x in file_list if file_list.count(x) > 1])
    print(dups)

if __name__ == '__main__':
    transfer_all_cut0_object_files()