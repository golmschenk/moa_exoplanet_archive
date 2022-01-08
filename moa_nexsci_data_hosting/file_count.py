from pathlib import Path

merged_version_directory = Path('merged_MOA2dia')
paths = merged_version_directory.glob('**/*.txt.gz')
print(f'{len(list(paths))}')
