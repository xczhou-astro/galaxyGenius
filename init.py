import os
import argparse
from pathlib import Path
from galaxyGenius.config import Configuration

def none(value):
    if value == 'None':
        return None
    return value

parser = argparse.ArgumentParser()

parser.add_argument('-w', '--workspace', type=str, default='workspace')
parser.add_argument('-s', '--surveys', type=none, default='None')

args = parser.parse_args()

workspace = args.workspace
surveys = args.surveys

os.makedirs(workspace, exist_ok=True)

os.chdir(f'{workspace}')

config_path = Path('.')
config_files = list(config_path.glob('config*.toml'))

if len(config_files) > 0:
    print('config files found, cleaning and create new ones.')

for file in config_files:
    os.remove(file)

configuration = Configuration(surveys=surveys)
configuration.init()

print(f'Configuration files are created in {workspace}. Please edit them!')
