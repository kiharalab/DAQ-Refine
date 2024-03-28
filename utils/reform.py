import logging
import subprocess
from pathlib import Path
import os
def reform(map_path: Path, reformed_path: Path):

    if reformed_path.exists():
        logging.info(f"Reformed map already exists.")
        pass
    else:
        reformed_path.parent.mkdir(parents=True, exist_ok=True)
    command = ['/bio/dragon/webdklab/.local/bin/julia', '/bio/kihara-web/www/em/emweb-jobscheduler/reform.jl', f'--map_path={map_path}', f'--reformed_path={reformed_path}']
    process = subprocess.run(command, capture_output=True, text=True)

    # Check if the command was successful
    if process.returncode != 0:
        logging.error(f'Julia Subprocess Error: {process.stderr}')
    else:
        logging.info(process.stdout)
import sys
reform(Path(sys.argv[1]),Path(sys.argv[2]))
