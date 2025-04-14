import logging
import subprocess
from pathlib import Path
import os
import sys
import platform

def download_and_install_julia(julia_install_path: Path):
    JULIA_VERSION = "1.9.1"  # Specify the desired version of Julia
    JULIA_VERSION_SHORT = "1.9"  # The version used in the app name
    
    # Determine the appropriate download URL based on the operating system
    system = platform.system().lower()
    if system == "darwin":  # macOS
        JULIA_URL = f"https://julialang-s3.julialang.org/bin/mac/x64/{JULIA_VERSION_SHORT}/julia-{JULIA_VERSION}-mac64.dmg"
    elif system == "linux":
        JULIA_URL = f"https://julialang-s3.julialang.org/bin/linux/x64/{JULIA_VERSION[:3]}/julia-{JULIA_VERSION}-linux-x86_64.tar.gz"
    else:
        raise OSError(f"Unsupported operating system: {system}")

    # Ensure the installation directory exists
    julia_install_path.parent.mkdir(parents=True, exist_ok=True)

    if system == "darwin":
        # For macOS, we need to handle DMG installation
        dmg_path = julia_install_path.with_suffix('.dmg')
        command_download = f"curl -L -o {dmg_path} {JULIA_URL}"
        command_mount = f"hdiutil attach {dmg_path}"
        command_copy = f"cp -R /Volumes/Julia-{JULIA_VERSION}/Julia-{JULIA_VERSION_SHORT}.app {julia_install_path.parent}/"
        command_unmount = f"hdiutil detach /Volumes/Julia-{JULIA_VERSION}"
        command_cleanup = f"rm {dmg_path}"
        
        subprocess.run(command_download, shell=True, check=True)
        subprocess.run(command_mount, shell=True, check=True)
        subprocess.run(command_copy, shell=True, check=True)
        subprocess.run(command_unmount, shell=True, check=True)
        subprocess.run(command_cleanup, shell=True, check=True)
        
        # Use the short version in the path
        julia_bin = julia_install_path.parent / f"Julia-{JULIA_VERSION_SHORT}.app/Contents/Resources/julia/bin/julia"
    else:
        # For Linux, use the existing tar.gz approach
        julia_tar_path = julia_install_path.with_suffix('.tar.gz')
        command_download = f"wget -O {julia_tar_path} {JULIA_URL}"
        command_extract = f"tar -xzf {julia_tar_path} -C {julia_install_path.parent} && rm {julia_tar_path}"

        subprocess.run(command_download, shell=True, check=True)
        subprocess.run(command_extract, shell=True, check=True)

        # Find the Julia binary inside the extracted folder
        julia_bin = next(julia_install_path.parent.glob(f"julia-{JULIA_VERSION}/bin/julia"))
    
    return julia_bin

def reform(map_path: Path, reformed_path: Path):
    # Define local Julia installation path
    local_julia_path = Path(__file__).parent / "julia_local" / "bin" / "julia"

    # Check if Julia is installed locally; if not, download and install it
    if not local_julia_path.exists():
        logging.info("Downloading and installing Julia locally...")
        local_julia_path = download_and_install_julia(local_julia_path.parent)

    if reformed_path.exists():
        logging.info("Reformed map already exists.")
    else:
        reformed_path.parent.mkdir(parents=True, exist_ok=True)

    command = [str(local_julia_path), '/bio/kihara-web/www/em/emweb-jobscheduler/reform.jl', f'--map_path={map_path}', f'--reformed_path={reformed_path}']
    process = subprocess.run(command, capture_output=True, text=True)

    # Check if the command was successful
    if process.returncode != 0:
        logging.error(f'Julia Subprocess Error: {process.stderr}')
    else:
        logging.info(process.stdout)

if __name__ == "__main__":
    reform(Path(sys.argv[1]), Path(sys.argv[2]))
