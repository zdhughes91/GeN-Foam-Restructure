#!/bin/bash
# Command to build an application for SerpentToFoamXS
# Must be executed in the current folder.
# Author: Thomas Guilbaud, EPFL/Transmutex SA
# Last update: 18/03/2022

# Create build directory
buildFolder=build
mkdir $buildFolder

# Copy the source file to the build directory
cp SerpentToFoamXSApp.py $buildFolder
cp -r SerpentToFoamXS $buildFolder
cd $buildFolder

# Install PyInstaller to generate the SerpentToFoamXS executable
python3 -m pip install pyinstaller

# Find PyInstaller
pyin=$(sudo find /home/ -name pyinstaller)
echo $pyin

# Build the Application SerpentToFoamXS
$pyin SerpentToFoamXSApp.py SerpentToFoamXS/*.py SerpentToFoamXS/*/*.py --onefile --distpath bin/ --clean --hidden-import='PIL._tkinter_finder'

# End
echo "End of building"
