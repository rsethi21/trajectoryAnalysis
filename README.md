# trajectoryAnalysis
Original trajectory data are available in 
 /home/pkekeneshuskey/data/molecular_dynamics/ph_binsun
 (from Bin Sun)

## Information about trajectory
- Currently simulates mitochondrial IDP Q...
- Has 100 copies and between 6-800 frames of protein movement simulated at pH 3, 7
- pH 3 is assumed to be unfolded and pH 7 is assumed to be folded
- ""_all.pdb includes ions/water
- ""_protein.pdb includes proteins only

# scripts under tools folder
## processSeries.py
- Currently will output time series Radgyr, RMSF, WaterShell, IonShell from pytraj utilities

## processContacts.py
- Currently will output contacts of protein backbone and sidechains over time

## clustering ipynb
- Currently clusters center of masses for amino acid residues and will be used to collect statistics on the clustering over time

# side notes
- needed to install pytraj from source since its setup.py was not functioning in any version of python3
