# Receiver-Deployment
Supplementary software for "Frontiers" journal paper "Optimal Deployment of Anchors for Underwater Acoustic Localization" 

This repository is supplementary software for "Frontiers" journal paper "Optimal Deployment of Anchors for Underwater Acoustic Localization".
The three programs are show_areas, find_areas and utils. Also, some reproduction data is given. 
show_areas:
For a given receiver deployment show the usable and coverage areas.
If 'bellhop=1' use the TL map from Bellhop. Else assume cylindrical spread
program adjusted for SD i.e.
grid size - according to the area's grid
TL - two options for the TL map
    a. TL of San-Diego for typical SVP
    b. TL of San-Diego for iso SVP
TL map was calculated using Bellhop package
Coverage files -
either set manually the receivers' positions by variable POP or read results of
'find_areas.py'. Here the results for 5 and 10 receivers for typical SVP and
for 10 receivers for iso SVP are provided.
TL map and Coverage results files are the
    'SD_tl_data_50x40_iso' and
    'cover_data_50x40.plk'
Respectively

find_areas
Find optimal receiver deployment search for N receivers using evolutionary
algorithms base on DEAP and SCOOP packages
Adjusted to San-Diego thus, uses 50X40 grid with BELLHOP TL calculation
use the 'bellhop' flag to decide if TL map is used or cylinder spread in use
Take care to adjust:
    - TL file name
    - grid size
    - area of interest size
    - bound for evolutionary algo
    - number of receivers N
    - uncomment the future.map line of wishing to run on parallel cores
    - input ga parameters e.g., number of generations

