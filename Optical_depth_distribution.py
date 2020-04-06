"""
Jack J. Window
3rd Year Physics Undergraduate
Blackett Laboratory, Imperial College London

07.04.2020

Data analysis script for Solar Radiation 3rd year experiment.
Data taken on 09.03.2020 during cloudy conditions to study the probability distribution
of cloud optical depths in London.
"""

# ----------------------------------- IMPORT MODULES -----------------------------------
import csv
import scipy as sp
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------

def readData(filename):
    """
    Function to read data from a pre-organised CSV of solar instrument measurements.
    
    CSV is assumed to be of the format ->

    NUM | LIGHTMETER [KLux] | SOLARIMETER [mV]  <- Headers here given as indication. 
    ----------------------------------             CSV must have NO HEADERS to process
     1  |        lm1        |      sol1             the data properly.
     .  |        .          |       .
     .  |        .          |       .
     .  |        .          |       .
     60 |        lmN        |       solN           <- For N data points
    ----------------------------------

    PARAMETERS
    -------------------------------------------------------------------------------------
    filename    :   str, the name of the csv from which to read data if stored in the 
                    same directory as the script. Otherwise this should be a full filepath.
    
    RETURNS
    -------------------------------------------------------------------------------------
    N           :   int, number of measurements in csv.
    raw_lm      :   list, contains all raw lightmeter measurements (no conversion factor)
    raw_sol     :   list, contails all raw solarimeter measurements (no conversion factor)
    """   

    # read raw measurements from CSV   
    raw_lm = []
    raw_sol = []
    with open(filename) as solar:
        rawData = csv.reader(solar, delimiter=',')
        for row in rawData:
            # Populate column arrays with data from CSV, using format in docstring
            raw_lm.append(row[1])
            raw_sol.append(row[2])

    N = len(raw_lm)

    # Perform conversions
    # Solarimeter: mV -> W/m^2 = x14.98
    # Lightmeter: Klux -> flux = x7.9 (0.0079 flux per lux)

    lm = [7.9*_lm for _lm in raw_lm]
    sol = [14.98*_sol for _sol in sol]

    return N, lm, sol

