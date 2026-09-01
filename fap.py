import random
import numpy as np
import pandas as pd
import lightkurve as lk

import time
from rich.progress import Progress

import warnings
warnings.filterwarnings('ignore')

def fap_local(pyriod, samples, days = True):
    """Determines the false alarm probability of a given periodogram when searching for pulsation periods.
    
    pyriod is a pyriod object
    
    samples is the number of samples you want the program to draw
    
    days is whether the time unit of the lightcurve is days. Default is True."""
    
    raw_data = np.vstack((pyriod.lc_resid.time.value, np.array(pyriod.lc_resid.flux.value))).T
    
    data = pd.DataFrame(raw_data, columns=['time', 'flux'])
    
    if days == False:
    
        data['time'] = data['time']/86400

    base_lc = lk.LightCurve(time = data['time'], flux = data['flux'])

    lc = base_lc.copy().remove_nans()

    lc.flux = lc.flux + 1

    max_amp_array = np.zeros(samples) # an array to store each calcualted fap in

    rng = np.random.default_rng()

    with Progress(auto_refresh = False) as p:
        
        t = p.add_task("Processing...", total=samples)

        for i in range(samples): # this is where the monte carlo simulation of it all happens
            lc.flux = rng.choice(lc.flux.value, len(lc.flux), replace = True)*lc.flux.unit

            pg = lc.to_periodogram(freq_unit = 'microHertz') # create a new periodogram from the new lc with the randomized indices
        
            max_amp_array[i] = pg.max_power*1000 # get our units correct for amplitude (make it into mma)
            
            if i%100 == 0:
                    
                p.update(t, advance=100, refresh = True)

    print("The " + str(10/samples*100) + "% FAP is " + str(np.flip(np.sort(max_amp_array))[9]))
        
def fap_helios(lc_file, samples, days = True):
    """Determines the false alarm probability of a given periodogram when searching for pulsation periods.
    
    lc_file is a .dat file that contains time and flux columns
    
    samples is the number of samples you want the program to draw
    
    days is if the time of the lightcurve is in days. Default is True."""
    
    data = pd.read_csv(str(lc_file), sep=r"\s+", comment = '#', header = None, names = ['time', 'flux'])

    if days == False:
    
        data['time'] = data['time']/86400

    base_lc = lk.LightCurve(time = data['time'], flux = data['flux'])

    lc = base_lc.copy().remove_nans()
    
    lc.flux = lc.flux + 1

    max_amp_array = np.zeros(samples) # an array to store each calcualted fap in

    rng = np.random.default_rng()

    with Progress(auto_refresh = False) as p:
        
        t = p.add_task("Processing...", total=samples)

        for i in range(samples): # this is where the monte carlo simulation of it all happens
            lc.flux = rng.choice(lc.flux.value, len(lc.flux), replace = True)*lc.flux.unit
            
            pg = lc.to_periodogram(freq_unit = 'microHertz') # create a new periodogram from the new lc with the randomized indices
        
            max_amp_array[i] = pg.max_power*1000 # get our units correct for amplitude (make it into mma)
            
            if i%100 == 0:
                    
                p.update(t, advance=100, refresh = True)

    print("The " + str(10/samples*100) + "% FAP is " + str(np.flip(np.sort(max_amp_array))[9]))
        
def pyriod_to_file(pyriod, filename):
    """Saves the time and flux columns of a pyriod object's residual lightcurve (thus allowing for prewhitening) to a .dat file.

    pyriod is the name of a pyriod object
    filename is what the saved .dat will be called"""
    
    np.savetxt('{}.dat'.format(filename), np.vstack((pyriod.lc_resid.time.value, np.array(pyriod.lc_resid.flux.value))).T) # take the residual time and flux columns from a pyriod lightcurve and save them in a .dat
