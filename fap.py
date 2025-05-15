import numpy as np
import random
import pandas as pd
import lightkurve as lk

import warnings
warnings.filterwarnings('ignore')

def pyriod_to_file(pyriod, filename):
    """Saves the time and flux columns of a pyriod object's residual lightcurve (thus allowing for prewhitening) to a .dat file.

    pyriod is the name of a pyriod object
    filename is what the saved .dat will be called"""
    
    np.savetxt('{}.dat'.format(filename), np.vstack((pyriod.lc.time.value,np.array(pyriod.lc.resid.value))).T) # take the residual time and flux columns from a pyriod lightcurve and save them in a .dat

def file_to_lc(filename):
    """Determines the false alarm probability of a given periodogram when searching for pulsation periods.

    filename is what the .dat file containing time and flux is called"""

    data = pd.read_csv(str(filename), sep="\s+", comment = '#', header = None, names = ['time', 'flux', 'flux_err'])

    data['time'] = data['time']/86400

    data['flux'] += 1 #lightkurve likes relative flux

    lc = lk.LightCurve(time = data['time'], flux = data['flux'], flux_err = data['flux_err'])

    return(lc)

def false_alarm_prob_newtest(lc, cycles):
    """Determines the false alarm probability of a given periodogram when searching for pulsation periods.
    
    lc is a lightcurve object with time and flux columns
    cycles is the number of cycles you want the program to run for - the tenth highest peak after 10000 cycles corresponds to a 1% chance of something above the cutoff being a false amplitude"""

    lc_time = lc['time']

    lc_flux = lc['flux']
    lc_flux_shifted = lc_flux + 1

    max_amp_array = np.zeros(cycles) # an array to store each calcualted fap in
    
    time = np.zeros(len(lc_time)) # taking the time column from your residuals array
    flux = np.zeros(len(lc_flux_shifted)) # taking the flux column from your residuals array

    for j in range(len(lc_time)): # cleaning up the arrays
        time[j] = str(lc_time[j]) # making time into a string
        element = str(lc_flux_shifted[j]) # making residual flux into a string
        tempflux = element.split(' ') # getting rid of any units tacked onto residual flux

        # this piece may not be needed anymore but i haven't checked yet
        if tempflux[0] == '———': # deal with the --- that appears in TESS data sometimes
            tempflux[0] = np.nan # replace with nan

        flux[j] = tempflux[0] # write the residual flux value to our flux array after dealing with issues like units and ---

    for i in range(cycles): # this is where the monte carlo simulation of it all happens

        newflux = flux
        
        random.shuffle(newflux)
            
        shuffled_lc = lk.LightCurve(time = time, flux = newflux) # create a new lc with the randomized indices

        pg = shuffled_lc.to_periodogram(freq_unit = 'microHertz') # create a new periodogram from the new lc with the randomized indices
        
        max_amp_array[i] = pg.max_power*1000 # get our units correct for amplitude (make it into mma)

    if(cycles <= 10):
        return np.mean(max_amp_array)
    else:
        return -np.sort(-max_amp_array)[9]
        #return max_amp_array[9] # return the false amplitude probability cutoff, which is the mean of the peak amplitude for however many cycles you asked for
