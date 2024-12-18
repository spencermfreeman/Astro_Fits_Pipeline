from calibration import create_master_frames, calibrate_light_frames
import matplotlib.pyplot as plt
from photometry import perform_photometry
import numpy as np

def calibrate(dir, transit_name):
    master_flat, master_bias = create_master_frames(dir)
    lights_calibrated = calibrate_light_frames(dir, transit_name, master_flat, master_bias)
    #in case we want to see master composite calibration frames
    return lights_calibrated

#this will become a measurement of realtive flux
def plot(time, target_flux, comp_flux, title, dir):
    fig, axs = plt.subplots()
    target_flux = np.array(target_flux)
    comp_flux = np.array(comp_flux)
    rel_flux = np.array(comp_flux/target_flux)
    axs.scatter(time, rel_flux)
    axs.set_title(title)
    plt.savefig(dir)
    plt.show()
    
if __name__ == "__main__":
    dir = '/Users/spencerfreeman/Desktop/stepUp/2024-09-4-skuban'
    transit_name = 'qatar-5b'
    lights_calibrated = calibrate(dir, transit_name)
    #lights_calibrated is a list of CCDData objects, arrays are accsessd by data attribute.
    print(lights_calibrated)     
    target_flux, comparison_flux, time = perform_photometry(lights_calibrated)
    plot(time, target_flux, comparison_flux, transit_name, dir)
    