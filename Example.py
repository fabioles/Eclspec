
import eclspec

#name of the planet
planet = 'HD 209458 b'

#type of atmosphere to be simulated (hotjupiter, earth, venus)
atmosphere_type = 'hotjupiter'

#transmission or emission spectroscopy
observation_type = 'transmission'

#File with paths to other required files
pathfile = 'Pathfile.txt'

system = eclspec.System(planet, pathfile, atmosphere_type,
        observation_type)

events = system.FindEvents()

obs = eclspec.Observation(system, order_setting = 'K/2/4')

#%%
#index of the event to be simulated
event_index = 0

save_directory = 'Example_Result'

# obs.Observe(event_index, save_directory)

#%%

        
import os
from astropy.io import fits
import numpy as np

def GetPhoenixSpectrum(self, spectrum_type = 'HiRes'):
    '''Options for spectrum type are HiRes and SpecInt'''
    
    self.FindClosestPhoenixParameters()
    self.GetPhoenixSpectrumLocation(spectrum_type)
    
    if self.phoenix_wave is None:
        path_phoenix_wave = os.path.join(self.paths['phoenixgrid'], 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
        
        with fits.open(path_phoenix_wave) as phxwave_hdulist:
            self.phoenix_wave = np.array(phxwave_hdulist[0].data)

    if spectrum_type == 'HiRes':
        with fits.open(self.filepath_phoenix_HiRes) as phxspec_hdulist:
            self.phoenix_HiRes = phxspec_hdulist[0].data
            
    elif spectrum_type == 'SpecInt':
        with fits.open(self.filepath_phoenix_SpecInt) as phxspec_hdulist:
            self.phoenix_SpecInt = phxspec_hdulist[0].data
    
GetPhoenixSpectrum(system, 'HiRes')
    
#%%
with fits.open(system.filepath_phoenix_SpecInt) as phxspec_hdulist:
    phoenix_SpecInt = phxspec_hdulist[0].data
    header = phxspec_hdulist[0].header