    
# EclSpec
Simulation of eclipsed spectra

EclSpec is used to simulate transmission and emission spectroscopy observations of CRIRES+. To get started, the target planet, the type of its atmosphere and the type of observation (transmission or emission) have to be declared. The *pathfile* is a text file containing paths to various other files such as the atmosphere models, stellar spectra models and a list of observation opportunities computed with ObsOp. 

### Disclaimer: The wavelength solution and blaze function used in the current version are outdated! They will hopefully be replaced in the future with more accurate ones.

```
import eclspec

#name of the planet
planet = 'HD 209458 b'

#type of atmosphere to be simulated (hotjupiter, earth, venus)
atmosphere_type = 'hotjupiter'

#transmission or emission spectroscopy
observation_type = 'transmission'

#File with paths to other required files
pathfile = 'Eclspec_par.txt'

system = eclspec.System(planet, pathfile, atmosphere_type,
        observation_type)
```

The *system* object contains all the information about the planet and its host star, such as the ephemerides and the planetary and stellar spectrum. This also includes a list of the transit or occultation events of this planet, which can be accessed with:

```
events = system.FindEvents()
```

Next, the observation is set up by specifying the desired spectral order setting and creating an *Observation* object:
```
obs = eclspec.Observation(planet, order_setting = 'K/2/4')
```

To simulate the observation of an event in the *events* list, specify the index of this event and a directory where the results are saved.

```
#index of the event to be simulated
event_index = 0

save_directory = 'HD209458b_Transmission_Simulation'

obs.Observe(event_index, save_directory)
```

The result consists of a collection of FITS-files containing the spectra, and an accompanying CSV-file with the exact time, airmass and barycentric velocity of each exposure.
