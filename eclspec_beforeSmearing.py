'''
Eclipsed spectrum simulator for CRIRES+

by Fabio Lesjak and Ansgar Reiners

July 2021
'''

import os
import wget
import pandas as pd
import numpy as np
from scipy import integrate
from scipy import ndimage
from astropy.io import fits
from PyAstronomy import pyasl
import pathlib
    
#constants and conversions
au = 1.495978707*10**11
au2sunr = 0.00465047 
c = 299792.458
G = 6.67430*10**(-11)
M_sun = 1.98847*10**30
M_jup = 1.89813*10**27
day2sec = 86400


def GetCoefficients(filename, order_setting):
    """
    returns an array with coefficients from Crires wiki files
    
    3 x n_orders x 3 with the coefficients for the correct setting (3 column is None for detector minmax)
    
    filename: where the coefficient table is stored
    order_setting: the crires order in the format 'Y/2/4'
    
    """
    wlcoeffs = pd.read_csv(filename, delim_whitespace=True, header=None)
    setting = order_setting[0]+'_'+order_setting[2]+'_'+order_setting[4] #correct formatting of setting
    setting_ind = np.where(wlcoeffs[2] == setting)[0][0]

    #Get the index where the following setting begins
    following_setting_ind = []
    for i in range(len(wlcoeffs[2])):
        if type(wlcoeffs[2][i]) == str:
            if len(wlcoeffs[2][i]) == 5:
                if i > setting_ind:
                    following_setting_ind.append(i)
    if len(following_setting_ind) == 0:
        following_setting_ind = len(wlcoeffs[2])
    else:
        following_setting_ind = following_setting_ind[0]

    #Get inidzes for each detector (the coefficients start one row further down)
    det1_ind = np.where(wlcoeffs[0] == 'DET1')[0]
    det1_ind = det1_ind[np.where(det1_ind > setting_ind)[0][0]]
    det2_ind = np.where(wlcoeffs[0] == 'DET2')[0]
    det2_ind = det2_ind[np.where(det2_ind > setting_ind)[0][0]]
    det3_ind = np.where(wlcoeffs[0] == 'DET3')[0]
    det3_ind = det3_ind[np.where(det3_ind > setting_ind)[0][0]]

    det1_array = np.array(wlcoeffs[:][det1_ind+1:det2_ind]).astype(float)
    det2_array = np.array(wlcoeffs[:][det2_ind+1:det3_ind]).astype(float)
    det3_array = np.array(wlcoeffs[:][det3_ind+1:following_setting_ind]).astype(float)

    if not det1_array.shape[0] == det2_array.shape[0] == det3_array.shape[0]:
        max_size = np.max([det1_array.shape[0], det2_array.shape[0], det3_array.shape[0]])
        while det1_array.shape[0] < max_size:
            det1_array = np.append(det1_array, [[np.nan, np.nan, np.nan]], axis=0)
        while det2_array.shape[0] < max_size:
            det2_array = np.append(det2_array, [[np.nan, np.nan, np.nan]], axis=0)
        while det3_array.shape[0] < max_size:
            det3_array = np.append(det3_array, [[np.nan, np.nan, np.nan]], axis=0)       
    
    return np.stack([det1_array, det2_array, det3_array], axis=0)

def ComputeSNR(m, t, airm, method='Crires', mband='K'):
    '''
    Computes SNR of stellar spectrum from magnitude, exposure time (in seconds) and airmass
    
    Parameters
    ----------
    m : float
        Magnitude in band specified with the parameter "mband"
    t : float
        Exposure time [s]
    airm : float
        Airmass
    method : string, optional, {"Crires", "Carmenes"}
        Specifiy the instrument for/ method with which to compute the SNR
    mband : string, optional, {"K, J"}
        Specifiy in which band the magnitude "m" is measured
        
    Returns
    -------
    SNR : float
        Resulting SNR (only a rough estimate)
    '''
    ###First compute SNR assuming airmass = 1
    extcof = 0.05 #extinction coefficient, see Lombardi et al., 2018
    #This is from old Carmenes documentation, factor 1.1774 so that it agrees better with Ansgars result
    if method == 'Carmenes':
        if mband == 'J':
            SNR_noairmass = 1.1774*100/np.sqrt(40*10**((m-4.2)/2.5))*np.sqrt(t)
        else:
            print('Use Jmag for calculation of Carmenes SNR.')
            SNR_noairmass = np.nan
    elif method == 'Crires_old':
        if mband == 'K':
            SNR_noairmass = 449.4241*np.sqrt(10**(-m/2.5))*np.sqrt(t)- 6.3144
        else:
            print('Use Kmag for calculation of old Crires SNR.')
            SNR_noairmass = np.nan
    elif method == 'Crires':
        if mband == 'K':
            snr_airm1_2 = 247.31342303*np.sqrt(10**(-m/2.5))*np.sqrt(t) - 3.20150241
            SNR_noairmass = snr_airm1_2 * 10**(extcof/5*(1.2 - 1))
        elif mband == 'J':
            snr_airm1_2 = 479.05726751*np.sqrt(10**(-m/2.5))*np.sqrt(t) - 2.90701871
            SNR_noairmass = snr_airm1_2 * 10**(extcof/5*(1.2 - 1))
        else:
            print('Use Kmag for calculation of Crires SNR.')
            SNR_noairmass = np.nan
        
    else:
        print('Method not recognized. Use Crires or Carmenes.')
        SNR_noairmass = np.nan
              
    #Scale to airmass = airm
    SNR = SNR_noairmass * 10**(extcof/5*(1 - airm))
    
    return SNR



def IntersectArea(rst, rpl, d):
  #area of intersection between two circles
  #radii of circles rst and rpl
  #distance between circles d
    
    if len(locals()) < 3:
        print('usage: aintersect(rad1, rad2, d)')
    else:
        Aintersect = np.copy(d)  # initialize

        ind = np.where(d > (rst+rpl))
        nind = len(ind[0])

        if nind > 0:
            #no intersection
            Aintersect[ind] = 0.

        ind = np.where(d < (rst-rpl))
        nind = len(ind[0])
        if nind > 0:
            #one circle fully contained in other
            Aintersect[ind] = np.pi*rpl**2

        ind = np.where(((d >= (rst-rpl)) & (d <= (rst+rpl))))

        nind = len(ind[0])
        if nind > 0:
            #partial overlap
            d1 = (rst**2 - rpl**2 + d[ind]**2) / (2 * d[ind])
            d2 = d[ind] - d1
            Aintersect[ind] = rst**2*np.arccos(d1/rst) - d1*np.sqrt(rst**2 - d1*d1) + rpl**2*np.arccos(d2/rpl) - d2*np.sqrt(rpl**2 - d2*d2)
        
        return Aintersect
    
def Mag_spectrum(spectrum, mag):
    #enhances normalized absorption spectrum by factor of mag according to curve of growth
    #continuum is set to 1, has to be rewritten if this changes
    #original idl pro has cont as variable, thats why some steps seem unnecessary
    x0 = np.copy(spectrum)
    y0 = 1 - np.copy(spectrum)
    ind_mask = (y0 >= 1)

    x0[ind_mask] = 0
    x0[~ind_mask] = -np.log(x0[~ind_mask])


    x1 = x0 * mag
    y1 = np.exp(-x1)
    y1[ind_mask] = 0
    
    return y1

def GaussPsf(npix, fwhm, normalize=True): 
    # Initialize PSF params 
    cntrd = (npix - 1.0) * 0.5 #centre
    st_dev = 0.5 * fwhm / np.sqrt(2.0 * np.log(2)) 
    # Make PSF 
    i = range(npix) 
    psf = np.array( [np.exp(-(((cntrd-x)/st_dev)**2/2)) for x in i] ) 
    # Normalize 
    if normalize: psf /= psf.sum() 
    return psf


def Instrbrod(wave, spec, R):
    #instrumental broadening
    dw = np.median(wave[1:]-wave[:-1]) #median wavelenght difference
    npix = 2*int(2*wave[int(wave.size/2)]/R/dw)+1
    FWHM = wave[int(wave.size/2)]/R/dw
    sigma = FWHM / (2*np.sqrt(np.log(2)*2))
    
    lsf = GaussPsf(npix, FWHM)
    
    if len(spec.shape) == 1:
        brodspec = ndimage.convolve(spec, lsf, mode='wrap')
    elif len(spec.shape) == 2:
        brodspec = [ndimage.convolve(spec[i,:], lsf, mode='wrap') for i in range(spec.shape[0])]
    else:
        raise Exception('Shape of spectrum not compatible for instrumental broadening function')

    return brodspec
    
def ld(r, coeff):
    #limb darkening as function of radial distance r to center
    #coeff is a dataframe created by GetLimbDarkCoefficients
    ld_set = [coeff.a1, coeff.a2, coeff.a3, coeff.a4] 
    mu = np.sqrt(1 - r**2)
    ld = 1
    
    for k in range(1,5):
        ld -= ld_set[k-1]*(1 - mu**(k/2))
    return ld

    

def ApplyBlaze(det_wave_org, det_spec_org, wave_coeffs, blaze_coeffs, normalization):
    """
    for one detector, applies the according blaze function to each order
    
    det_wave: contains wavelengths for each order of the detector
    det_spec: contains spectrum for each order and each timestep
    wave_coeffs: contains all available coefficients for the wavelength polynomials for that detector
    blaze_coeffs: contains all available coefficients for the blaze function polynomials for that detector
    normalization: highest value of all blaze functions (across all 3 detectors) for the setting
    """
    det_wave = det_wave_org.copy()
    det_spec = det_spec_org.copy()
    pix = np.arange(2048)+1  #the fit to the blaze function is defined on this pixel array
    
    n = det_spec[0].shape[0]
  
    blaze_w = []  #wavelengths of the blaze function fit

    for coeff in wave_coeffs:
        #*10 to convert from nm to Angstrom
        blaze_w.append(np.polyval(coeff,pix)*10)
    
    #determine which blaze function corresponds to which detector order by comparing mean wavelenghts
    det_mean = []
    blaze_mean = []
    for wave in det_wave:
        det_mean.append(np.mean(wave))
    for wave in blaze_w:
        blaze_mean.append(np.mean(wave))

    #for each det order, the result is the corresponding index in the blaze coefficient array
    detind_to_blazeind = np.empty(len(det_mean))*np.nan

    for i in range(len(det_mean)):
        if not np.isnan(det_mean[i]):
            #use arbitrary value of 20 as minimal required distance between means to result in a match
            closest_ind = np.where(abs(blaze_mean - det_mean[i]) < 20)[0]

            if len(closest_ind) > 0:
                detind_to_blazeind[i] = closest_ind[0]

    #get blaze functions for all orders and normalize so that highest value is 1 across entire detector
    #also resample the spectrum on the grid of 2048 pixels for each order

    blaze_spec = np.ones([len(det_spec),2048]) * np.nan

    #a new array where the resampled wavelenghts for each order are saved
    det_spec_resampled = np.ones([len(det_spec),n,2048]) * np.nan
    for i in range(len(det_mean)):
        if np.isnan(detind_to_blazeind[i]):
            #set orders to nan if there is no corresponding blaze function
            det_wave[i] = pix * np.nan
        else:
            #add blaze function of this order to list
            blaze_spec[i] = np.polyval(blaze_coeffs[int(detind_to_blazeind[i])],pix)

            #resample
            for j in range(det_spec[i].shape[0]):
                det_spec_resampled[i][j] = np.interp(blaze_w[int(detind_to_blazeind[i])], det_wave[i], det_spec[i][j])
            
    blaze_spec /= normalization

    #apply blaze function to the spectrum for each order
    for i in range(len(blaze_spec)):
        det_spec_resampled[i] *= blaze_spec[i]
    
    #arange wavelengths in the correct format with nans
    det_wave_resampled = np.ones([len(det_spec),2048]) * np.nan
    for i in range(len(det_spec)):
        if not np.isnan(detind_to_blazeind[i]):
            det_wave_resampled[i] = blaze_w[int(detind_to_blazeind[i])]
        
    return det_wave_resampled, det_spec_resampled


class System:
    '''
    Class for holding information about planet and its host star
    '''
    def __init__(self, name, pathfile, atmosphere_type='hotjupiter', obs_type = 'transmission'):

        self.name = name
        self.transits = None
        self.atmosphere_type = atmosphere_type
        self.obs_type = obs_type
        
        self.phoenix_wave = None
        
        #Read file containing paths
        self.paths = {}
        for line in open(pathfile):
            if len(line) > 1:
                key = line.split(':')[0].replace(' ','')
                cval = ':'.join(line.split(':')[1:]).replace(' ','').replace('\n','')
                self.paths[key] = cval
       
        #Get planet and star parameters from nexa
        self.GetPlanetData(self.paths['nexa'])
        
        #Read stellar and planetary spectra
        self.GetPlanetSpectrum()
        self.GetStellarSpectrum()
        
        #Set up keplarian orbits
        self.SetupOrbits()
        
        
        
    def GetPlanetData(self, path):
        '''
        Get parameters of planet and star from Nexa catalog
        '''
        if os.path.exists(path) == False:
            print('This file does not exist: {}'.format(path))
        else:
            nexaData = pd.read_csv(path)
    
            pl_index = np.where(nexaData['pl_name'] == self.name)[0]
    
            if pl_index.size == 0:
                print(self.name + ' not found in Nexa database')
            else:
                pl_index = pl_index[0]

                self.Period = nexaData['pl_orbper'][pl_index]
                self.Inc = nexaData['pl_orbincl'][pl_index]
                self.Ecc = nexaData['pl_orbeccen'][pl_index]
                self.Periastron = nexaData['pl_orblper'][pl_index]
                ProjObliquity = nexaData['pl_projobliq'][pl_index]
                if pd.isnull(ProjObliquity):
                    self.ProjObliquity = 0
                else:
                    self.ProjObliquity = nexaData['pl_projobliq'][pl_index]
                self.Rs = nexaData['st_rad'][pl_index]
                self.Rp = nexaData['pl_radj'][pl_index] / self.Rs / 9.946  #planet radius [stellar radii]
                self.Mp = nexaData['pl_bmassj'][pl_index]
                self.Ms = nexaData['st_mass'][pl_index]
                self.Tdur = nexaData['pl_trandur'][pl_index] / 24 #convert from hours to days
                self.Tdepth = nexaData['pl_trandep'][pl_index]
                self.ImpactPar = nexaData['pl_imppar'][pl_index]
                self.Teq = nexaData['pl_eqt'][pl_index]
                self.Teff = nexaData['st_teff'][pl_index]
                self.Vsini = nexaData['st_vsin'][pl_index]
                self.RadVelStar = nexaData['st_radv'][pl_index]
                self.Hmag = nexaData['sy_hmag'][pl_index]
                self.Jmag = nexaData['sy_jmag'][pl_index]
                self.Kmag = nexaData['sy_kmag'][pl_index]
                self.Vmag = nexaData['sy_vmag'][pl_index]
                self.Ra = nexaData['ra'][pl_index]
                self.Dec = nexaData['dec'][pl_index]
                                
                self.logg = nexaData['st_logg'][pl_index]
                self.Z = nexaData['st_mass'][pl_index]
                self.H = self.Calc_ScaleHeight()
                
                sma_au = ((self.Period*day2sec)**2 * G * (self.Ms*M_sun + self.Mp*M_jup) / (4*np.pi**2))**(1/3) / au
                self.Sma = sma_au/au2sunr/self.Rs
                
    def Calc_ScaleHeight(self, mean_mol_weight=None):
        '''
        Calculate scale height of planetary atmosphere
        '''
        #scale height
        grav_constant = 6.67430*10**-11
        k_b = 1.380649*10**-23
        m_H = 1.6735575*10**-27 #mass of hydrogen atom [kg]
        surface_gravity = grav_constant * self.Mp * 1.898*10**27 / (self.Rp * self.Rs * 696340000)**2

        if pd.isnull(mean_mol_weight):
            if self.Rp*11.2095 < 1.5:
                #water atmosphere for small planets with R < 1.5 R_earth
                mean_mol_weight = 18
            else:
                #for hot jupiter atmosphere
                mean_mol_weight = 2.3 

        if np.isnan(self.Teq):
            #if no value for Teq, calculate it using formula from Kempton 2018
            Teq = self.Teff * np.sqrt(self.Rs*696340 / (self.Sma * 1.496e+8)) * 0.25**0.25
            
        H = k_b * self.Teq / (mean_mol_weight * m_H * surface_gravity)

        return H   
    
    def GetPlanetSpectrum(self):
        '''
        Load model spectrum for planet
        '''
        valid_atmosphere_types = ['hotjupiter', 'venus', 'earth', 'gj1214']
        valid_observation_types = ['transmission', 'emission']
        
        if not self.atmosphere_type in valid_atmosphere_types:
            raise Exception('The atmosphere type "{}" was not recognized.'.format(self.atmosphere_type))
        if not self.obs_type in valid_observation_types:
            raise Exception('The observation type "{}" was not recognized.'.format(self.obs_type))
        

        #load as pickle, and if that doesnt exist create new pickle file       
        if not self.atmosphere_type + '_' + self.obs_type in self.paths:
            raise Exception('No {} model availabel for {}-like atmosphere'.format(self.obs_type, self.atmosphere_type))
        else:
            atmosphere_path = self.paths[self.atmosphere_type + '_' + self.obs_type]
            
        pickle_path = atmosphere_path[:-3]+'pkl'
        try:
            spectrum_pl = pd.read_pickle(pickle_path)
        except:
            print('Planetary Atmosphere model: No pickled file found. Reading txt-file and creating new pickle.')
            spectrum_pl = pd.read_csv(atmosphere_path,
                      delim_whitespace=True, skiprows = 1, names = ['Wavelength', 'Spectrum'])
            spectrum_pl.to_pickle(pickle_path)

        if self.obs_type == 'transmission':
            self.wavelength_p = spectrum_pl.Wavelength * 10**4
            #For transmission, spectrum_pl.Spectrum is H in km
            H_lambda = spectrum_pl.Spectrum / 695700 / self.Rs

            #scale to scale height of planet
            reference_H = {'hotjupiter': 189064, #scale height of HD189733 in meter
                           'venus': 11000, #scale height of venus-like with logg=3, Teq=240K, in meter
                           'earth': 8500, #height of Earth's atmosphere
                           'gj1214': 247909 #scale height of GJ1214b
                          }
            reference_mol_weight = {'hotjupiter': 2.3,
                                    'venus': 18,
                                    'earth': 18,
                                    'gj1214': 2.3
                                   }
            scalefactor = (self.Calc_ScaleHeight(mean_mol_weight=reference_mol_weight[self.atmosphere_type])
                           / reference_H[self.atmosphere_type])

            H_lambda *= scalefactor

            self.spectrum_p = H_lambda + self.Rp
            
        elif self.obs_type == 'emission':
            self.wavelength_p = spectrum_pl.Wavelength * 10**8

            #scale from flux/Hz to flux/lambda
            c_cm = 2.998*10**10
            self.spectrum_p = spectrum_pl.Spectrum * c_cm / spectrum_pl.Wavelength**2
            
    def GetStellarSpectrum(self):
        '''
        Load stellar spectrum from Phoenix synthetic spectral library
        '''
        ##Find the best fitting phoenix spectrum

        #Get parameters of available phoenix files
        phoenix_T = []
        phoenix_logg = []
        phoenix_Z = []
        phoenix_files = []
        files_in_phoenixfolder = os.listdir(self.paths['phoenixgrid'])
        for file in files_in_phoenixfolder:
            if file[-38:] == 'PHOENIX-ACES-AGSS-COND-2011-HiRes.fits':
                phoenix_T.append(float(file[3:8]))
                phoenix_logg.append(float(file[9:13]))
                phoenix_Z.append(float(file[13:17]))
                phoenix_files.append(file)

        phoenix_T = np.array(phoenix_T)
        phoenix_logg = np.array(phoenix_logg)
        phoenix_Z = np.array(phoenix_Z)

        close_T = abs(phoenix_T - self.Teff) == np.min(abs(phoenix_T - self.Teff))
        close_logg = abs(phoenix_logg[close_T] - self.logg) == np.min(abs(phoenix_logg[close_T] - self.logg))
        close_Z = abs(phoenix_Z[close_T][close_logg] - self.Z) == np.min(abs(phoenix_Z[close_T][close_logg] - self.Z))

        chosen_index = np.arange(len(phoenix_files))[close_T][close_logg][close_Z][0]
        path_phxfile = phoenix_files[chosen_index]

        path_phxwavefile = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'

        #read stellar spectrum
        with fits.open(os.path.join(self.paths['phoenixgrid'], path_phxwavefile)) as phxwave_hdulist:
            phxwave = np.array(phxwave_hdulist[0].data)

        with fits.open(os.path.join(self.paths['phoenixgrid'], path_phxfile)) as phxspec_hdulist:
            phxspec = phxspec_hdulist[0].data

        #resample phx spectrum so it is evenly spaced
        stepsize = 0.02
        self.wavelength_s = np.arange(phxwave.min(), phxwave.max(), stepsize)
        self.spectrum_s = np.interp(self.wavelength_s, phxwave, phxspec)

    
    
    def GetPhoenixSpectrum(self, spectrum_type = 'HiRes'):
        '''Options for spectrum type are HiRes and SpecInt'''
        
        self.FindClosestPhoenixParameters()
        self.GetPhoenixSpectrumLocation(spectrum_type)
    
    def FindClosestPhoenixParameters(self):
        
        Teff_Phoenix_Grid = np.append(np.arange(2300, 7000, 100), np.arange(7000, 12001, 200))
        logg_Phoenix_Grid = np.arange(0, 6.1, 0.5)
        Z_Phoenix_Grid = np.append(np.arange(-4, -2, 1), np.arange(-2, 1.1, 0.5))
        
        self.Teff_Phoenix = Teff_Phoenix_Grid[np.argmin(abs(self.Teff - Teff_Phoenix_Grid))]
        self.logg_Phoenix = logg_Phoenix_Grid[np.argmin(abs(self.logg - logg_Phoenix_Grid))]
        self.Z_Phoenix = Z_Phoenix_Grid[np.argmin(abs(self.Z - Z_Phoenix_Grid))]
        
        
    def GetPhoenixFilename(self, spectrum_type = 'HiRes'):
        if spectrum_type == 'HiRes':
            filename_template = 'HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z{z:+3.1f}/lte{teff:05d}-{logg:1.2f}{z:+3.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
            
        elif spectrum_type == 'SpecInt':
            filename_template = 'SpecIntFITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/Z{z:+3.1f}/lte{teff:05d}-{logg:1.2f}{z:+3.1f}.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
            
        filename = filename_template.format(teff = self.Teff_Phoenix, logg = self.logg_Phoenix, z = self.Z_Phoenix)
        
        return filename
            
    def DownloadPhoenixSpectrum(self, url, savepath):
        
            os.makedirs(os.path.dirname(savepath), exist_ok=True)
            try:
                wget.download(url, savepath)
            except:
                print('File not found on phoenix server: {}'.format(url))
                
                
    def GetPhoenixSpectrumLocation(self, spectrum_type = 'HiRes'):
        
        urlbase = 'ftp://phoenix.astro.physik.uni-goettingen.de/'
        filename = self.GetPhoenixFilename(spectrum_type = spectrum_type)
        
        url = urlbase + filename
        filepath = os.path.join(self.paths['phoenixgrid'], filename)
        
        if spectrum_type == 'HiRes':
            self.filepath_phoenix_HiRes = filepath
        elif spectrum_type == 'SpecInt':
            self.filepath_phoenix_SpecInt = filepath
        
        if not os.path.exists(filepath):
            self.DownloadPhoenixSpectrum(url, filepath)
            
    def SetupOrbits(self):
        '''
        Set up keplerian orbit objects (from PyAstronomy) for star and planet
        '''
        #planet orbit
        semimajor_a = self.Sma/au2sunr/self.Rs   #semimajor axis (stellar radius)

        self.Orbit_p = pyasl.KeplerEllipse(self.Sma, self.Period, e=0,
                                           Omega=self.ProjObliquity, i=self.Inc, w=0)

        #stellar orbit
        semimajor_a_s = self.Sma / (1 + self.Ms*M_sun/(self.Mp*M_jup))   #distance between star and barycenter (AU)
        self.Orbit_s = pyasl.KeplerEllipse(semimajor_a_s, self.Period, e=0,
                                           Omega=self.ProjObliquity, i=self.Inc, w=0)
        
    def FindEvents(self):
        '''
        Get available observation opportunities from list (computed using Obsop)
        '''
        if self.obs_type == 'transmission':
            list_path = self.paths['transitlist']
        elif self.obs_type == 'emission':  
            list_path = self.paths['eclipselist']

        try:
            eventlist = pd.read_csv(list_path)

        except FileNotFoundError:
            print('Error: File was not found')
            raise

        else:   
            #Find events that lie within the timespan
            planet_events = eventlist[eventlist.System == self.name]
            if 'Unnamed: 0' in planet_events.columns:
                planet_events = planet_events.drop('Unnamed: 0', axis=1)

            good_cond = planet_events.GoodCond
            available_opp = planet_events[good_cond]

            available_opp = available_opp.reset_index(drop=True)

            if len(available_opp) == 0:
                print('No events found in the requested timespan!')
                    
            self.events = available_opp
            
            return available_opp
        
    

class Observation:
    '''
    Class for observation
    '''
    def __init__(self, planet, observatory = 'esoparanal', order_setting = 'K/2/4'):
        
        self.planet = planet
        self.paths = self.planet.paths
        
        self.readout_time = 14 / day2sec
        self.exposure_time = 300 / day2sec
        self.timestep = self.readout_time + self.exposure_time
        
        self.resolution = 100000  #instrument resolution
        self.sampling = 3.0 #instrument pixel sampling (pixel per resolution element)
        self.order_setting = order_setting #CRIRES order setting
        
        self.extcof = 0.034
        
        self.SetObservatory(observatory)
        
        self.opportunities = self.planet.FindEvents()
        
        #Load coefficients for detector wavelengths and blaze functions, compute blaze normalization
        self.GetDetectorCoefficients()
        
        self.wave_range = [np.nanmin(self.crformat)*(1-200/300000), np.nanmax(self.crformat)*(1+200/300000)]
        
        #Select regions in spectra that are in wave_range
        in_range = np.where((self.planet.wavelength_s < self.wave_range[1]) &
                            (self.planet.wavelength_s > self.wave_range[0]))[0]
        self.wavelength_s = self.planet.wavelength_s[in_range]
        self.spectrum_s = self.planet.spectrum_s[in_range]
        
        in_range = np.where((self.planet.wavelength_p < self.wave_range[1]) &
                            (self.planet.wavelength_p > self.wave_range[0]))[0]
        self.wavelength_p = self.planet.wavelength_p[in_range]
        self.spectrum_p = self.planet.spectrum_p[in_range]
        
        #Load telluric spectrum
        self.GetTelluricSpectrum()
        
        #Read limb darkening coefficients
        self.GetLimbDarkCoefficients(self.planet.Teff, self.planet.logg, self.planet.Z, self.order_setting[0])
        
        #compute total intensity
        nint = 1000
        rgrid = np.arange(nint)/(nint-1)
        self.total_int = integrate.simps(2*np.pi*rgrid*ld(rgrid, self.ld_coeff), rgrid)
        
    def SetObservatory(self, observatory):
        '''
        Set location of observatory using string conforming to pyasl names
        '''
        self.observatory = observatory
        observatory_data = pyasl.observatory(observatory)
        self.lon = observatory_data["longitude"]
        self.lat = observatory_data["latitude"]
        self.alt = observatory_data["altitude"]
        
    def GetTelluricSpectrum(self):
        with fits.open(self.paths['telluric']) as telldat_hdulist:
            telldat = telldat_hdulist[1].data
        in_range = np.where((telldat.wave > self.wave_range[0]) & (telldat.wave < self.wave_range[1]))
        self.wavelength_tell = telldat[in_range].wave
        self.spectrum_tell = telldat[in_range].flux
        
    def GetDetectorCoefficients(self):
        '''
        Reads the coefficients of the wavelength solution and blaze function
        '''
        #get table of start and end of each order for each detector
        self.crformat = GetCoefficients(self.paths['detector_wl'], self.order_setting) * 10

        #get blaze function coefficients
        self.blaze_wave = GetCoefficients(self.paths['blaze_wl'], self.order_setting)
        self.blaze_coeff = GetCoefficients(self.paths['blaze_coeff'], self.order_setting)
        
        #compute total normalization of blaze function
        blaze_intensity = np.empty(0) #intensity
        pixel = np.arange(2048)+1  #the fit to the blaze function is defined on this pixel array
        for det in range(self.blaze_coeff.shape[0]):
            for order in range(self.blaze_coeff.shape[1]):
                blaze_intensity = np.append(blaze_intensity,np.polyval(self.blaze_coeff[det, order],pixel))

        self.blaze_norm = np.max(blaze_intensity)
        
    def GetLimbDarkCoefficients(self, Teff, logg, Z, band):
        #Load  coefficients from Claret et al. 2011
        ldc_all = pd.read_csv(self.planet.paths['limbdark_coeff'], sep=';', skiprows=[i for i in range(0,47)]+[48,49,50])

        #select entries with correct spectral band
        ldc = ldc_all.loc[np.where(np.array([filt.strip() for filt in ldc_all.Filt]) == band)].reset_index()

        if len(ldc) == 0:
            print('Warning: No limb darkening coeff. found for spectral band {}!'.format(band))

        close_T = abs(ldc.Teff - Teff) == np.min(abs(ldc.Teff - Teff))
        close_logg = abs(ldc.logg[close_T] - logg) == np.min(abs(ldc.logg[close_T] - logg))
        close_Z = abs(ldc.Z[close_T][close_logg] - Z) == np.min(abs(ldc.Z[close_T][close_logg] - Z))

        chosen_index = np.arange(len(ldc))[close_T][close_logg][close_Z]

        self.ld_coeff = ldc.loc[chosen_index[0]]

    def Observe(self, nr, savepath = None):
        if len(self.opportunities)-1 < nr:
            raise Exception('Opportunity not found. Try an index between 0 and {}.'.format(len(self.opportunities)-1))
        
        tmid = self.opportunities.Tmid[nr]
        
        #Estimate SNR of one exposure
        if self.order_setting[0] == 'J':
            snr = ComputeSNR(self.planet.Jmag, self.exposure_time * day2sec, self.opportunities.Airm_mid[nr],
               method='Crires', mband='J')            
        elif self.order_setting[0] == 'K':
            snr = ComputeSNR(self.planet.Kmag, self.exposure_time * day2sec, self.opportunities.Airm_mid[nr],
               method='Crires', mband='K')
        else:
            raise Exception('For this spectral band no SNR estimation is available.')
        
        #Set up times
        if self.planet.obs_type == 'transmission':
            #desired out of transit observation time: 1h for short transits, 0.5Tdur for long transits
            if self.planet.Tdur > 2 / 24:
                outofTransit = self.planet.Tdur / 2
            else:
                outofTransit = 1 /24

            if self.opportunities.GoodCond_before[nr]/24 > outofTransit:
                obs_start = self.opportunities.Trans_start[nr] - outofTransit
            else:
                obs_start = self.opportunities.Trans_start[nr] - self.opportunities.GoodCond_before[nr]/24

            if self.opportunities.GoodCond_after[nr]/24 > outofTransit:
                obs_end = self.opportunities.Trans_end[nr] + outofTransit
            else:
                obs_end = self.opportunities.Trans_end[nr] + self.opportunities.GoodCond_after[nr]/24

            t = np.arange(obs_start, obs_end, self.timestep)
            n = len(t)

            #Shift so that the observation starts at a random time around obs_start
            t += np.random.rand()*self.timestep
            t_rel = t - tmid

            #Concrete orbit of planet
            #advance the orbit such that t=0 faces us for any omega
            omega = 0 #argument of periapsis, can be changed in the future
            t_sim = t_rel - self.planet.Period/4 - self.planet.Period*omega/360
        
        elif self.planet.obs_type == 'emission':
            t_ecl = self.opportunities.T_sececl[nr]
        
            obs_start = tmid - self.opportunities.GoodCond_before[nr]/24
            obs_end = tmid + self.opportunities.GoodCond_after[nr]/24
            t = np.arange(obs_start, obs_end, self.exposure_time)
            n = len(t)

            #Shift so that the observation starts at a random time around obs_start
            t += np.random.rand()*self.timestep
            t_rel = t - tmid

            #Concrete orbit of planet
            #advance the orbit such that t=0 faces us for any omega
            omega = 0 #argument of periapsis, can be changed in the future
            t_sim = t_rel + self.planet.Period/4 - self.planet.Period*omega/360 + (self.opportunities.Tmid[nr] - t_ecl) 
            
        self.t = t
        
        pos_pl = self.planet.Orbit_p.xyzPos(t_sim)
        vel_pl = self.planet.Orbit_p.xyzVel(t_sim)
        pos_pl = np.swapaxes(pos_pl, 0, 1)    #swap axis so that first component is spatial dimension
        vel_pl = np.swapaxes(vel_pl, 0, 1)

        rv_pl = vel_pl[2]*self.planet.Rs*au2sunr*au/1000/day2sec  #radial velocity in km/s

        #Concrete orbit of star
        pos_st = -self.planet.Orbit_s.xyzPos(t_sim)
        vel_st = -self.planet.Orbit_s.xyzVel(t_sim)    #minus so that star moves opposite to planet
        pos_st = np.swapaxes(pos_st, 0, 1)    #swap axis so that first component is spatial dimension
        vel_st = np.swapaxes(vel_st, 0, 1)    

        rv_st = vel_st[2]*self.planet.Rs*au2sunr*au/1000/day2sec  #radial velocity in km/s

        #velocity semi-amplitude K in km/s, calculated from absolute velocity and inclination
        K_pl = (np.sqrt(vel_pl[0][0]**2+vel_pl[1][0]**2+vel_pl[2][0]**2)*self.planet.Rs*au2sunr*au/1000/day2sec * 
               np.sin(self.planet.Inc/180*np.pi))
        K_st = (np.sqrt(vel_st[0][0]**2+vel_st[1][0]**2+vel_st[2][0]**2)*self.planet.Rs*au2sunr*au/1000/day2sec * 
               np.sin(self.planet.Inc/180*np.pi))
        
        #compute altitude of star at each time
        altaz = pyasl.eq2hor(t, np.ones(t.size) * self.planet.Ra, np.ones(t.size) * self.planet.Dec,
                                      lon=self.lon, lat=self.lat,alt=self.alt)
        altitude = altaz[0]

        #calculate airmass with plane parallel approx.
        airm = pyasl.airmassPP(90-altitude)

        #compute lightcurve
        lightcurve = np.ones(n)
        r = np.sqrt(pos_pl[0]**2 + pos_pl[1]**2)

        intersect_pl_st = IntersectArea(1, self.planet.Rp, r) #area of planet in front of star

        planetBehindStar = np.where(pos_pl[2] > 0)[0]
        planetBeforeStar = np.where(pos_pl[2] <= 0)[0]

        ld_realisitic = []

        for ind in planetBeforeStar:    
            if r[ind] > (1 + self.planet.Rp):
                ld_realisitic.append(0)
            else:
                n_points = 5000
                x = np.random.rand(n_points) * 2 * self.planet.Rp + pos_pl[0][ind] - self.planet.Rp
                y = np.random.rand(n_points) * 2 * self.planet.Rp + pos_pl[1][ind] - self.planet.Rp

                dist_to_planet = np.sqrt((pos_pl[0][ind] - x)**2 + (pos_pl[1][ind] - y)**2)
                dist_to_star = np.sqrt(x**2 + y**2)

                in_intersection = np.where((dist_to_planet <= self.planet.Rp) & (dist_to_star <= 1))

                if len(in_intersection[0]) == 0:
                    #intersection is so small that no points fall in the region
                    ld_realisitic.append(ld(1, self.ld_coeff))
                else:
                    ld_realisitic.append(np.mean(ld(dist_to_star[in_intersection], self.ld_coeff)))
        ld_realisitic = np.array(ld_realisitic)

        #remove intersected portion if planet is in front of star (not at secondary eclipse)
        lightcurve[planetBeforeStar] -= (intersect_pl_st[planetBeforeStar] * ld_realisitic) / self.total_int
        self.lightcurve = lightcurve
        
        #relative amount of planetary disk area visible (not oculted by star)
        planetdiskVisible = np.ones(n)
        planetdiskVisible[planetBehindStar] -= intersect_pl_st[planetBehindStar] / (self.planet.Rp**2*np.pi)
        
        #fraction of planetary disk that is in front of stellar disk
        ratio_intersected_area = intersect_pl_st[planetBeforeStar] / (np.pi * self.planet.Rp**2)
    
    
        ###Calculate barycentric redshift
        v_bary = np.array([pyasl.helcorr(self.lon, self.lat, self.alt, self.planet.Ra, self.planet.Dec, time)[0] for time in t])

        #v_bary is positive if the earth moves towards the star, so it has to be subtracted
        v_total = -v_bary + self.planet.RadVelStar

        if abs(np.mean(v_total)) < 2:
            print('Warning: the total velocity difference between system and observer is only {} km/s!'.format(round(np.mean(v_total),2)))
            
        ###Add redshift due to stars motion around barycenter and barycentric velocity
        spectrum_s_shift = np.ones([n, len(self.spectrum_s)])
        for i in range(n):
            spectrum_s_shift[i,:] = np.interp(self.wavelength_s, (1 + (v_total[i]+rv_st[i])/c)*self.wavelength_s, self.spectrum_s)
            
        ###Add redshift due to planets motion around barycenter and barycentric velocity,
        #and resample to wl grid of stellar spectrum
        spectrum_p_shift = np.ones([n, len(self.spectrum_s)])
        for i in range(n):
            spectrum_p_shift[i,:] = np.interp(self.wavelength_s, (1 + (v_total[i]+rv_pl[i])/c)*self.wavelength_p, self.spectrum_p)
            
            
        if self.planet.obs_type == 'transmission':
            #calculate eclipsed spectrum of stellar region behind planet
            #(but it is not yet scaled for the size of the planet)
            spectrum_ecl = np.dot((ratio_intersected_area * ld_realisitic)[:, None], self.spectrum_s[None, :])

            #compute Doppler velocity (from stellar rotation) and shift stellar spectrum
            vproj = pos_pl[0] * self.planet.Vsini
            for i in range(n):
                spectrum_ecl[i,:] = np.interp(self.wavelength_s, (1 + vproj[i]/c)*self.wavelength_s, spectrum_ecl[i,:])

                
            #imprint planetary spectrum
            #here the size of opaque core and wavelength dependent atmosphere is in spectrum_p_shift,
            #which is a radius ratio
            #np.pi/ti acounts for the difference in limbdarkening (mean LD of entire stellar disk)
            #spectrum_p_shift is the radius ratio
            spectrum_comb = spectrum_s_shift - spectrum_ecl * np.pi/self.total_int * spectrum_p_shift**2
            
        elif self.planet.obs_type == 'emission':
            #scale with total iluminated area and hide planet behind star
            planet_emission = np.pi * spectrum_p_shift * self.planet.Rp**2 * planetdiskVisible[:, None]

            #scale with sin to account for day/night side
            planet_emission *= ((np.sin(2*np.pi/self.planet.Period * (t-t_ecl) + 1/2*np.pi)+1)*0.5)[:, None]

            #add planetary emission to stellar spectrum
            spectrum_comb = spectrum_s_shift + planet_emission
            
        #apply extinction
        spectrum_comb *= np.dot(10**(-self.extcof/2.5 * (airm))[:, None], np.ones(self.wavelength_s.size)[None, :])

        
        #compute telluric spectra
#        #old way
#         telluric_array = np.zeros([n, self.wavelength_s.size])
#         for i in range(n):
#             telluric_amplified = Mag_spectrum(self.spectrum_tell, airm[i])
#             telluric_array[i,:] = np.interp(self.wavelength_s, self.wavelength_tell, telluric_amplified)
            
        #new way
        telluric_array = np.zeros([n, self.wavelength_s.size])
        telluric_resample = np.interp(self.wavelength_s, self.wavelength_tell, self.spectrum_tell)
        for i in range(n):
            telluric_array[i,:] = Mag_spectrum(telluric_resample, airm[i])

        #apply telluric atmosphere
        spectrum_comb *= telluric_array
        
        #instrument broadening
        spectrum_comb = Instrbrod(self.wavelength_s, spectrum_comb, self.resolution)

        #convert from energy per wavelength bin to photons per pixel
        spectrum_comb *= (np.dot(np.ones(n)[:, None], self.wavelength_s[None, :]) / np.min(self.wavelength_s))**2

        #RESAMPLING
        ws = np.min(self.wavelength_s)
        we = np.max(self.wavelength_s)

        dwlgrid = np.log(1 + 1/(self.resolution*self.sampling))
        nwlgrid = int((np.log(we) - np.log(ws))/dwlgrid)
        wlgrid = np.log(ws) + np.arange(nwlgrid)*dwlgrid
        spectrum_res = np.empty([n, wlgrid.size])
        for i in range(n):
            spectrum_res[i,:] = np.interp(wlgrid, np.log(self.wavelength_s), spectrum_comb[i,:])
        wavelength_res = np.exp(wlgrid)
        

        #Divide spectrum for each detector
        det1_wave = []
        det1_spec= []
        det2_wave = []
        det2_spec= []
        det3_wave = []
        det3_spec= []

        n_orders = self.crformat.shape[1]
        for i in range(n_orders):
            index1 = np.where((wavelength_res > self.crformat[0,i,0]) & (wavelength_res < self.crformat[0,i,1]))[0]
            if index1.size != 0:
                det1_wave.append(wavelength_res[index1])
                det1_spec.append(spectrum_res[:, index1])
            else: #if the order is not included in the chosen wavelength range, set it to nan
                det1_wave.append(np.ones([1500])*np.nan)
                det1_spec.append(np.ones([n,1500])*np.nan)

            index2 = np.where((wavelength_res > self.crformat[1,i,0]) & (wavelength_res < self.crformat[1,i,1]))[0]
            if index2.size != 0:
                det2_wave.append(wavelength_res[index2])
                det2_spec.append(spectrum_res[:, index2])
            else:
                det2_wave.append(np.ones([1500])*np.nan)
                det2_spec.append(np.ones([n,1500])*np.nan)

            index3 = np.where((wavelength_res > self.crformat[2,i,0]) & (wavelength_res < self.crformat[2,i,1]))[0]
            if index3.size != 0:
                det3_wave.append(wavelength_res[index3])
                det3_spec.append(spectrum_res[:, index3])
            else:
                det3_wave.append(np.ones([1500])*np.nan)
                det3_spec.append(np.ones([n,1500])*np.nan)
                

        #Add blaze function
        det1_wave_blazed, det1_spec_blazed = ApplyBlaze(det1_wave, det1_spec, self.blaze_wave[0], self.blaze_coeff[0], self.blaze_norm)
        det2_wave_blazed, det2_spec_blazed = ApplyBlaze(det2_wave, det2_spec, self.blaze_wave[1], self.blaze_coeff[1], self.blaze_norm)
        det3_wave_blazed, det3_spec_blazed = ApplyBlaze(det3_wave, det3_spec, self.blaze_wave[2], self.blaze_coeff[2], self.blaze_norm)
        
        #Scale with SNR value
        #scale such that the mean snr of the central 3 orders is equal to the calculated snr value
        central_exposure_id = int(n/2)
        #Get mean intensity of 3 central detector orders
        central_orders = []
        for i in range(1,4):
            central_orders.extend(det1_spec_blazed[i][central_exposure_id])
            central_orders.extend(det2_spec_blazed[i][central_exposure_id])
            central_orders.extend(det3_spec_blazed[i][central_exposure_id])

        mean_snr_central = np.mean(np.sqrt(central_orders))
        snr_scalefactor = snr**2 / mean_snr_central**2

        det1_spec_blazed *= snr_scalefactor
        det2_spec_blazed *= snr_scalefactor
        det3_spec_blazed *= snr_scalefactor
        
        self.test_wavelength = np.stack([det1_wave_blazed, det2_wave_blazed, det3_wave_blazed])
        self.test_spectrum = np.stack([det1_spec_blazed, det2_spec_blazed, det3_spec_blazed])
        
        #add noise
        for order in range(n_orders):
            for n_step in range(n):
                det1_spec_blazed[order][n_step] += np.sqrt(det1_spec_blazed[order][n_step])*np.random.normal(size=det1_spec_blazed[order][n_step].shape)
                det2_spec_blazed[order][n_step] += np.sqrt(det2_spec_blazed[order][n_step])*np.random.normal(size=det2_spec_blazed[order][n_step].shape)
                det3_spec_blazed[order][n_step] += np.sqrt(det3_spec_blazed[order][n_step])*np.random.normal(size=det3_spec_blazed[order][n_step].shape)

        self.obs_wavelength = np.stack([det1_wave_blazed, det2_wave_blazed, det3_wave_blazed])
        self.obs_spectrum = np.stack([det1_spec_blazed, det2_spec_blazed, det3_spec_blazed])
        
        #Save results
        if savepath != None:

            path_save_hdu = os.path.join(savepath, 'Spectrum_{}'.format(str.zfill(str(nr),2)))
            pathlib.Path(path_save_hdu).mkdir(parents=True, exist_ok=True)

            for i in range(n):
                wave_stacked = np.stack([det1_wave_blazed, det2_wave_blazed, det3_wave_blazed])
                spectrum_stacked = np.stack([det1_spec_blazed[:, i, :], det2_spec_blazed[:, i, :], det3_spec_blazed[:, i, :]])


                hdu_list = fits.HDUList([
                    fits.PrimaryHDU(),
                    fits.ImageHDU(wave_stacked),
                    fits.ImageHDU(spectrum_stacked), 
                ])

                header = hdu_list[0].header
                header['exposure'] = (self.exposure_time*24*3600, 'exposure time for each step [s]')
                header['readout'] = (self.readout_time*24*3600, 'readout time between exposures [s]')
                header['K_star'] = (np.round(K_st,4), 'radial vel semiamplitude of the star [km/s]')
                header['incl'] = (self.planet.Inc, 'inclination of planetary orbit [degrees]')
                header['o_set'] = (self.order_setting, 'Crires order setting')
                header['P'] = (self.planet.Period, 'Orbital period [days]')
                header['res'] = (self.resolution, 'spectral resolving power')
            #    header['vsini'] = (vsini, 'vsini of host star [km/s]')
                header['v_sys'] = (self.planet.RadVelStar, 'rad. v of stellar system barycenter [km/s]')

                hdu_list.writeto(os.path.join(path_save_hdu, self.planet.name.replace(" ", "") + '_' + str.zfill(str(i),3)+'.fits'), overwrite=True)

            #Save additional information
            path_save_adddata = savepath
            additional_data = pd.DataFrame(np.transpose([airm, t, -v_bary]), columns=['airmass', 'time', 'barycentric velocity'])
            additional_data.to_csv(os.path.join(path_save_adddata, self.planet.name.replace(" ", "") + '_additional_data_' + str.zfill(str(nr),2)+'.csv'), index=False)