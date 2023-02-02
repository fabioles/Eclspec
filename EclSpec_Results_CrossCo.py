# For Eclspec simulations, the orbit has to be circular!

#%% Planet, model and nights
import sys

sys.path.append('/Users/fabio/Data/PhD/Projects/AtmAnaTools')

import numpy as np
import matplotlib.pyplot as plt

import CrossCo as crossco
import excali as exc

data_directory = '/Users/fabio/Data/Master/Eclspec/Simulations/Results/220810/LTT1445Ab/LTT1445Ab_SNR100_EarthAtmosphere_Smear20'
nexa_path = '/Users/fabio/Data/PhD/Projects/Data/NexaComposite.csv'


planet = exc.Planet('LTT 1445 A b', nexa_path)
planet.ecc = 0
planet.periastron_arg_degrees = 0
planet.trans_epoch = 2459489.827186169
# planet.trans_epoch = 2458423.426290
planet.period = 5.35876570
# planet.period = 5.35882000


# path_model = ('/Users/fabio/Data/Master/Eclspec/Simulations/'
#               + 'LTT1445Ab_Trans_H2O_CH4_CO2.txt')
path_model = ('/Users/fabio/Data/Master/Eclspec/Atmosphere_models/'
              + 'Earth-transmisison-CO2-CH4-H2O.txt')
model = crossco.Model()
# model.LoadModelSpectrum(path_model, 'angstrom')
model.LoadModelSpectrum(path_model, 'micron', 'km', planet)
model.NormalizeModel('CRIRES')


def RunCC(cleaner, model, at_negative_Kp=False, plotAll=False):
    cc = crossco.CrossCorrelator(cleaner, model)
    
    if at_negative_Kp:
        cc.AnalyzeAllSysremIterations(-200, 200, 1.0, -400, 0, True)
        if plotAll:
            cc.plot.DetMapAllIterations(-75,75,-400,0,True)
    
    else:
        cc.AnalyzeAllSysremIterations(-200, 200, 1.0, 0, 400)
        if plotAll:
            cc.plot.DetMapAllIterations()
            
    return cc


#%%
directory = ('/Users/fabio/Data/Master/Eclspec/Simulations/Results/' 
             +'221128/LTT1445Ab_SNR100_EarthAtmosphere_Smear10_strongerInjected30')
# directory = data_directory
index = 0
night0 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 1
night1 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 2
night2 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 3
night3 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 4
night4 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 5
night5 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 6
night6 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 7
night7 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 8
night8 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)
index = 9
night9 = exc.Observation(planet, directory, data_type = 'EclSpec', index = index)

# from Data_extract import CARMENES_extract, CRIRES_extract, EclSpec_extract
# data = EclSpec_extract(directory, index = index)
#%%%

# night = night0
# data = np.empty([3,18,night.n_files,2048])
# data[0] = night.wave
# data[1] = night.spectrum
# data[2] = night.error
# np.save('/Users/fabio/Desktop/LTT1445Ab_300Injected/'+
#         'LTT1445Ab_SNR100_EarthAtmosphere_Smear10_strongerInjected300_0.npy',data)  


# import pickle
# with open('/Users/fabio/Desktop/LTT1445Ab_300Injected/LTT1445Ab_SNR100_EarthAtmosphere_Smear10_strongerInjected300.pkl', 'wb') as f:
#     pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
    
# import pickle

# df = dict()
# with open('/Users/fabio/Desktop/LTT1445Ab_SNR100_EarthAtmosphere_Smear10_strongerInjected300.pkl', 'rb') as f:
#     df = pickle.load(f)
    
# test = np.load('/Users/fabio/Desktop/LTT1445Ab_SNR100_EarthAtmosphere_Smear10_strongerInjected300.npy') 
#%% Setup

injectionStrength = 0
at_negative_Kp = False
plotAll = True
binning_method = 'smoothing'
#%%
cleaner0 = crossco.SpectrumCleaner(night0)
cleaner0.Normalize(binning_method)
cleaner0.RunSysrem(5)

cc0 = RunCC(cleaner0, model, at_negative_Kp, plotAll = plotAll)
# cc0.plot.DetMap(-30, 30, 0, 250)
#%%
cleaner1 = crossco.SpectrumCleaner(night1)
cleaner1.Normalize(binning_method)
cleaner1.RunSysrem(5)

cc1 = RunCC(cleaner1, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner2 = crossco.SpectrumCleaner(night2)
cleaner2.Normalize(binning_method)
cleaner2.RunSysrem(5)

cc2 = RunCC(cleaner2, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner3 = crossco.SpectrumCleaner(night3)
cleaner3.Normalize(binning_method)
cleaner3.RunSysrem(5)

cc3 = RunCC(cleaner3, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner4 = crossco.SpectrumCleaner(night4)
cleaner4.Normalize(binning_method)
cleaner4.RunSysrem(5)

cc4 = RunCC(cleaner4, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner5 = crossco.SpectrumCleaner(night5)
cleaner5.Normalize(binning_method)
cleaner5.RunSysrem(5)

cc5 = RunCC(cleaner5, model, at_negative_Kp, plotAll = plotAll)
# cc0.plot.DetMap(-30, 30, 0, 250)
#%%
cleaner6 = crossco.SpectrumCleaner(night6)
cleaner6.Normalize(binning_method)
cleaner6.RunSysrem(5)

cc6 = RunCC(cleaner6, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner7 = crossco.SpectrumCleaner(night7)
cleaner7.Normalize(binning_method)
cleaner7.RunSysrem(5)

cc7 = RunCC(cleaner7, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner8 = crossco.SpectrumCleaner(night8)
cleaner8.Normalize(binning_method)
cleaner8.RunSysrem(5)

cc8 = RunCC(cleaner8, model, at_negative_Kp, plotAll = plotAll)
#%%
cleaner9 = crossco.SpectrumCleaner(night9)
cleaner9.Normalize(binning_method)
cleaner9.RunSysrem(5)

cc9 = RunCC(cleaner9, model, at_negative_Kp, plotAll = plotAll)

#%% Combining 
cc_list  = [cc0, cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8, cc9]
comb = crossco.CrossCorrelatorMultipleNights(cc_list, model)

comb.MergeData([5, 5, 5, 5, 5, 5, 5, 5, 5, 5])
comb.CalcDetMap(Kp_trial_min = 0, Kp_trial_max = 400)
comb.FindDetMapMax()
comb.plot.DetMap(-60, 60, 0, 400)
#%% Old stuff
path_model = ('/Users/fabio/Data/PhD/Projects/Data/AtmosphereModels/Hat-P-11b/'
              + 'Trans_H2O_CH4_740K_Chachan2019.txt')
model = crossco.Model()
model.LoadModelSpectrum(path_model, 'angstrom')

cc1 = crossco.CrossCorrelator(cleaner1, model)
cc1.AnalyzeAllSysremIterations(-300, 300, 1.0, 0, 400, at_negative_Kp = False)
cc0.plot.DetMapAllIterations(-30, 30, 0, 250, at_negative_Kp = False)

for i in range(night1.n_orders):
    cc4.plot.CCF(i, 10)
    
#%% Wasp 107
directory = '/Users/fabio/Data/Master/Eclspec/Simulations/Results/220718/Wasp107b/Wasp107b_HotJup_smearing_test/Wasp107b_SNR100_HotJupAtmosphere_withSmear10'
nexa_path = '/Users/fabio/Data/PhD/Projects/Data/NexaComposite.csv'


planet = exc.Planet('WASP-107 b', nexa_path)
planet.ecc = 0
planet.periastron_arg_degrees = 0
path_model = ('/Users/fabio/Data/Master/Eclspec/Atmosphere_models/'
              + 'Transmission-HD189733-CO-H2O.txt')
model = crossco.Model()
# model.LoadModelSpectrum(path_model, 'angstrom')
model.LoadModelSpectrum(path_model, 'micron', 'km', planet)
model.NormalizeModel('CRIRES')

night0w = exc.Observation(planet, directory, data_type = 'EclSpec', index = 0)

cleaner0w = crossco.SpectrumCleaner(night0w)
cleaner0w.Normalize(binning_method)
cleaner0w.RunSysrem(10)

cc0w = RunCC(cleaner0w, model, at_negative_Kp, plotAll = plotAll)
plt.plot([cc0w.Cohen_d(i) for i in range(2, 11)])
#%% Using corrected data from Ansgar

path = '/Users/fabio/Data/Master/Eclspec/Comparison_Ansgars_CC/SNR Paper-selected/corrected_data.npz'
data = np.load(path)

# data['night'].shape
# for item in data:
#     print(item)
    
index = np.where(data['night'] == 0)
wave = data['wave'][index]
flux = data['flux'][index]
errors = data['uncs'][index]

n_files = night0.n_files
wave = np.swapaxes(wave.reshape((n_files, 18, 2048)), 0, 1)
flux = np.swapaxes(flux.reshape((n_files, 18, 2048)), 0, 1)
errors = np.swapaxes(errors.reshape((n_files, 18, 2048)), 0, 1)

cleaner0.spectrumNorm = flux
cleaner0.errorNorm = errors

cleaner0.RunSysrem(20)
cc0 = RunCC(cleaner0, model, at_negative_Kp, plotAll = plotAll)

# plt.plot(wave[0], flux[0])

# for i in range(1):
#     plt.plot(wave1[i,0], flux1[i,0])
    # plt.plot(night0.wave[i,0], night0.spectrum[i,0])

#%%
from __future__ import division, print_function

import warnings

import numpy as np
from scipy import constants
from tqdm import tqdm

c_light = constants.c / 1e3

class Sysrem:
    def __init__(
        self,
        input_data: np.ndarray,
        errors: np.ndarray = None,
        iterations: int = 100,
        tolerance: float = 1e-10,
        spec: np.ndarray = None,
        return_spec: bool = False,
    ):

        self.epoch_dim, self.stars_dim = input_data.shape
        self.input_data = input_data
        if errors is None:
            errors = np.sqrt(np.abs(input_data))
        self.errors = errors

        self.iterations = iterations
        self.tolerance = tolerance
        self.return_spec = return_spec
        self.spec = spec

    def run(self, num: int):
        errors_squared = np.require(self.errors ** 2, requirements="C")
      
        # Create all the memory we need for SYSREM
        if self.spec is None:
            c = np.ones(self.stars_dim)
            c_loc = np.zeros(self.stars_dim)
            c_numerator = np.zeros(self.stars_dim)
            c_denominator = np.zeros(self.stars_dim)
        else:
            c = c_loc = self.spec

        a = np.ones(self.epoch_dim)
        a_loc = np.zeros(self.epoch_dim)
        a_numerator = np.zeros(self.epoch_dim)
        a_denominator = np.zeros(self.epoch_dim)

        # remove the median as the first component
        median = np.nanmedian(self.input_data, axis=0)
        residuals = self.input_data - median

        syserrors = [None] * (num + 1)
        syserrors[0] = median

        for n in tqdm(range(num), desc="Removing Systematic #", leave=False):
            previous_diff = np.inf

            with warnings.catch_warnings():
                warnings.simplefilter(action="ignore", category=RuntimeWarning)
                # minimize a and c values for a number of iterations, iter
                for i in tqdm(range(self.iterations), desc="Converging", leave=False):
                    # Using the initial guesses for each a value of each epoch, minimize c for each star
                    if self.spec is None:
                        np.nansum(
                            a[:, None] * residuals / errors_squared,
                            axis=0,
                            out=c_numerator,
                        )
                        np.nansum(
                            a[:, None] ** 2 / errors_squared, axis=0, out=c_denominator
                        )
                        np.divide(c_numerator, c_denominator, out=c_loc)

                    # Using the c values found above, minimize a for each epoch
                    np.nansum(
                        c_loc[None, :] * residuals / errors_squared,
                        axis=1,
                        out=a_numerator,
                    )
                    np.nansum(
                        c_loc[None, :] ** 2 / errors_squared,
                        axis=1,
                        out=a_denominator,
                    )
                    np.divide(a_numerator, a_denominator, out=a_loc)

                    diff = np.nanmean((c_loc - c) ** 2) + np.nanmean((a_loc - a) ** 2)
                    # Swap the pointers to the memory
                    c, c_loc = c_loc, c
                    a, a_loc = a_loc, a
                    if (
                        self.tolerance is not None
                        and diff < self.tolerance
                        or diff > previous_diff
                    ):
                        break
                    previous_diff = diff

            # Create a matrix for the systematic errors:
            # syserr = np.zeros((stars_dim, epoch_dim))
            syserr = c[None, :] * a[:, None]
            # Remove the systematic error
            residuals -= syserr

            if self.return_spec:
                syserrors[n + 1] = np.copy(c)
            else:
                syserrors[n + 1] = syserr

        return residuals, syserrors
    
sysrem = Sysrem(flux, errors)
residuals, syserrors = sysrem.run(10)
#%%
residuals = np.swapaxes(residuals.reshape((n_files, 18, 2048)), 0, 1)


cleaner0.spectrumRemoved[10] = residuals
cc0 = RunCC(cleaner0, model, at_negative_Kp, plotAll = plotAll)

#%%
plt.plot([37.6, 42.7, 46.1, 46.1, 46.1, 45.6, 45.4, 44.4, 43.3, 43.3])
# plt.plot([2.6, 4.0, 3.1, 4.7, 5.4, 4.9, 5.0, 5.1, 5.6, 5.9])
plt.xlabel('Number of added transits (added from best to worst)')
plt.ylabel('S/N')

# plt.plot(cc9.detmap[-1][:,200])
# plt.vlines(77.3, 16, 29, ls='dashed', alpha=0.5)
# plt.xlabel('Kp')
# plt.ylabel('S/N')
# plt.title('''Cut through Det-map in Kp direction. \n Dashed line is literature Kp.''')

# cc9.plot.CCF('sum', align=True)
# v_rv = cc9.planet.RadialVelocityPlanet(night9.time, Kp=400)
# v_total = v_rv + night9.planet.v_sys - night9.v_bary
# plt.plot(v_total, night9.phase, align=True)
import plotter
plotter.Plotter().MakePlotTransparent(plt.gcf())
#%%

plt.imshow(cc1.detmap_unscaled[-1], origin='lower', aspect='auto')
plt.plot(cc3.detmap_unscaled[-1][97])
plt.plot(cc1.detmap_unscaled[-1][97]+cc3.detmap_unscaled[-1][97])

detmap = np.append(cc1.detmap_unscaled[-1][:, std_mask], cc3.detmap_unscaled[-1][:, std_mask])
std_mask = (abs(cc1.rv_list) > 50)
std = np.nanstd(detmap)
mean = np.nanmean(detmap)

plt.plot(cc1.detmap_unscaled[-1][:, std_mask].flatten())
plt.plot(cc3.detmap_unscaled[-1][:, std_mask].flatten())

#%%
cc_list  = [cc0, cc1, cc]
std_mask = (abs(cc1.rv_list) > 50) & (abs(cc1.rv_list) < 150)
cc1.ccf_sum[-1] -= np.nanmean(cc1.ccf_sum[-1][:,std_mask])
std_mask = (abs(cc0.rv_list) > 50) & (abs(cc0.rv_list) < 150)
cc0.ccf_sum[-1] -= np.nanmean(cc0.ccf_sum[-1][:,std_mask])


comb = crossco.CrossCorrelatorMultipleNights(cc_list, model)

comb.MergeData([5, 5])
comb.CalcDetMap(Kp_trial_min = 0, Kp_trial_max = 400)
comb.FindDetMapMax()
comb.plot.DetMap(-60, 60, 0, 400)

v_total = cc0.obs.v_total
ccf_unaligned = np.copy(cc0.ccf_sum[-1])
ccf_align_0 = np.copy(cc0.ccf_sum[-1])
for spectrumId in range(cc0.n_files):
    ccf_align_0[spectrumId] = np.interp(cc0.rv_list, cc0.rv_list - v_total[spectrumId],
                                        ccf_unaligned[spectrumId])  
v_total = cc9.obs.v_total
ccf_unaligned = np.copy(cc9.ccf_sum[-1])
ccf_align_1 = np.copy(cc9.ccf_sum[-1])
for spectrumId in range(cc9.n_files):
    ccf_align_1[spectrumId] = np.interp(cc9.rv_list, cc9.rv_list - v_total[spectrumId],
                                        ccf_unaligned[spectrumId])  

v_total = comb.v_total
ccf_unaligned = np.copy(comb.ccf_sum[-1])
ccf_align_2 = np.copy(comb.ccf_sum[-1])
for spectrumId in range(comb.n_files):
    ccf_align_2[spectrumId] = np.interp(comb.rv_list, comb.rv_list - v_total[spectrumId],
                                        ccf_unaligned[spectrumId])  

std_mask = (abs(cc1.rv_list) > 50) & (abs(cc1.rv_list) < 150)

plt.plot(cc1.rv_list, np.nanmean(ccf_align_0, axis=0)-0.0171635, label='One transit 0')
plt.plot(cc1.rv_list, np.nanmean(ccf_align_1, axis=0)-0.276920, label='One transit')
plt.plot(cc1.rv_list, np.nanmean(ccf_align_2, axis=0), label='Two transits')
plt.legend()
plotter.Plotter().MakePlotTransparent(plt.gcf())

std_0 = std = np.nanstd(np.nanmean(ccf_align_1, axis=0)[std_mask])
std_1 = std = np.nanstd(np.nanmean(ccf_align_1, axis=0)[std_mask])
std_2 = std = np.nanstd(np.nanmean(ccf_align_2, axis=0)[std_mask])

max_1 = np.nanmax(np.nanmean(ccf_align_1, axis=0))
max_2 = np.nanmax(np.nanmean(ccf_align_2, axis=0))

for i in range(18):
    print(i)
    for cc in cc_list:
        # print(cc.detmap_unscaled[-1].max())
        print(cc.ccf[-1][i].max())

for cc in cc_list:
    print(cc.ccf_sum[-1].std())
    
    
#%% Welch t test

cc0.obs.v_total

inTransit_grid = np.repeat(cc0.obs.InTransit[:, np.newaxis], cc0.rv_list.size, axis = 1)
rv_grid = np.repeat(cc0.rv_list[np.newaxis], cc0.obs.n_files, axis = 0)
rv_grid -= cc0.obs.v_total[:,np.newaxis]

in_trail = (abs(rv_grid) < 5) & inTransit_grid
out_of_trail = (abs(rv_grid) > 25) & inTransit_grid

N_sysrem = 2

in_sample = cc0.ccf_sum[N_sysrem][in_trail]
out_sample = cc0.ccf_sum[N_sysrem][out_of_trail]

plt.hist(in_sample, density=True, alpha = 0.5, label = 'in trail')
plt.hist(out_sample, alpha = 0.5, density=True, label = 'out of trail')
plotter.Plotter().MakePlotTransparent(plt.gcf())
plt.legend()
#%%

def welch_t(a, b, ua=None, ub=None):
    # t = (mean(a) - mean(b)) / sqrt(std(a)**2 + std(b)**2)
    if ua is None:
        ua = a.std() / np.sqrt(a.size)
    if ub is None:
        ub = b.std() / np.sqrt(a.size)

    xa = a.mean()
    xb = b.mean()
    t = (xa - xb) / np.sqrt(ua ** 2 + ub ** 2)
    return t


def cohen_d(a, b):
    sa = a.std()
    sb = b.std()
    s = ((a.size - 1) * sa ** 2 + (b.size - 1) * sb ** 2) / (a.size + b.size - 2)
    s = np.sqrt(s)
    d = np.abs(a.mean() - b.mean()) / s
    return d

welch_t(in_sample, out_sample)

#%%
plt.plot([cc0.Welch_t(i) for i in range(11)])
plotter.Plotter().SetStyle(True, False)
plt.ylim(-10, 10)
plt.ylabel('Welch_t')
    


