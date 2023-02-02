import numpy as np
from PyAstronomy import pyasl
from scipy import integrate

def ld(r, coeff):
    #limb darkening as function of radial distance r to center
    #coeff is a dataframe containing the Limb Darkening Coefficients
    ld_set = [coeff.a1, coeff.a2, coeff.a3, coeff.a4] 
    mu = np.sqrt(1 - r**2)
    ld = 1
    
    for k in range(1,5):
        ld -= ld_set[k-1]*(1 - mu**(k/2))
    return ld

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
    
    
# distance planet disk centre - stellar disk centre
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

#fraction of planetary disk that is in front of stellar disk
ratio_intersected_area = intersect_pl_st[planetBeforeStar] / (np.pi * self.planet.Rp**2)
    
#compute total intensity
nint = 1000
rgrid = np.arange(nint)/(nint-1)
self.total_int = integrate.simps(2*np.pi*rgrid*ld(rgrid, self.ld_coeff), rgrid)
    
#calculate eclipsed spectrum of stellar region behind planet
#(but it is not yet scaled for the size of the planet)
spectrum_ecl = np.dot((ratio_intersected_area * ld_realisitic)[:, None], self.spectrum_s[None, :])

#np.pi/total_int acounts for the difference in limbdarkening (mean LD of entire stellar disk)
spectrum_ecl *= np.pi/self.total_int

#compute Doppler velocity (from stellar rotation) and shift stellar spectrum
vproj = pos_pl[0] * self.planet.Vsini
for i in range(n_divided):
    spectrum_ecl[i,:] = np.interp(self.wavelength_s, (1 + vproj[i]/c)*self.wavelength_s, spectrum_ecl[i,:])

    
#imprint planetary spectrum
#here the size of opaque core and wavelength dependent atmosphere is in spectrum_p,
#which is a radius ratio
#spectrum_p is the radius ratio
spectrum_comb = spectrum_s - spectrum_ecl * spectrum_p**2
