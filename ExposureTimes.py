import numpy as np
from scipy import integrate

# ===========================================
# MIKE Exposure Time Calculator
#
# This code is still under development! I often have to multiply by exptimes by 3
# ===========================================

# Just an interpolation that I don't want to recalculate every time
MIKE_mag = None

# Find fraction of light that's actually in the slit.
# By fault, assume a 5.0" slit in length
def gauss2d(slit, seeing, slit_length=5.0):
	sigma = seeing/2.355
	f = lambda y, x: np.exp(-(x**2 + y**2)/(2.*(sigma**2)))
	norm = 2.*np.pi*(sigma**2)
	return integrate.dblquad(f, -slit/2., slit/2, lambda x: -slit_length/2.0, lambda x: slit_length/2.0)[0]/norm
	
def read_in_efficiencies():
	import os
	from scipy.interpolate import interp1d
	path = os.path.dirname(os.path.realpath(__file__))
	MIKE_EFF_wave,MIKE_EFF_mag = np.genfromtxt(path+'/data/MIKE_efficiency.csv', delimiter=',').T
	MIKE_mag = interp1d(MIKE_EFF_wave,MIKE_EFF_mag)
	return MIKE_mag

def MIKE_etc(snr, centroid, vmag, slit, seeing=0.5, binning=1.0, gain=0.47, **kwargs):
	# Reset globals so we don't have to do this every time
	global MIKE_mag
	if MIKE_mag is None: MIKE_mag = read_in_efficiencies()
	
	# TODO: add an empirical conversion since the seeing is not the same vs. wavelength and airmass
	#seeing_at_centroid = np.max([10.**(-1.3227e-3*centroid + 6.1686)-0.47, 0])+np.max([seeing, 0.47])
	frac = gauss2d(slit,seeing, **kwargs)
	
	slit_res = 65000.
	if centroid < 5000: slit_res = 83000.
	count_rate = gain * binning * (centroid/slit_res) * 10**(0.4*(MIKE_mag(centroid)-vmag))*frac
	return (snr**2) / (count_rate)

# ===========================================
# Query to get Vmags
# Useful for MIKE exposure times

def get_vmags(targets_by_name):
	from astroquery.simbad import Simbad
	Simbad.add_votable_fields("flux(V)")
	results = Simbad.query_objects(targets_by_name)
	return list(results['FLUX_V'])

# ===========================================
# Query to get pms, Vmags, dec, and ra
# Useful for writing catalogs

def get_Simbad_data(targets_by_name):
	from astroquery.simbad import Simbad
	Simbad.add_votable_fields("flux(V)")
	Simbad.add_votable_fields("pmdec")
	Simbad.add_votable_fields("pmra")
	Simbad.add_votable_fields("flux(g)")
	Simbad.add_votable_fields("pm_bibcode")
	Simbad.add_votable_fields("fe_h")
	results = Simbad.query_objects(targets_by_name)
	return list(results['FLUX_V'])

# ===========================================

