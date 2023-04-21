from ObservingStrategy import Plan

# To use, simply add your target list (Simbad-resolvable names) and a corresponding list of exposure times
targets = [#'my_star', 
		   'HR1190', '2MASSJ05514213-3327337', 
		   '2MASSJ05521578-3953184', '2MASSJ06290787-6709523',
		   '2MASSJ10160492-3955512',
		   'HD122563']

exptimes = [120, 3600, 1200, 900, 1800, 60, 1200]
coords = ['12:45:19.5 -47:51:06.9']

site = 'lco'
year = 2023
month = 1
day = 29

# Or use the exposure-time calculator - STILL UNDER DEVELOPMENT
from ExposureTimes import MIKE_etc,get_vmags

snr = 150
centroid = 4000
slit = 0.5
seeing = 0.7
binning = 1.0

exptimes = [2*MIKE_etc(snr, centroid, vmag, slit, binning) for vmag in get_vmags(targets)]

# ===============================================
obs = Plan(site='lco', year=year, month=month, day=day, targets=targets, exptimes=exptimes, coords=coords)
