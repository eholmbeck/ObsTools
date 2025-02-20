"""
To use, call this file (or your own, following 
the basic steps below.

	> python example.py
"""

# ===============================================
# There are two ways to use this:
#  1. Input your own exposure times -or-
#  2. Use my MIKE one
#
# get_vmags queries simbad to get the V magnitude,
# but you can also just input your own vmags.
# ===============================================

targets = [#'my_star', 
		   'HD122563', 'HR5724', '2MASS J15582962-1224344', '2MASSJ17273886-5133298',
		   'RAVE J192819.9-633935', '2MASS J22190836-2333467', 
		   'BPS CS 22945-0017', '2MASSJ01425445-0904162'
		   ]

exptimes = [120, 3600, 1200, 900, 1800, 60, 1200]

# Or use the exposure-time calculator - STILL UNDER DEVELOPMENT
from ExposureTimes import MIKE_etc,get_vmags

snr = 120
centroid = 4000
slit = 0.5
seeing = 0.7
binning = 2.0

# In my experience, I need twice as much exposure as MIKE calculates
exptimes = [2*MIKE_etc(snr, centroid, vmag, slit, binning) for vmag in get_vmags(targets)]

# ===============================================
#           Site and date
# ===============================================

#coords = ['12:45:19.5 -47:51:06.9']

site = 'lco'
year = 2023
month = 7
day = 18

# ===============================================
#  Interactively plan observations throughout the night
#  (highly calibrated for MIKE)
# ===============================================
from ObservingStrategy import Plan

obs = Plan(site='lco', year=year, month=month, day=day, targets=targets, exptimes=exptimes)#, coords=coords)
