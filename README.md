# ObsTools

For now, just a code for me to approximately plan my observing night.

### What it needs:
  - A list of targets (Simbad-resolvable names)
  - Corresponding list of exposure times (stay tuned for exposure time calculators)
  - Site location or name, e.g., 'lco'
  - Year, month, and day (local) that your observations start
  
### What it gives:
  - A visibility plot of your targets that you can interact with on the command line
  - Experiment with your targets' placement throughout the night to maximize your time efficiency
  - Add targets, undo, print, and even go on to the next consecutive night
  - Prints an approximate observing plan based on your input
  
### How to use it:
First install these:
  - `matplotlib`
  - `numpy`
  - `astropy` (there may be a version problem here)
  - `timezonefinder`
  - `pytz`
  
Then call the planner with:
  `obs = Plan(site='lco', year=2023, month=1, day=28, targets=targets, exptimes=exptimes)`.

There is an example file included to show you how it works. Just call `python example.py`!

**Disclaimer:** I haven't yet coded this to be a nice Python package, so you may have to contend with relative/full paths.
