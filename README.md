# ObsTools

For now, just a code for me to approximately plan my observing night.

---

**Disclaimer:** Error handling is very premature. The program may quit suddenly if the wrong input is selected. Please let me know about these, and I will eventually code in all the right catches. 

---
### What it needs:
  - A list of targets (Simbad-resolvable names)
  - Corresponding list of exposure times (stay tuned for exposure time calculators)
  - Site location or name, e.g., 'lco'
  - Year, month, and day (local) that your observations start
  
### What it gives:
  - A visibility plot of your targets that you can interact with on the command line
  - Experiment with your targets' placement throughout the night to maximize your time efficiency
  - Add targets, undo, print, and even go on to the next consecutive night interactively
  - Prints an approximate observing plan based on your input

![Here is a demo image:](https://raw.githubusercontent.com/eholmbeck/ObsTools/main/data/planner_demo.png)

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
