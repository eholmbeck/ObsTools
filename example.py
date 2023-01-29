from ObservingStrategy import Plan

# To use, simply add your target list (Simbad-resolvable names) and a corresponding list of exposure times
targets = ['HR1190', '2MASSJ06290787-6709523', 'HD122563']
exptimes = [120, 3600, 600]

# ===============================================
obs = Plan(site='lco', year=2023, month=1, day=28, targets=targets, exptimes=exptimes)
