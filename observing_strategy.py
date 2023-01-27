#from matplotlib import use
#use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time, TimeDelta

from scipy.interpolate import interp1d
from scipy import integrate

try: plt.style.use('/Users/eholmbeck/.config/matplotlib/stylelib/erika.mplstyle')
except: plt.style.use('/Users/holmbeck/.config/matplotlib/stylelib/erika.mplstyle')

# ===========================================
# MIKE Exposure Time Calculator

def gauss2d(slit, seeing):
	sigma = seeing/2.355
	f = lambda y, x: np.exp(-(x**2 + y**2)/(2.*(sigma**2)))
	norm = 2.*np.pi*(sigma**2)
	return integrate.dblquad(f, -slit/2., slit/2, lambda x: -2.5, lambda x: 2.5)[0]/norm
	
path = '/Users/holmbeck/Library/Mobile Documents/com~apple~CloudDocs/Observing'
MIKE_EFF_wave,MIKE_EFF_mag = np.genfromtxt(path+'/MIKE_efficiency.csv', delimiter=',').T
MIKE_mag = interp1d(MIKE_EFF_wave,MIKE_EFF_mag)

def etc(snr, vmag, slit, binning, gain):
	frac = gauss2d(slit,seeing)
	slit_res = 65000.
	if centroid < 5000: slit_res = 83000.
	count_rate = gain * binning * (centroid/slit_res) * 10**(0.4*(MIKE_mag(centroid)-vmag))*frac
	return (snr**2) / (count_rate)

# ===========================================

def get_vmags(targets_by_name):
	from astroquery.simbad import Simbad
	Simbad.add_votable_fields("flux(V)")
	results = Simbad.query_objects(targets_by_name)
	return list(results['FLUX_V'])

# ===========================================

# Useful to print out naut. and astr. twilight to the minute
def twilight(year=2000, month=1, day=1, site='lco', alt=-15):
	from scipy.optimize import fsolve
	if type(site) is str:
		site = EarthLocation.of_site(site)
	
	# Try at 6:00 pm
	run_start = {'year': year, 'month': month, 'day': day, 'hour': 18}
	time = Time(run_start, format='ymdhms', scale='local').jd
	
	def find_twilight(time, alt=alt):
		run_start = Time(time, format='jd')
		altazframe = AltAz(obstime=run_start, location=site)
		sunaltaz = get_sun(run_start).transform_to(altazframe)
		return sunaltaz.alt.value - alt
	
	dusk = fsolve(find_twilight, time, args=(alt))[0]
	dawn = fsolve(find_twilight, dusk+0.5, args=(alt))[0]
	
	return Time(dusk, format='jd'), Time(dawn, format='jd')

# ===========================================
'''
class Plan:
	def __init__(target):
		self.target = target
		self.airmas = self.
	
	def twilight(year=2000, month=1, day=1, site='lco'):
		if type(site) is str:
			site = EarthLocation.of_site(site)
	
		run_start = {'year': year, 'month': month, 'day': day}
		time = Time(run_start, format='ymdhms').jd
	
		def solve_dawn(time):
			run_start = Time(time, format='jd')
			altazframe = AltAz(obstime=run_start, location=site)
			sunaltaz = get_sun(run_start).transform_to(altazframe)
			return sunaltaz.alt.value - 15.0

		def solve_dusk(time):
			run_start = Time(time, format='jd')
			altazframe = AltAz(obstime=run_start, location=site)
			sunaltaz = get_sun(run_start).transform_to(altazframe)
			return sunaltaz.alt.value + 15.0
	
		dawn = fsolve(solve_dawn, time)[0]
		dusk = fsolve(solve_dusk, dawn+0.4)[0]
	
		return Time(dawn, format='jd'), Time(dusk, format='jd')
'''

targets = ['2MASSJ04520910-6209121', '2MASSJ05381700-7516207', '2MASSJ05514213-3327337', 
		   '2MASSJ05521578-3953184', '2MASSJ06290787-6709523', '2MASSJ07114252-3432368',
		   '2MASSJ10160492-3955512', '2MASSJ10191573-1924464', '2MASSJ10401894-4106124']

snr = 180
centroid = 4000
slit = 0.35
seeing = 0.5
binning = 1.0
gain = 0.47
site = 'lco'
year = 2023
month = 1
day = 28

exptimes = [3.*etc(snr, vmag, slit, binning, gain) for vmag in get_vmags(targets)]

site_loc = EarthLocation.of_site(site)

dusk_jd, dawn_jd = twilight(year=year, month=month, day=day, site=site_loc, alt=-12)
times = np.linspace(dusk_jd.value, dawn_jd.value, 100)
airmass = []

fig = plt.figure(figsize=(9,6))
ax = fig.add_axes([0.08,0.1,0.55,0.8])
#thismanager = plt.get_current_fig_manager()
#thismanager.window.wm_geometry("+900+0")

ax.set_ylim(0,90)
ax.set_yticks(range(0,91,10))
ax.set_yticklabels([r'{:2.0f}$^\circ$'.format(y) for y in ax.get_yticks()])

ax.set_ylabel('Altitude')
ax.set_xlabel('UTC')
ax2 = ax.twinx()
ax2.set_ylim(0,90)
ax2.set_yticks(range(10, 91,10))
ax2.set_yticklabels(['{:.2f}'.format(1.0/np.sin(y*np.pi/180.)) for y in ax2.get_yticks()])
ax2.text(1.01,0,'Airmass', va='center', ha='left', rotation=90, transform=ax2.transAxes)

night_length = (dawn_jd - dusk_jd).sec/3600.
hms = np.array(dusk_jd.iso.split()[-1].split(':'), dtype=float)
t0 = hms[0] + hms[1]/60. + hms[2]/3600.
decimal_times = np.linspace(t0, t0+night_length, 100)

ticks = np.arange(int(decimal_times[0])-1, int(decimal_times[-1])+2)
ax.set_xticks(ticks)
ax.set_xticklabels(['{:.0f}'.format(t) for t in ticks%24.])
ax.grid(True, ls=':', color='gray')

ax.set_xlim(decimal_times[0]-1, decimal_times[-1]+1)
ax.axvline(decimal_times[0], color='gray', alpha=0.6, ls='--')
ax.axvline(decimal_times[-1], color='gray', alpha=0.6, ls='--')

equal_times = len(decimal_times)/len(targets)
all_stars = []

for i,target in enumerate(targets):
	coords = SkyCoord.from_name(target)
	
	altaz = [coords.transform_to(AltAz(obstime=Time(time, format='jd'),location=site_loc))\
			 for time in times]
	
	#airmass.append([alt.secz.value for alt in altaz])
	airmass.append([alt.alt.value for alt in altaz])
	l = ax.plot(decimal_times, airmass[-1], color='gray', lw=0.7, alpha=0.7)
	all_stars.append(l)
	label_loc = int(equal_times*i)
	ax.text(decimal_times[label_loc], airmass[-1][label_loc], i, color='gray', fontsize='x-small', va='bottom')

plt.draw()
plt.pause(0.01)

for i in range(len(targets)):
	print('{:>3.0f}: {:<24} {:>5.0f} s ({:.2f}hr)'.format(i, targets[i], exptimes[i], exptimes[i]/3600.))


break_flag = False
plan = []

choice = input('Select target: ')
start = dusk_jd
header = ' # | UTC   | Star_ID                  | Exp. (s)'
fmt = '{:>2.0f} | {:02.0f}:{:02.0f} | {:<24} |{:>8.0f}'
plan_text = [header]

while type(choice) is str:
	if choice.lower() not in 'ucpn':
		try: choice = int(choice)
		except: break

		begin_exp = np.where(times<=start.value)[0][-1]
		start += TimeDelta(exptimes[choice]+300.0, format='sec').jd # Include a 5-minute delay
		exptime = exptimes[choice]
		try: end_exp = np.where(times>=start.value)[0][0]
		except:
			extra_time = (start-dawn_jd).sec
			start -= TimeDelta(extra_time, format='sec')
			print('WARNING: {:} requires {:.0f}s more exposure past dawn!'.format(targets[choice], extra_time))
			end_exp = len(times)-1
			exptime -= extra_time
			break_flag = True
	
		ax.plot(decimal_times[begin_exp:end_exp], airmass[choice][begin_exp:end_exp], label=targets[choice], color='C%i' %(len(plan)))
		ax.legend(loc='upper left', bbox_to_anchor=(1.1,1), fontsize='small', handlelength=1.0, frameon=False)
	
		plan.append([choice, int(decimal_times[begin_exp]), decimal_times[begin_exp]%1 * 60,
					 targets[choice], exptime])
		plan_text.append(fmt.format(*plan[-1]))

	elif choice in 'uU':
		start -= TimeDelta(plan[-1][-1]+300.0, format='sec')
		plan = plan[:-1]
		plan_text = plan_text[:-1]
		ax.lines.pop(-1)

	elif choice in 'cC':
		for p in plan: ax.lines.pop(-1)
		start = dusk_jd
		plan = []
		plan_text = [header]
		
	elif choice in 'pP':
		print('\n'.join(plan_text))
	
	elif choice in 'nN':
		plan_text.append('-'*len(header))
		start = dusk_jd
	
	else: break
	
	plt.draw()
	plt.pause(0.01)
	if break_flag:
		plan_text[-1] += '+'
		break_flag = False
		print('- End of night -')
		choice = input('Select target or option (u=undo, c=clear, p=show plan, n=next night): ')
		if choice.lower() not in 'ucpn':
			plan_text.append('-'*len(header))
			start = dusk_jd
		continue

	choice = input('Select target or option (u=undo, c=clear, p=show plan, n=next night): ')

print('\n'.join(plan_text))
plt.close()





	